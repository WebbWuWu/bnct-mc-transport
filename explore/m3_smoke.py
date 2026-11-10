# -*- coding: utf-8 -*-
# ============================================================
#  m3_smoke.py  --  M3 小跑通：水模体，S(α,β) / 真实角分布 开、关各跑一次（小 N）
#  Claude 2026-11-05 写（跑的脚本，不默写）
#
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5（只用 h5py + numpy）。
#  三种跑法，同一个种子、同样的输入（和 water.py 交付那次一样：KUR 源、z_src = 0、隐式、15 cm、30 层）：
#    M2   ：thermal = None、angle = None（= 老程序）
#    +S   ：只开 S(α,β)
#    +S+A ：S(α,β) + 弹性真实角分布（H1、O16）
#  N 小（每种 3000 条），只看两件事：
#    ① 跑得通、开关真的接上了（开了以后结果不再和 M2 逐位相同 —— 不同才对）
#    ② 往哪个方向变、大概多少（涨落大，R 列照实看；定量结论等 11/6 水模体对 OpenMC）
#  验证（2 条）：
#    [M3-S1] +S 和 M2 不逐位相同      抓：thermal 没传到 lookup_xs / collide（开关没接上）
#    [M3-S2] +S+A 和 +S 不逐位相同    抓：angle 没传到 collide、或 mu_cm 没传给 free_gas_kinematics
# ============================================================

import os
import sys
import time
import math

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # my_code 根目录：本文件在 my_code\explore\ 里

from transport import run_slab
from water import build_water
from dose import kur_table, KEEP_FACTOR, D, NL
from thermal import load_thermal, SAB_PATH, T_TAG
from elastic_angle import load_elastic_angle, LIB, NUCS

N_SMOKE = 3000
SEED_SMOKE = 20261105
LAYERS = (0, 3, 6, 10, 16, 22, 29)                  # 要打印的层（每层 0.5 cm）
QTY = (("track", "通量"), ("Br", "D_B 反应"), ("DH", "D_H 反冲"))


def mean_R(r, key, k, N):
    """第 k 层的每条历史均值和相对误差 R = SE / 均值。"""
    s = r[key + "_sum"][k]
    q = r[key + "_sq"][k]
    m = s / N
    var = max(q / N - m * m, 0.0)
    se = math.sqrt(var / (N - 1))
    return m, (se / m if m != 0.0 else float("inf"))


def run_one(name, nuclides, order, edges, shares, thermal, angle):
    t0 = time.time()
    r = run_slab(N_SMOKE, SEED_SMOKE, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                 src=(edges, shares), thermal=thermal, angle=angle)
    print("%-5s 跑完：%.1f s · 每条历史碰撞 %.2f 次 · 透射 %.4f · 反射 %.4f"
          % (name, time.time() - t0, r["coll_per_hist"], r["frac_T"], r["frac_R"]))
    return r


if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    t0 = time.time()
    thermal = load_thermal(SAB_PATH, T_TAG)
    angle = load_elastic_angle(LIB, NUCS)
    print("S(α,β) + 角分布读好（%.1f s）" % (time.time() - t0))

    res = {}
    res["M2"] = run_one("M2", nuclides, order, edges, shares, None, None)
    res["+S"] = run_one("+S", nuclides, order, edges, shares, thermal, None)
    res["+S+A"] = run_one("+S+A", nuclides, order, edges, shares, thermal, angle)

    # [M3-S1] 开了 S(α,β) 却和 M2 逐位相同 → thermal 没接上
    assert res["+S"]["track_sum"] != res["M2"]["track_sum"], "[M3-S1] 开了 S(α,β)，通量和 M2 逐位相同 —— thermal 没传下去"
    # [M3-S2] 再开角分布却和只开 S(α,β) 逐位相同 → angle / mu_cm 没接上
    assert res["+S+A"]["track_sum"] != res["+S"]["track_sum"], "[M3-S2] 开了角分布，通量和只开 S(α,β) 逐位相同 —— angle 或 mu_cm 没传下去"
    print("M3-S1 ✅ M3-S2 ✅ 两个开关都接上了（结果都变了）\n")

    for key, label in QTY:
        print("【%s】每条历史均值（括号里是 R）· 比值 = 开 / M2" % label)
        print("  层   深度cm   %-20s %-20s %-20s  +S/M2   +S+A/M2" % ("M2", "+S", "+S+A"))
        for k in LAYERS:
            row = [mean_R(res[c], key, k, N_SMOKE) for c in ("M2", "+S", "+S+A")]
            m0 = row[0][0]
            print("  %2d   %5.2f    %s  %6.3f  %6.3f" % (
                k, (k + 0.5) * D / NL,
                "  ".join("%11.4e (%.2f)" % (m, R) for m, R in row),
                row[1][0] / m0 if m0 else float("nan"), row[2][0] / m0 if m0 else float("nan")))
        print()
    print("R > 0.2 的格子不可信（规则 14），比值别拿去下结论。")
