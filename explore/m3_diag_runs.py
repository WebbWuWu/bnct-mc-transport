# -*- coding: utf-8 -*-
# ============================================================
#  m3_diag_runs.py  --  M3 判据 3 C3 没过（B10 Q = 56.47）之后的两个诊断跑（11/6 用户定：甲乙都跑）
#  Claude 2026-11-06 写（跑的脚本，不默写）
#
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5。两样先后跑，共约 25 分钟，中间只在每样开头结尾打印。
#  输入和 m3_water_mine.py 一模一样（KUR 源、z_src = 0、隐式俘获、15 cm、30 层、N = 10 万），只换下面两处：
#    甲 M2new ：thermal = None、angle = None（= M2 老程序，但用的是 10/31 换过的新发生器），种子 1010
#              → 存 outputs\water_M2new_N100000_s1010.json（不盖 10/24 旧发生器的 water_N100000_s1010.json）
#    乙 SA1011：S(α,β) + 角分布都开，种子 1011（换一个种子）
#              → 存 outputs\water_SA_N100000_s1011.json
#  已经存在的就跳过。跑完 Ctrl+F5 跑 my_code\check_M3_diag.py。
# ============================================================

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # my_code 根目录：本文件在 my_code\explore\ 里

from transport import run_slab
from water import build_water
from dose import kur_table, KEEP_FACTOR, D, NL, N, SEED, OUT_DIR, save_result
from thermal import load_thermal, SAB_PATH, T_TAG
from elastic_angle import load_elastic_angle, LIB, NUCS

if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    thermal = load_thermal(SAB_PATH, T_TAG)
    angle = load_elastic_angle(LIB, NUCS)

    # (名字, 种子, thermal, angle, 存的文件名)
    CASES = [
        ("甲 M2new（S、A 都关）", SEED,     None,    None,  "water_M2new_N%d_s%d" % (N, SEED)),
        ("乙 SA1011（S、A 都开）", SEED + 1, thermal, angle, "water_SA_N%d_s%d" % (N, SEED + 1)),
    ]
    for name, seed, th, an, fname in CASES:
        stem = os.path.join(OUT_DIR, fname)
        if os.path.exists(stem + ".json"):
            print("%s：已经有了 %s.json，跳过" % (name, fname))
            continue
        print("%s：开始跑 %d 条（种子 %d）……" % (name, N, seed))
        t0 = time.time()
        r = run_slab(N, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                     src=(edges, shares), thermal=th, angle=an)
        print("  跑完：%.1f 分钟 · 每条历史碰撞 %.2f 次 · 透射 %.4f · 反射 %.4f"
              % ((time.time() - t0) / 60, r["coll_per_hist"], r["frac_T"], r["frac_R"]))
        save_result(r, stem)
        print("  存好：%s.json" % fname)
    print("两样都好了 → Ctrl+F5 跑 my_code\\check_M3_diag.py")
