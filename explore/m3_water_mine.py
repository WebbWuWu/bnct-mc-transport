# -*- coding: utf-8 -*-
# ============================================================
#  m3_water_mine.py  --  M3 对拍：你的水模体 10 万条，S(α,β) + 弹性真实角分布 都开
#  Claude 2026-11-06 写（跑的脚本，不默写）
#
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5（只用 h5py + numpy）。约 15–25 分钟，中间不打印，别关。
#  输入和 M2 交付那次（water.py 主程序）一模一样：KUR 源、z_src = 0、隐式俘获、15 cm、30 层、N = 10 万、种子 1010；
#  只多两个参数 thermal、angle（和 11/5 m3_smoke.py 的「+S+A」那一种一样，只是 N 从 3000 换成 10 万）。
#  存：outputs\water_SA_N100000_s1010.json / .npz（SA = S(α,β) + Angle）。已经存在就不重跑。
#  下一步：Ctrl+F5 跑 my_code\check_M3_3.py（要先在 WSL 跑完 explore/m3_water_openmc.py）。
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
    stem = os.path.join(OUT_DIR, "water_SA_N%d_s%d" % (N, SEED))
    if os.path.exists(stem + ".json"):
        print("已经有了：%s.json —— 不重跑。要重跑先把它删掉。" % stem)
        sys.exit(0)

    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    t0 = time.time()
    thermal = load_thermal(SAB_PATH, T_TAG)
    angle = load_elastic_angle(LIB, NUCS)
    print("S(α,β) + 角分布读好（%.1f s）。开始跑 %d 条（种子 %d），约 15–25 分钟……" % (time.time() - t0, N, SEED))

    t0 = time.time()
    r = run_slab(N, SEED, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                 src=(edges, shares), thermal=thermal, angle=angle)
    print("跑完：%.1f 分钟 · 每条历史碰撞 %.2f 次 · 透射 %.4f · 反射 %.4f"
          % ((time.time() - t0) / 60, r["coll_per_hist"], r["frac_T"], r["frac_R"]))
    save_result(r, stem)
    print("存好：%s.json / .npz" % stem)
