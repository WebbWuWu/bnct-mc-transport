# -*- coding: utf-8 -*-
# ============================================================
#  m3_sa_seeds.py  --  M3：S+A 开着跑 10 个种子 × 2 万条，用你自己的种子标定误差和门槛（11/6 用户定甲）
#  Claude 2026-11-06 写（跑的脚本，不默写）
#
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5。每个种子约 3 分钟，10 个约 30 分钟。
#  输入和 m3_water_mine.py 一样（水模体、KUR 源、z_src = 0、隐式俘获、15 cm、30 层，thermal + angle 都开），
#  只是 N = 2 万、种子 3001–3010（和 10/31 那 30 个 M2 种子 2001–2030 错开）。
#  存：outputs\water_SA_N20000_s3001.json … s3010.json。已经存在的跳过（中途断了重开会接着跑）。
#  跑完 Ctrl+F5 跑 my_code\check3_m3_seeds.py。
# ============================================================

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # my_code 根目录：本文件在 my_code\explore\ 里

from lesson07_freegas import run_slab
from lesson11_water import build_water
from lesson08_dose import kur_table, KEEP_FACTOR, D, NL, OUT_DIR, save_result
from thermal import load_thermal, SAB_PATH, T_TAG
from elastic_angle import load_elastic_angle, LIB, NUCS

N_SEED = 20000
SEEDS = range(3001, 3011)

if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    thermal = load_thermal(SAB_PATH, T_TAG)
    angle = load_elastic_angle(LIB, NUCS)
    t_all = time.time()
    for seed in SEEDS:
        fname = "water_SA_N%d_s%d" % (N_SEED, seed)
        stem = os.path.join(OUT_DIR, fname)
        if os.path.exists(stem + ".json"):
            print("种子 %d：已经有了，跳过" % seed)
            continue
        t0 = time.time()
        r = run_slab(N_SEED, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                     src=(edges, shares), thermal=thermal, angle=angle)
        save_result(r, stem)
        print("种子 %d 存好（%.1f 分钟，累计 %.1f 分钟）" % (seed, (time.time() - t0) / 60, (time.time() - t_all) / 60))
    print("10 个种子都好了 → Ctrl+F5 跑 my_code\\check3_m3_seeds.py")
