# -*- coding: utf-8 -*-
# ============================================================
#  m3_A_seeds.py  --  诊断：只开角分布（A）、不开 S，跑 10 个种子 × 2 万条（11/6 18:45 用户定甲）
#  Claude 2026-11-06 写（跑的脚本，不默写）。由 m3_sa_seeds.py 改一处：thermal=None（不读 c_H_in_H2O）
#  为什么：S+A 合起来 vs OpenMC 开 S，H1 / O16 kerma 浅处低 2–3%、深处高 5–15%（形状 ✗）；
#    M2（都不开）vs OpenMC 不开 S，H1 / O16 深处反而低 5–9%。快中子不受 S 影响 → 先把 A 单独拎出来比。
#
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5（不是 F5）。每个种子约 2 分钟，10 个约 20 分钟。
#  输入和 m3_sa_seeds.py 一样（水模体、KUR 源、z_src = 0、隐式俘获、15 cm、30 层、种子 3001–3010），只是 thermal 关、angle 开。
#  存：outputs\water_A_N20000_s3001.json … s3010.json。已经存在的跳过（中途断了重开会接着跑）。
#  跑完 Ctrl+F5 跑 my_code\check_M3_seeds.py（默认就是 A 这一组）。
# ============================================================

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # my_code 根目录：本文件在 my_code\explore\ 里

from transport import run_slab
from water import build_water
from dose import kur_table, KEEP_FACTOR, D, NL, OUT_DIR, save_result
from elastic_angle import load_elastic_angle, LIB, NUCS

N_SEED = 20000
SEEDS = range(3001, 3031)   # 11/7 Claude 改（用户定甲）：3011–3030 补 20 个，已有的 3001–3010 自动跳过

if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    angle = load_elastic_angle(LIB, NUCS)
    t_all = time.time()
    for seed in SEEDS:
        fname = "water_A_N%d_s%d" % (N_SEED, seed)
        stem = os.path.join(OUT_DIR, fname)
        if os.path.exists(stem + ".json"):
            print("种子 %d：已经有了，跳过" % seed)
            continue
        t0 = time.time()
        r = run_slab(N_SEED, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                     src=(edges, shares), thermal=None, angle=angle)
        save_result(r, stem)
        print("种子 %d 存好（%.1f 分钟，累计 %.1f 分钟）" % (seed, (time.time() - t0) / 60, (time.time() - t_all) / 60))
    print("10 个种子都好了 → Ctrl+F5 跑 my_code\\check_M3_seeds.py")
