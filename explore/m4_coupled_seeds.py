# -*- coding: utf-8 -*-
# ============================================================
#  m4_coupled_seeds.py  --  M4 耦合对拍（P7）你这边的数：run_slab(gamma=True) 跑 20 个种子 × 2 万条
#  Claude 2026-11-25 写（跑的脚本，不默写）。由 m3_sa_seeds.py 改三处：gamma=True · 种子 4001–4020 · 结果里多存一个 "N"
#
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5。
#  用时：M3 同样的设置每个种子约 3 分钟；开光子后每条历史平均多跑约 2 个光子，预计每个种子 3–4 分钟，
#        20 个约 1–1.5 小时（Claude 没在真库上量过，每跑完一个种子打印一次用时）。可以挂着去做别的。
#  输入和 M3 一样：水模体、KUR 源、z_src = 0、隐式俘获、15 cm、30 层，thermal + angle 都开；只是多开 gamma。
#  存：outputs\coupled_SA_N20000_s4001.json … s4020.json。已经存在的跳过（中途断了重开会接着跑）。
#  跑完 + WSL 里 openmc_coupled.py 也跑完 → Ctrl+F5 跑 my_code\checks\check_p7.py。
# ============================================================

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # my_code 根目录：本文件在 my_code\explore\ 里

from transport import run_slab
from water import build_water
from dose import kur_table, KEEP_FACTOR, D, NL, OUT_DIR, save_result
from thermal import load_thermal, SAB_PATH, T_TAG
from elastic_angle import load_elastic_angle, LIB, NUCS

N_SEED = 20000
SEEDS = range(4001, 4021)   # 和 M2（2001–2030）、M3（3001–3030）错开

if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    thermal = load_thermal(SAB_PATH, T_TAG)
    angle = load_elastic_angle(LIB, NUCS)
    t_all = time.time()
    for seed in SEEDS:
        stem = os.path.join(OUT_DIR, "coupled_SA_N%d_s%d" % (N_SEED, seed))
        if os.path.exists(stem + ".json"):
            print("种子 %d：已经有了，跳过" % seed)
            continue
        t0 = time.time()
        r = run_slab(N_SEED, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                     src=(edges, shares), thermal=thermal, angle=angle, gamma=True)
        r["N"] = N_SEED          # check_p7 要用 N 把 _sum / _sq 换成每个源中子的均值和 SE
        save_result(r, stem)
        print("种子 %d 存好（%.1f 分钟，累计 %.1f 分钟）" % (seed, (time.time() - t0) / 60, (time.time() - t_all) / 60), flush=True)
    print("种子都好了 → WSL 里 openmc_coupled.py 也跑完后，Ctrl+F5 跑 my_code\\checks\\check_p7.py")
