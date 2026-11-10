# -*- coding: utf-8 -*-
# ============================================================
#  m3_tissue_angle_runs.py  --  M3 收尾：软组织开真实角分布（H1、O16）跑一次
#  Claude 2026-11-08 写（跑的脚本，不默写）。照 explore\m1_tissue_runs.py，只改两处：
#    ① angle = load_elastic_angle(LIB, NUCS)（H1、O16 真实角分布；C12、N14 不在 angle 里 →
#       collide 第 704 行 `nam in angle` 为 False → 照旧质心系各向同性，不会 KeyError）
#    ② 存新文件名 tissue_kermaA_N100000_s1010.json、..._s1200.json（不盖掉 11/6 的 tissue_kerma_*）
#  其余和 11/6 那次一模一样：ICRU-44 软组织 + 15 ppm B10、KUR 源、z_src = 0、隐式俘获、15 cm、30 层、
#  10 万条、新发生器、种子 1010 / 1200、每个核素加 MT=301 heat；S(α,β) 关（M1 的 OpenMC 没开 S）。
#  要回答的：11/6 看到的「浅处 0–3 cm H1 / C12 / O16 kerma 高 4–8%」开了角分布还在不在；
#           C12（N14 OpenMC 没打分，看不到）要不要也接角分布。
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5（不是 F5）。每个种子约 13 分钟（11/8 实测 12.8′）。
#  11/8 补种子（用户定甲）：再跑 2001–2008 共 8 个，约 1 h 45；已有的 json 自动跳过。跑完 Ctrl+F5 跑 check_tissue_angle10.py。
#  跑完：Ctrl+F5 跑 my_code\check_tissue_angle.py。
# ============================================================

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # my_code 根目录：本文件在 my_code\explore\ 里

from transport import run_slab
from dose import build_tissue, kur_table, KEEP_FACTOR, D, NL, N, T, OUT_DIR, save_result
from heating import load_heating
from elastic_angle import load_elastic_angle, LIB, NUCS

SEEDS = (1010, 1200) + tuple(range(2001, 2009))     # 1010、1200 = 和 11/6 那次（角分布关）同种子（已跑完，自动跳过）；
                                                    # 2001–2008 = 11/8 用户定甲补的 8 个 → 凑 10 个，用种子散布当误差（check_tissue_angle10.py）

if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_tissue(KEEP_FACTOR * edges[-1])
    for nuc in nuclides:
        E2, heat2 = load_heating(nuc["name"], T)
        # [Claude assert] 抓：heat 和截面不是同一张能量网格（温度标签 / 核素读错）—— 同 m1_tissue_runs 那条
        assert len(E2) == len(nuc["E"]) and np.array_equal(E2, nuc["E"]), \
            "%s：load_heating 的能量网格和 load_nuclide 的不是同一张" % nuc["name"]
        nuc["heat"] = heat2
    print("核素：" + "、".join("%s（dens %.4e）" % (n["name"], n["dens"]) for n in nuclides))

    angle = load_elastic_angle(LIB, NUCS)
    names = [n["name"] for n in nuclides]
    # [Claude assert] 抓：NUCS 里的核素不在软组织里（名字对不上 → 角分布悄悄没用上）
    assert all(nm in names for nm in angle), "angle 里有软组织没有的核素：%s" % list(angle)
    print("角分布：" + "、".join(angle) + " 真实；" + "、".join(nm for nm in names if nm not in angle) + " 各向同性")

    for seed in SEEDS:
        fname = "tissue_kermaA_N%d_s%d" % (N, seed)
        stem = os.path.join(OUT_DIR, fname)
        if os.path.exists(stem + ".json"):
            print("种子 %d：已经有了 %s.json，跳过" % (seed, fname))
            continue
        print("种子 %d：开始跑 %d 条……" % (seed, N))
        t0 = time.time()
        r = run_slab(N, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                     src=(edges, shares), thermal=None, angle=angle)
        print("  跑完：%.1f 分钟 · 每条历史碰撞 %.2f 次 · 透射 %.4f · 反射 %.4f"
              % ((time.time() - t0) / 60, r["coll_per_hist"], r["frac_T"], r["frac_R"]))
        save_result(r, stem)
        print("  存好：%s.json" % fname)
    print("两个种子都好了 → Ctrl+F5 跑 my_code\\check_tissue_angle.py")
