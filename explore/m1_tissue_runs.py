# -*- coding: utf-8 -*-
# ============================================================
#  m1_tissue_runs.py  --  追查 ① 第 5 步 + 技术报告 4.5（11/6 用户定甲：两件事合并跑）
#  Claude 2026-11-06 写（跑的脚本，不默写）
#
#  运行：Windows 的 openmc (3.10.20) 环境，打开本文件 Ctrl+F5。两个种子先后跑，共约 20–30 分钟。
#  跑什么：M1 的软组织问题（dose 的 ICRU-44 软组织 H/C/N/O + 15 ppm B10、ρ = 1.06、KUR 源、z_src = 0、
#          隐式俘获、15 cm、30 层、10 万条），和 M1 只差三处：
#    ① 新发生器（10/31 乙方案）       → 回答技术报告 4.5「D_B 深处低 2–3% 是不是随机数重叠」
#    ② 每个核素加上 MT=301「heat」     → run_slab 会按核素记 kerma（M2 的记账），回答追查 ①「换成 kerma 还差不差 20%」
#    ③ 种子和 M1 一样：1010、1200     → 新旧发生器同种子
#  S(α,β)、角分布都关（thermal = None、angle = None）：和 M1 的 OpenMC（不开 S(α,β)）一一对应。
#  核素表用你 dose.py 的 build_tissue（不改它），这里只给每个核素字典加一个 "heat" 键 ——
#  和 water.build_water 里 heat 的来源（heating.load_heating）一模一样。
#  存：outputs\tissue_kerma_N100000_s1010.json、..._s1200.json。已经存在就跳过。跑完 Ctrl+F5 跑 my_code\check_tissue.py。
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

SEEDS = (1010, 1200)                                # M1 用过的两个种子

if __name__ == "__main__":
    edges, shares = kur_table()
    nuclides, order = build_tissue(KEEP_FACTOR * edges[-1])
    for nuc in nuclides:
        E2, heat2 = load_heating(nuc["name"], T)
        # [Claude assert] 抓：heat 和截面不是同一张能量网格（温度标签 / 核素读错）—— 同 build_water 里那条
        assert len(E2) == len(nuc["E"]) and np.array_equal(E2, nuc["E"]), \
            "%s：load_heating 的能量网格和 load_nuclide 的不是同一张" % nuc["name"]
        nuc["heat"] = heat2
    print("核素：" + "、".join("%s（dens %.4e）" % (n["name"], n["dens"]) for n in nuclides))

    for seed in SEEDS:
        fname = "tissue_kerma_N%d_s%d" % (N, seed)
        stem = os.path.join(OUT_DIR, fname)
        if os.path.exists(stem + ".json"):
            print("种子 %d：已经有了 %s.json，跳过" % (seed, fname))
            continue
        print("种子 %d：开始跑 %d 条……" % (seed, N))
        t0 = time.time()
        r = run_slab(N, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                     src=(edges, shares), thermal=None, angle=None)
        print("  跑完：%.1f 分钟 · 每条历史碰撞 %.2f 次 · 透射 %.4f · 反射 %.4f"
              % ((time.time() - t0) / 60, r["coll_per_hist"], r["frac_T"], r["frac_R"]))
        save_result(r, stem)
        print("  存好：%s.json" % fname)
    print("两个种子都好了 → Ctrl+F5 跑 my_code\\check_tissue.py")
