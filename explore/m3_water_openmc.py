# -*- coding: utf-8 -*-
# ============================================================
#  m3_water_openmc.py  --  M3 对拍：OpenMC 水模体，开 S(α,β)（c_H_in_H2O）；弹性真实角分布 OpenMC 本来就开着
#  Claude 2026-11-06 写（跑的脚本，不默写）。材料用你的 mater_water(sab=True)。
#
#  只在 WSL 里跑（不要 F5 / Ctrl+F5）：
#    先输 wsl → conda activate /home/siyiwu/miniconda3/envs/openmc → cd /mnt/e/BNCT_GPU/my_code
#    → python explore/m3_water_openmc.py            （默认只跑种子 1，约 2–3 分钟）
#  和 openmc_water.py 主程序一模一样（同一张 tally 清单、100 批 × 10 万、同一个 run()），只差两处：
#    ① 材料 mater_water(sab=True)  ② 文件夹 / json 名字带 sab（不盖掉 10/24 的 openmc_water_s*.json）
#  存：outputs/openmc/water_sab_s1/（OpenMC 输入输出）· outputs/openmc_water_sab_s1.json（check_M3_3.py 读）
#  要补量 ρ̄ 才带种子号：python explore/m3_water_openmc.py 2 3 4 5 6 7 8 9 10（check_M3_3.py 文件头说什么时候要）
# ============================================================

# ---- 路径（2026-11-10 整理文件夹，Claude 加）：本文件在 my_code 的子文件夹里。
#      下面几行把 my_code 根目录和各子文件夹加进搜索路径，后面的 import 才找得到别的文件。----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import os
import sys
import json

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # my_code 根目录：本文件在 my_code/explore/ 里

import openmc
from openmc_tissue import kur, run, tallies, path1
from openmc_water import mater_water, water_tally, water_total_tally, WATER_TALLIES, WATER_TOTAL_TALLIES
from dose import OUT_DIR

if __name__ == "__main__":
    seeds = [int(a) for a in sys.argv[1:]] or [1]
    mat = mater_water(sab=True)

    # [Claude assert] 抓：sab 参数没接上 / 名字拼错 —— 材料里没有热散射表，OpenMC 就照旧把水里的 H 当自由气体，对拍白跑
    #   看的是 OpenMC 真正拿到的输入（写出来的 materials.xml），不是 Python 里的变量
    chk = os.path.join(path1, "water_sab_check")
    os.makedirs(chk, exist_ok=True)
    openmc.Materials([mat]).export_to_xml(os.path.join(chk, "materials.xml"))
    with open(os.path.join(chk, "materials.xml"), "r") as f:
        xml = f.read()
    assert 'name="c_H_in_H2O"' in xml, "materials.xml 里没有 c_H_in_H2O 热散射表 —— mater_water(sab=True) 没加上 S(α,β)"
    print("材料过：materials.xml 里有 c_H_in_H2O")

    for seed in seeds:
        tal = water_tally(WATER_TALLIES)
        for t in water_total_tally(WATER_TOTAL_TALLIES):     # 两个集合合成一个：run() 只收一个
            tal.append(t)
        inp = os.path.join(path1, "water_sab_s%d" % seed)
        stat_path = run(kur(), tal, inp, batch=100, n_bat=100000, seed=seed, mat=mat)
        out = {}
        for a, b, c in WATER_TALLIES + WATER_TOTAL_TALLIES:
            m, s = tallies(stat_path, a)
            out[a] = {"mean": m, "se": s}
        out["batches"] = 100
        out["particles"] = 100000
        out["seed"] = seed
        out["sab"] = "c_H_in_H2O"
        path = os.path.join(OUT_DIR, "openmc_water_sab_s%d.json" % seed)
        with open(path, "w") as f:
            json.dump(out, f)
        print("种子 %d 存好：%s" % (seed, path))
