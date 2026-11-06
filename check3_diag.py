# -*- coding: utf-8 -*-
# ============================================================
#  check3_diag.py  --  M3 判据 3 C3 没过之后的诊断：三次对拍并排看
#  判据归 Claude，2026-11-06 写。运行：Windows Ctrl+F5（只读 json，几秒钟）。先跑完 explore\m3_diag_runs.py。
#
#  三对（C1–C4 和门槛都不改，直接调 lesson11_check3.run_check3）：
#    ① SA1010 ：你的 S+A 种子 1010（下午那次） vs OpenMC 开 S 种子 1      —— 已知：B10 Q = 56.47 ✗
#    ② M2new  ：你的 M2 新发生器 种子 1010     vs OpenMC 不开 S 种子 1    —— 甲：偏差 S(α,β) 之前就有吗？
#    ③ SA1011 ：你的 S+A 种子 1011            vs OpenMC 开 S 种子 1      —— 乙：1010 是不是碰巧？
#  怎么读（11/6 跑之前写下）：
#    ② 也不过、而且 B10 的 z 也是「浅处负、5–10 cm 正」 → 偏差和 S(α,β) 无关，是 M2 就有的（新发生器下第一次看到）
#    ② 过、③ 也不过且样子像 ①                        → S(α,β) / 角分布带进来的，去查热能区（M3 框架 T3 盲区：修正 E′）
#    ② 过、③ 过                                        → ① 是 0.17% 的碰巧；照实记两次，不算「判据过」，M3 小结里写明
#  ③ 和 ① 用的是同一个 OpenMC 结果：OpenMC 那一侧的涨落两次是同一份（它的 SE 是你的 1/√10，影响小），这一点写进结论。
# ============================================================
import os

from lesson11_check3 import run_check3, OUT_DIR, NL

PAIRS = [
    ("① SA1010", "water_SA_N100000_s1010.json",    "openmc_water_sab_s1.json", "check3_diag_SA1010.json"),
    ("② M2new",  "water_M2new_N100000_s1010.json", "openmc_water_s1.json",     "check3_diag_M2new.json"),
    ("③ SA1011", "water_SA_N100000_s1011.json",    "openmc_water_sab_s1.json", "check3_diag_SA1011.json"),
]
SHOW = ("通量", "H1 heating", "O16 heating", "B10 heating")

if __name__ == "__main__":
    for tag, mine, omc, out in PAIRS:
        for f in (mine, omc):
            assert os.path.exists(os.path.join(OUT_DIR, f)), "%s：找不到 outputs\\%s —— 先跑 explore\\m3_diag_runs.py" % (tag, f)

    res = {}
    for tag, mine, omc, out in PAIRS:
        print("\n" + "=" * 30 + " " + tag + " " + "=" * 30)
        res[tag] = run_check3(mine, omc, out_name=out)

    # ---- 并排总结 ----
    print("\n\n" + "#" * 70)
    print("并排总结（Q 门槛：30 层 [13.12, 52.34]，29 层 [12.46, 50.99]；总量 |z| ≤ 3）")
    print("%-12s" % "量" + "".join("%26s" % tag for tag, _, _, _ in PAIRS))
    for name in SHOW:
        row = []
        for tag, _, _, _ in PAIRS:
            it = res[tag][1]["items"][name]
            row.append("%26s" % ("Q %6.2f %s · z总 %+.2f" % (it["Q"], "过" if it["C3"] else "✗", it["z_tot"])))
        print("%-12s" % name + "".join(row))

    print("\nB10 heating 逐层 z（看样子像不像：浅处负、5–10 cm 正？）")
    for tag, _, _, _ in PAIRS:
        z = res[tag][1]["items"]["B10 heating"]["z"]
        for row0 in range(0, NL, 10):
            zs = " ".join("%6s" % ("--" if z[k] is None else "%+.2f" % z[k]) for k in range(row0, row0 + 10))
            print("  %-9s %4.1f–%4.1f cm  %s" % (tag if row0 == 0 else "", row0 * 0.5, (row0 + 10) * 0.5, zs))
    print("\n三对结论：" + " · ".join("%s %s" % (tag, "全过" if res[tag][0] else "有 ✗") for tag, _, _, _ in PAIRS))
