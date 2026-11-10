# -*- coding: utf-8 -*-
# ============================================================
#  check_M3_3.py  --  M3 判据 3（方案丙）：你的水模体（开 S(α,β) + 真实角分布）vs OpenMC 水模体（开 S(α,β)）
#  判据归 Claude，2026-11-06 写。运行：Windows Ctrl+F5（只读 json，几秒钟）。
#
#  被测：outputs\water_SA_N100000_s1010.json   ← explore\m3_water_mine.py（你的 run_slab，10 万条，种子 1010，thermal + angle 都开）
#    vs  outputs\openmc_water_sab_s1.json       ← explore\m3_water_openmc.py（OpenMC 种子 1，100 批 × 10 万，材料加 c_H_in_H2O）
#
#  C1–C4 和 M2 完全一样，直接调 check_M2_3.run_check3（同一个函数、同一张对应表、同一套门槛）：
#    C1 前置：OpenMC 100 批 × 10 万；两边 4 个总量 = 各自 30 层之和        抓：读错文件 / 键配错核素
#    C2 总量：通量、H1 / O16 / B10 heating 总量 |z| ≤ 3                    抓：整体偏（30 层一起高或低）
#    C3 形状：Σ(z − z̄)² / (1 − ρ̄) ∈ χ²(层数 − 1) 的 [0.5%, 99.5%]           抓：曲线形状不一样（比如 S(α,β) 只在一边起作用 → 浅处 / 深处反着偏）
#    C4 试响：埋 ±6 SE 倾斜 → C3 响；埋 +6 SE 整体偏 → C2 响               抓：判据写成了怎么都过
#
#  门槛 11/6 在看数据之前定（规则 14）：
#    · |z| ≤ 3、χ² 99% 区间：和 M2 一样，不改
#    · ρ̄ 沿用 10/24 用 OpenMC（无 S(α,β)）10 个种子量的、往小取过的值（通量 0.2、H1 / O16 heating 0.3、其余 0.2）
#      理由：层间相关来自「快中子源这一抽抽到高能还是低能」，30 层一起涨落；S(α,β) 只改 10 eV 以下，不碰源抽样
#      ⚠️ 如果 C3 在两个以上的量上贴着上限（Q > 45）却过了 / 刚好不过：先回头看这条 ——
#         WSL 里补跑 OpenMC 种子 2–10（python explore/m3_water_openmc.py 2 3 4 5 6 7 8 9 10，约 20 分钟）重量 ρ̄
#
#  另外打一张表（不是判据，只看）：OpenMC 自己开 / 不开 S(α,β) 的比值（同种子 1、同 1000 万）
#    = S(α,β) 在水里到底改了多少 —— 11/5 小跑通（3000 条）说「D_B 峰值附近 ↑ 14–25%、深处 ×0.5–0.6」，这里给定量
#    两边都是 OpenMC、都开真实角分布，只差 S(α,β) 一个开关 → 比的是同一个量（纪律 2）
# ============================================================
import json
import os

from check_M2_3 import run_check3, OUT_DIR, DZ, NL

MY_JSON = "water_SA_N100000_s1010.json"
OMC_SAB = "openmc_water_sab_s1.json"
OMC_NOSAB = "openmc_water_s1.json"          # 10/24 M2 那次（无 S(α,β)），只用来打下面那张比值表
LAYERS = (0, 3, 6, 10, 16, 22, 29)          # 打印的层（和 m3_smoke.py 同一组）
SHOW = (("通量", "water flux"), ("B10 (n,α)  = D_B", "water B10 n,a"), ("H1 heating = D_H", "water H1 heating"))

if __name__ == "__main__":
    for name in (MY_JSON, OMC_SAB):
        assert os.path.exists(os.path.join(OUT_DIR, name)), \
            "找不到 outputs\\%s —— 先跑 explore\\%s" % (name, "m3_water_mine.py（Windows）" if name == MY_JSON else "m3_water_openmc.py（WSL）")

    all_ok, report = run_check3(MY_JSON, OMC_SAB)

    # ---- 只看：S(α,β) 在 OpenMC 里改了多少（开 / 不开）----
    with open(os.path.join(OUT_DIR, OMC_SAB), "r") as f:
        sab = json.load(f)
    with open(os.path.join(OUT_DIR, OMC_NOSAB), "r") as f:
        nosab = json.load(f)
    print("\n—— 只看，不是判据 ——  OpenMC 开 S(α,β) / 不开（种子 1，各 1000 万；括号里是比值的相对误差，按两次独立算的 —— 同种子有正相关，真实误差更小，这是上限）")
    print("  %-18s" % "层（深度 cm）" + "".join("%13s" % ("%d（%.2f）" % (k, (k + 0.5) * DZ)) for k in LAYERS))
    for label, t in SHOW:
        cells = []
        for k in LAYERS:
            a, sa = sab[t]["mean"][k], sab[t]["se"][k]
            b, sb = nosab[t]["mean"][k], nosab[t]["se"][k]
            r = a / b
            rr = ((sa / a) ** 2 + (sb / b) ** 2) ** 0.5
            cells.append("%13s" % ("%.3f(%.3f)" % (r, rr)))
        print("  %-18s" % label + "".join(cells))
    print("\n判据 3（丙）M3：%s" % ("全过 → M3 水模体对上 OpenMC" if all_ok else "有不过的量（✗）—— 看上面逐层 z，发我"))
