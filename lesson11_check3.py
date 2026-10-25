# -*- coding: utf-8 -*-
# ============================================================
#  lesson11_check3.py  --  M2 判据 3（方案丙）：你的水模体 vs OpenMC 水模体，逐层 + 总量
#  判据归 Claude，2026-10-25 写。门槛 10/24 在看对比数据之前定死（判据3门槛重标方案 §6）。
#  被测：outputs\water_N100000_s1010.json（你的 run_slab，10 万条） vs  outputs\openmc_water_s1.json（OpenMC 种子 1，1000 万）
#  运行：Windows Ctrl+F5（默认比 OpenMC 种子 1）；或在终端里带文件名比别的 OpenMC 结果，例如方案乙：
#        python lesson11_check3.py openmc_water_iso_s1.json
#        只读 json，不 import openmc，几秒钟。
# ============================================================
#
#  比的是同一个量（纪律 2；对应表见 OpenMC水模体用法_2026-10-24.md 第 0 节）—— 都是「每个源中子」
#    你的键                              OpenMC tally（× 系数）
#    track_sum / track_sq               water flux                    cm
#    kerma_sum[H1 / O16 / B10]          water H1 / O16 / B10 heating  eV
#    gH_sum / gH_sq                     water H1 n,gamma              个
#    Bc_sum / Bc_sq                     water B10 n,a × Q_val         eV（Q_val = lesson07_freegas.py 第 188 行，2341569 eV）
#    track_tot_* / kerma_tot_*[核素]    water total flux / heating    （总量，只有这 4 个）
#
#  记号
#    m, se     均值、均值标准误；你的用 mcstat.mean_se(Σx, Σx², N) 算，OpenMC 的是它自报的
#    R         相对误差 = se / m。R > 0.2 的层 SE 不可信（规则 14），不进形状检验
#    z_k       第 k 层：(m你 − mOpenMC) / sqrt(se你² + seOpenMC²)
#    z̄         留下来的层 z 的平均 —— 公共模式涨落（30 层一起偏）大部分在这里
#    ρ̄         层间相关（10/24 用 OpenMC 10 个种子量的，往小取到 0.1）
#
#  判据（4 条；每条一行「抓什么」）
#    C1 前置：OpenMC 是 100 批 × 10 万；两边各自的「总量」= 自己 30 层之和（相对差 ≤ 1e-9）
#       抓：读错文件 / 键配错核素 / 方案乙的 json 少存了总量 tally
#    C2 总量：4 个总量 |z| ≤ 3
#       抓：整体偏（30 层一起高或一起低）—— 比如截面、密度、源归一化差了一截
#    C3 形状：Σ(z_k − z̄)² / (1 − ρ̄) 落在 χ²(留下的层数 − 1) 的 [0.5%, 99.5%] 分位之间（30 层全留 = [13.12, 52.34]）
#       抓：曲线形状不一样（比如深处偏低、浅处偏高）—— D_H 深处偏低就是这一类
#    C4 试响（纪律 4）：在 H1 heating 上故意埋两种错，C2、C3 必须各响一次
#       ① 你的每层均值加一个「一头 −6、另一头 +6 个合并 SE」的线性倾斜 → C3 必须不过
#       ② 你的总量加 6 个合并 SE → C2 必须不过
#       抓：判据写成了怎么都过
#
#  C2、C3 不用 assert 停下来：不过的量照样把 30 层 z 打出来，你要拿去查偏差。C1、C4 是 assert。
# ============================================================
import json
import math
import os
import sys

from mcstat import mean_se

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")
NL = 30
DZ = 0.5
MY_JSON = "water_N100000_s1010.json"
N_HIST = 100000
OMC_DEFAULT = "openmc_water_s1.json"
Q_VAL = 2341569.0
R_MAX = 0.2
Z_TOT = 3.0

# (显示名, 你的 sum 键, 你的 sq 键, 子键, OpenMC 30 层 tally, OpenMC 总量 tally, 你的总量 sum 键, 你的总量 sq 键, OpenMC 系数, ρ̄)
TABLE = [
    ("通量",        "track_sum", "track_sq", None,  "water flux",        "water total flux",        "track_tot_sum", "track_tot_sq", 1.0,   0.2),
    ("H1 heating",  "kerma_sum", "kerma_sq", "H1",  "water H1 heating",  "water total H1 heating",  "kerma_tot_sum", "kerma_tot_sq", 1.0,   0.3),
    ("O16 heating", "kerma_sum", "kerma_sq", "O16", "water O16 heating", "water total O16 heating", "kerma_tot_sum", "kerma_tot_sq", 1.0,   0.3),
    ("B10 heating", "kerma_sum", "kerma_sq", "B10", "water B10 heating", "water total B10 heating", "kerma_tot_sum", "kerma_tot_sq", 1.0,   0.2),
    ("H1 (n,γ)",    "gH_sum",    "gH_sq",    None,  "water H1 n,gamma",  None,                      None,            None,           1.0,   0.2),
    ("B10 (n,α)",   "Bc_sum",    "Bc_sq",    None,  "water B10 n,a",     None,                      None,            None,           Q_VAL, 0.2),
]

# χ² 分布的 0.5% / 99.5% 分位（自由度 = 留下的层数 − 1；scipy.stats.chi2.ppf 算的，Claude 容器里）
CHI2_995 = {15: (4.60, 32.80), 16: (5.14, 34.27), 17: (5.70, 35.72), 18: (6.26, 37.16), 19: (6.84, 38.58),
            20: (7.43, 40.00), 21: (8.03, 41.40), 22: (8.64, 42.80), 23: (9.26, 44.18), 24: (9.89, 45.56),
            25: (10.52, 46.93), 26: (11.16, 48.29), 27: (11.81, 49.64), 28: (12.46, 50.99), 29: (13.12, 52.34)}


def pick(d, key, sub):
    v = d[key]
    return v[sub] if sub is not None else v


def mine_layers(d, sk, qk, sub):
    su, sq = pick(d, sk, sub), pick(d, qk, sub)
    m, s = [], []
    for k in range(NL):
        mk, _, sk_ = mean_se(su[k], sq[k], N_HIST)
        m.append(mk)
        s.append(sk_)
    return m, s


def shape_test(m1, s1, m2, s2, rho):
    z = [None] * NL
    keep = []
    for k in range(NL):
        if m1[k] > 0 and m2[k] > 0 and s1[k] / m1[k] <= R_MAX and s2[k] / m2[k] <= R_MAX:
            z[k] = (m1[k] - m2[k]) / math.sqrt(s1[k] ** 2 + s2[k] ** 2)
            keep.append(k)
    zbar = sum(z[k] for k in keep) / len(keep)
    Q = sum((z[k] - zbar) ** 2 for k in keep) / (1 - rho)
    lo, hi = CHI2_995[len(keep) - 1]
    return {"z": z, "keep": keep, "zbar": zbar, "Q": Q, "lo": lo, "hi": hi, "ok": lo <= Q <= hi}


def total_z(m1, s1, m2, s2):
    return (m1 - m2) / math.sqrt(s1 ** 2 + s2 ** 2)


if __name__ == "__main__":
    omc_name = sys.argv[1] if len(sys.argv) > 1 else OMC_DEFAULT
    with open(os.path.join(OUT_DIR, MY_JSON), "r") as f:
        mine = json.load(f)
    with open(os.path.join(OUT_DIR, omc_name), "r") as f:
        omc = json.load(f)
    print("你的：%s（N = %d）   vs   OpenMC：%s\n" % (MY_JSON, N_HIST, omc_name))

    # ---- C1 前置 ----
    assert omc["batches"] == 100 and omc["particles"] == 100000, \
        "C1：%s 是 %r 批 × %r，应该是 100 × 10 万" % (omc_name, omc["batches"], omc["particles"])
    for name, sk, qk, sub, t30, ttot, mtk, mtq, coe, rho in TABLE:
        if ttot is None:
            continue
        a, b = omc[ttot]["mean"][0], sum(omc[t30]["mean"])
        assert abs(a - b) <= 1e-9 * abs(b), "C1：OpenMC 的 %s = %.12g，30 层之和 %.12g" % (ttot, a, b)
        a, b = pick(mine, mtk, sub), sum(pick(mine, sk, sub))
        assert abs(a - b) <= 1e-9 * abs(b), "C1：你的 %s[%s] = %.12g，30 层之和 %.12g" % (mtk, sub, a, b)
    print("C1 过：OpenMC 100 批 × 10 万；两边 4 个总量 = 各自 30 层之和\n")

    # ---- C2 + C3 ----
    report = {"mine": MY_JSON, "openmc": omc_name, "items": {}}
    print("%-12s %9s %5s   %9s %4s %16s %6s   %s" % ("量", "总量 z", "C2", "形状 Q", "层数", "门槛", "C3", "z̄（30 层一起偏多少）"))
    all_ok = True
    for name, sk, qk, sub, t30, ttot, mtk, mtq, coe, rho in TABLE:
        m1, s1 = mine_layers(mine, sk, qk, sub)
        m2 = [x * coe for x in omc[t30]["mean"]]
        s2 = [x * coe for x in omc[t30]["se"]]
        sh = shape_test(m1, s1, m2, s2, rho)
        item = {"z": sh["z"], "ratio": [m1[k] / m2[k] if m2[k] > 0 else None for k in range(NL)],
                "R_mine": [s1[k] / m1[k] if m1[k] > 0 else None for k in range(NL)],
                "Q": sh["Q"], "keep": sh["keep"], "zbar": sh["zbar"], "rho": rho, "C3": sh["ok"]}
        if ttot is not None:
            tm1, _, ts1 = mean_se(pick(mine, mtk, sub), pick(mine, mtq, sub), N_HIST)
            zt = total_z(tm1, ts1, omc[ttot]["mean"][0], omc[ttot]["se"][0])
            ok2 = abs(zt) <= Z_TOT
            item.update({"z_tot": zt, "C2": ok2, "ratio_tot": tm1 / omc[ttot]["mean"][0]})
            c2 = "%9.2f %5s" % (zt, "过" if ok2 else "✗")
        else:
            ok2 = True
            c2 = "%9s %5s" % ("—", "—")
        all_ok = all_ok and ok2 and sh["ok"]
        print("%-12s %s   %9.2f %4d   [%6.2f, %6.2f] %6s   %+.2f"
              % (name, c2, sh["Q"], len(sh["keep"]), sh["lo"], sh["hi"], "过" if sh["ok"] else "✗", sh["zbar"]))
        report["items"][name] = item

    # ---- 逐层：z 和 你 / OpenMC（查偏差用）----
    print("\n逐层（第 k 层 = 深度 [k·0.5, k·0.5+0.5) cm；「--」= R > 0.2 没进形状检验）")
    for name in report["items"]:
        it = report["items"][name]
        print("\n%s   ρ̄ = %.1f" % (name, it["rho"]))
        for row0 in range(0, NL, 10):
            zs = " ".join("%6s" % ("--" if it["z"][k] is None else "%+.2f" % it["z"][k]) for k in range(row0, row0 + 10))
            rs = " ".join("%6s" % ("--" if it["ratio"][k] is None else "%.3f" % it["ratio"][k]) for k in range(row0, row0 + 10))
            print("  %4.1f–%4.1f cm  z     %s" % (row0 * DZ, (row0 + 10) * DZ, zs))
            print("  %12s  你/OMC %s" % ("", rs))

    # ---- C4 试响（H1 heating）----
    name, sk, qk, sub, t30, ttot, mtk, mtq, coe, rho = TABLE[1]
    m1, s1 = mine_layers(mine, sk, qk, sub)
    m2, s2 = omc[t30]["mean"], omc[t30]["se"]
    tilt = [m1[k] + (-6 + 12 * k / (NL - 1)) * math.sqrt(s1[k] ** 2 + s2[k] ** 2) for k in range(NL)]
    assert not shape_test(tilt, s1, m2, s2, rho)["ok"], "C4①：H1 heating 埋了 ±6 SE 的倾斜，C3 居然过了"
    tm1, _, ts1 = mean_se(pick(mine, mtk, sub), pick(mine, mtq, sub), N_HIST)
    tm1_bad = tm1 + 6 * math.sqrt(ts1 ** 2 + omc[ttot]["se"][0] ** 2)
    assert abs(total_z(tm1_bad, ts1, omc[ttot]["mean"][0], omc[ttot]["se"][0])) > Z_TOT, "C4②：总量埋了 +6 SE，C2 居然过了"
    print("\nC4 试响过：埋倾斜 → C3 响；埋整体偏 → C2 响")

    # ---- 存 + 总结论 ----
    stem = os.path.splitext(omc_name)[0]
    path = os.path.join(OUT_DIR, "check3_" + stem + ".json")
    with open(path, "w") as f:
        json.dump(report, f)
    print("\n逐层结果存好：%s" % path)
    print("判据 3（丙）：%s" % ("全过" if all_ok else "有不过的量（✗）—— 看上面逐层 z，去查偏差"))
