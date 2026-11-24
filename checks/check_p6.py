# -*- coding: utf-8 -*-
# ============================================================
#  check_p6.py  --  M4 判据 P6：光子单独对拍。2.2233 MeV 光子束垂直射进 15 cm 水，逐层沉积能量 vs OpenMC
#  判据归 Claude。2026-11-24 写。运行：Windows Ctrl+F5（openmc (3.10.20) 环境），要读光子库。
#  约 2–5 分钟（Claude 这边没有真库，没量过速度；每跑完一个种子打印一次用时）。
#  ⚠️ 写的时候 Claude 那边没有 h5py 和真库，只在假数据上跑通过。第一次跑的完整输出发给 Claude；
#     它自己报错（不是判据响）多半是 Claude 的锅。
#
#  前提：已经在 WSL 跑过 openmc_compare/openmc_photon.py → outputs\openmc_photon_led.json、_ttb.json（11/24 16:24 跑过）
#
#  比什么（两条独立的路；先确认比的是同一个量）：
#    路 A（你的程序）  photon.run_beam：Dg_t（径迹长度 × μ_tr × E）、Dg_c（碰撞处记交出的能量），每层「每个源光子沉积多少 eV」
#    路 B（OpenMC）   heating tally，同一份光子库、同样 1 keV 截断、电子 LED（= 你的近似 #40），每层「每个源光子 eV」
#    两边几何一样：15 cm、30 层、x / y 无限大（OpenMC 用反射面）、z 两头真空；源一样：z = 0 起、沿 +z、单能 2.2233 MeV。
#    剩下唯一的差别 = 近似 #32：你的康普顿用自由电子 KN；OpenMC 用束缚电子（S(q) 角分布 + 多普勒）。
#    总截面两边都是库里那一列（含束缚效应），所以 #32 只影响「往哪散、交出多少」，预计很小。
#
#  门槛（2026-11-24 20:30 写进文件；还没看过你程序跑出的任何数。OpenMC 两个 json 只看过总量和 R）：
#    一个量「过」= 下面两条满足一条：
#      ① |z| ≤ 3，z = (你 − OpenMC) / √(SE_你² + SE_OpenMC²)        —— 统计上分辨不出差别
#      ② |相对差| < 容差：总量 0.5%、逐层 1%                         —— 分辨出了，但在近似 #32 预计的范围里（照实打印大小）
#    照实交代：定容差前 Claude 用「假截面」（只有 KN、没有相干、光电随手编）把本文件试跑通过，
#    看到 Dg_t 是径迹估计、方差很小，30 万条时逐层分辨力约 0.3–0.4%，比 #32 可能的大小还细。
#    #32 改的是散射方向和交出比例，会让能量在层与层之间挪一点，逐层比总量敏感 → 逐层容差单独放到 1%。
#    [P6-1] 整块板总量（Dg_t、Dg_c 各判一次）：必须过。
#           总量是 30 层合成的一个数 → 误差用种子间散布（纪律 10）：你这边分 K = 20 个种子跑，
#           20 个总量的均值和 SE（t 分布 19 个自由度，|t| > 3 的误报约 0.7%）
#    [P6-2] 逐层（Dg_t、Dg_c 各判一次）：R = SE / 均值 > 0.2 的层不判（纪律 9）；不过的层 ≤ 1 层算过。
#           为什么允许 1 层：30 层各判一次，全对时至少一层 |z| > 3 的概率 1 − 0.9973³⁰ ≈ 8%；允许 1 层后误报约 0.3%。
#  只看、不判：你 vs OpenMC TTB；OpenMC TTB vs LED（= 近似 #40 的大小，16:24 已量：总量 −0.28%）。
#
#  抓什么：μ_tr / ⟨T⟩/E 用错（Dg_t 偏）、交出能量记错（Dg_c 偏）、截面或密度错、湮灭光子没跑、
#          出板那段没记、方向或拆层错（逐层形状偏）
#  抓不到什么（照实写）：
#    · 生光子（中子 → γ）—— 归 G1–G5 和明天的耦合对拍
#    · 两边共用的东西一起错：光子库本身、LED 近似本身
#    · 比分辨力小的偏差（每次打印 3 × SE ÷ 均值）
#  自检（纪律 4；每个 K = 10 × 6000 条，只看总量 Dg_t 那一条必须不过）：
#    S1  ⟨T⟩/E 用康普顿边 2α/(1 + 2α)（最大值当平均，11/22 预测错过；同 P4 的 W4）
#    S2  四种截面都放大 10%（「密度 / 单位错」那一类）
# ============================================================
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json
import math
import time

import photon as PH
from water import build_water
from dose import kur_table, KEEP_FACTOR, OUT_DIR

D_SLAB, NL = 15.0, 30                  # 和 OpenMC 那边一样（dose.D、dose.NL）
K_MAIN, N_PER, SEED0 = 20, 15000, 81731   # 20 个种子 × 1.5 万 = 30 万条
K_SELF, N_SELF, SEED_SELF = 10, 6000, 27719
Z_GATE = 3.0
REL_TOL_TOT = 0.005
REL_TOL_LAYER = 0.01
R_MAX = 0.2
MAX_BAD_LAYERS = 1


# ---------------- 读 OpenMC ----------------
def load_openmc(e):
    path = os.path.join(OUT_DIR, "openmc_photon_%s.json" % e)
    assert os.path.exists(path), "找不到 %s —— 先在 WSL 跑 openmc_compare/openmc_photon.py" % path
    with open(path, "r", encoding="utf-8") as f:
        r = json.load(f)
    # [Claude assert] 抓：读错文件、OpenMC 那边能量 / 层数和这边对不上
    assert r["electron"] == e, "%s 里记的电子处理是 %r，不是 %r" % (path, r["electron"], e)
    assert abs(r["E0"] / PH.E_G_H - 1.0) < 1e-12, "OpenMC 用的光子能量 %.6g ≠ 你的 E_G_H %.6g" % (r["E0"], PH.E_G_H)
    m = r["photon heating"]["mean"]
    s = r["photon heating"]["se"]
    assert len(m) == NL and len(s) == NL, "OpenMC 分了 %d 层，这边是 %d 层" % (len(m), NL)
    return {"m": m, "s": s, "tot": r["photon total heating"]["mean"][0], "tot_s": r["photon total heating"]["se"][0],
            "n": r["batches"] * r["particles"]}


# ---------------- 跑你这边：K 个种子 ----------------
def run_mine(K, n_per, seed0, pho, nuclides, verbose=True):
    """返回每层的均值和 SE（30 万条合起来），以及每个种子的总量列表（算总量的种子散布用）。"""
    per = {"t": [], "c": []}               # 每个种子的结果字典
    t0 = time.time()
    for k in range(K):
        res = PH.run_beam(n_per, seed0 + 7919 * k, PH.E_G_H, pho, nuclides, D_SLAB, NL)
        per["t"].append((res["t_mean"], res["t_se"]))
        per["c"].append((res["c_mean"], res["c_se"]))
        if verbose:
            print("    种子 %2d / %d 跑完（%d 条，累计 %.0f 秒）" % (k + 1, K, n_per, time.time() - t0), flush=True)
    out = {}
    for key in ("t", "c"):
        m = [sum(per[key][k][0][j] for k in range(K)) / K for j in range(NL)]                       # 各种子条数相同 → 直接平均
        s = [math.sqrt(sum(per[key][k][1][j] ** 2 for k in range(K))) / K for j in range(NL)]      # 独立种子：SE 平方相加
        tots = [sum(per[key][k][0]) for k in range(K)]                                              # 每个种子一个总量
        tm = sum(tots) / K
        ts = math.sqrt(sum((x - tm) ** 2 for x in tots) / (K - 1) / K)                              # 种子间散布 → 均值的 SE
        out[key] = {"m": m, "s": s, "tot": tm, "tot_s": ts}
    out["n"] = K * n_per
    return out


# ---------------- 判 ----------------
def judge(x, sx, y, sy, tol=REL_TOL_TOT):
    d = x - y
    se = math.sqrt(sx * sx + sy * sy)
    z = d / se if se > 0 else float("inf")
    rel = d / y if y != 0 else float("inf")
    if abs(z) <= Z_GATE:
        return True, z, rel, "✅"
    if abs(rel) < tol:
        return True, z, rel, "✅（分辨出了，< %.1f%%，记作 #32）" % (100 * tol)
    return False, z, rel, "❌"


def compare(mine, ref, label, verbose=True):
    """mine = {"m","s","tot","tot_s"}；ref = OpenMC。返回 (总量过没过, 逐层过没过)。"""
    ok_t, z, rel, tag = judge(mine["tot"], mine["tot_s"], ref["tot"], ref["tot_s"])
    if verbose:
        print("  [P6-1] %s 总量：你 %.5e ± %.2e · OpenMC %.5e ± %.2e · 差 %+.3f%% · z = %+.2f  %s"
              % (label, mine["tot"], mine["tot_s"], ref["tot"], ref["tot_s"], 100 * rel, z, tag))
        print("         分辨力（3 × 合成 SE ÷ OpenMC）= %.3f%%"
              % (300 * math.sqrt(mine["tot_s"] ** 2 + ref["tot_s"] ** 2) / ref["tot"]))
    bad, skipped, lines = [], [], []
    for j in range(NL):
        Rm = mine["s"][j] / mine["m"][j] if mine["m"][j] > 0 else float("inf")
        Ro = ref["s"][j] / ref["m"][j] if ref["m"][j] > 0 else float("inf")
        if Rm > R_MAX or Ro > R_MAX:
            skipped.append(j + 1)
            lines.append("     %2d   %.4e  %.4e    —      —   不判（R 你 %.2f / OpenMC %.2f）" % (j + 1, ref["m"][j], mine["m"][j], Rm, Ro))
            continue
        ok, z, rel, tag = judge(mine["m"][j], mine["s"][j], ref["m"][j], ref["s"][j], REL_TOL_LAYER)
        if not ok:
            bad.append(j + 1)
        lines.append("     %2d   %.4e  %.4e  %+6.2f%%  %+6.2f  %s" % (j + 1, ref["m"][j], mine["m"][j], 100 * rel, z, tag))
    ok_l = len(bad) <= MAX_BAD_LAYERS
    if verbose:
        print("  [P6-2] %s 逐层（层 · OpenMC · 你 · 差 · z）：" % label)
        print("\n".join(lines))
        print("         不过的层：%s（允许 ≤ %d 层）；不判的层：%s  %s"
              % (bad or "无", MAX_BAD_LAYERS, skipped or "无", "✅" if ok_l else "❌"))
    return ok_t, ok_l


def look_only(mine, ref, label):
    """只看不判：逐层相对差一行打印。"""
    rel_tot = mine["tot"] / ref["tot"] - 1.0
    print("  %s：总量差 %+.3f%%；逐层差（%%）：" % (label, 100 * rel_tot))
    print("    " + " ".join("%+.2f" % (100 * (a / b - 1.0)) for a, b in zip(mine["m"], ref["m"])))


# ---------------- 自检 ----------------
def selftest(pho, nuclides, led):
    orig_tab, orig_sig = PH.kn_TE_tab, PH.photon_sigma

    def edge_tab(photon, E):                      # S1：康普顿边当平均
        a = E / PH.ME_C2
        return 2 * a / (1 + 2 * a)

    def sig_up(photon, nuclides_, E):             # S2：四种截面都 × 1.1（tot 跟着重算）
        sig = orig_sig(photon, nuclides_, E)
        for k in ("coh", "incoh", "pe", "pair"):
            sig[k] *= 1.1
        sig["tot"] = sig["coh"] + sig["incoh"] + sig["pe"] + sig["pair"]
        return sig

    all_ok = True
    for name, attr, fake in (("S1 ⟨T⟩/E 用康普顿边", "kn_TE_tab", edge_tab),
                             ("S2 截面都 × 1.1", "photon_sigma", sig_up)):
        setattr(PH, attr, fake)
        try:
            mine = run_mine(K_SELF, N_SELF, SEED_SELF, pho, nuclides, verbose=False)
        finally:
            PH.kn_TE_tab, PH.photon_sigma = orig_tab, orig_sig
        ok_t, z, rel, _ = judge(mine["t"]["tot"], mine["t"]["tot_s"], led["tot"], led["tot_s"])
        fired = not ok_t
        all_ok = all_ok and fired
        print("  %-22s Dg_t 总量差 %+.2f%%，z = %+.1f  %s" % (name, 100 * rel, z, "✅ 响了" if fired else "❌ 没响（判据不够灵，告诉 Claude）"))
    return all_ok


if __name__ == "__main__":
    t_start = time.time()
    led = load_openmc("led")
    ttb = load_openmc("ttb")
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    pho = PH.load_photon()

    print("P6：2.2233 MeV 光子束 → 15 cm 水，逐层沉积。OpenMC LED %d 个光子；你这边 %d 个种子 × %d 条" % (led["n"], K_MAIN, N_PER))
    mine = run_mine(K_MAIN, N_PER, SEED0, pho, nuclides)
    print()
    ok = True
    for key, label in (("t", "Dg_t（径迹）"), ("c", "Dg_c（碰撞）")):
        ok_t, ok_l = compare(mine[key], led, label)
        ok = ok and ok_t and ok_l
        print()

    print("只看、不判：")
    look_only(mine["t"], ttb, "你的 Dg_t vs OpenMC TTB")
    look_only(ttb, led, "OpenMC TTB vs LED（= 近似 #40）")
    print()

    print("自检：2 个故意写错的版本必须响")
    ok_s = selftest(pho, nuclides, led)

    ok = ok and ok_s
    print("\n  → %s    用时 %.0f 秒" % ("✅ P6 全过" if ok else "❌ P6 没全过", time.time() - t_start))
    assert ok, "P6 没全过：看上面哪一行是 ❌（总量 ❌ = 整体偏；逐层 ❌ = 形状偏；自检 ❌ = 判据不够灵）"
