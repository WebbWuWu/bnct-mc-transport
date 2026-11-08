# -*- coding: utf-8 -*-
# ============================================================
#  check_emit.py  --  M4 第 ② 段判据：emit_photons（中子碰撞时生光子、压栈）
#  判据归 Claude，2026-11-08 写（门槛在你写代码之前定）。运行：Windows Ctrl+F5（G1–G2 几秒；G4 约 1 分钟）。
#  用法清单：项目文档\emit_photons用法_2026-11-08.md
#
#  G1  假 xs：gH = 0.0411、gB = 0.0093、Σt = 0.036788 → y = 1.3700、H 份额 = 0.81548（都不是特殊值）
#      调 N = 100000 次（每次一个新的 hist），中子 p：z = 3.37、wt = 0.6173、方向随便
#      ① 每次个数 ∈ {1, 2}（y = 1.37 → ⌊y + ξ⌋ 只能是 1 或 2）
#      ② 平均个数 vs 1.37：z = (均值 − 1.37) / √(0.37·0.63/N)，|z| < 4
#      ③ H 线份额 vs 0.81548：二项 SE，|z| < 4；能量只能是 E_G_H 或 E_G_B
#      ④ ⟨w⟩ vs 0：SE = √(1/3 / M)，|z| < 4 · ⟨w²⟩ vs 1/3：SE = √(4/45 / M)，|z| < 4 · 每个方向 u² + v² + w² = 1（1e-12）
#      ⑤ 每个光子 z、wt 和中子一模一样；返回值 n = 这次压进栈的个数
#  G2  自检：5 个故意写错的版本（不取整 / 四舍五入 / 选线反了 / 方向只往前 / 权重恒 1），每个都必须被 G1 判错
#  G3  不在这里：关着 gamma 跑 check_split.py（框架判据 P5，逐位）
#  G4  只看：水模体 + KUR 源，开 gamma 跑 2000 条，打印平均每条历史生几个光子、2.2 MeV 占比（框架估约 2 个 / 历史）
# ============================================================
import math
import os
import sys

from rng import lcg, isotropic_direction
from photon import emit_photons, E_G_H, E_G_B

GH, GB, ST = 0.0411, 0.0093, 0.036788
Y = (GH + GB) / ST
FH = GH / (GH + GB)
N = 100000
Z_GATE = 4.0


def run_g1(emit, seed=48173):
    """对一个 emit 函数跑 G1，返回 (是否全过, 说明列表)。"""
    rng = lcg(seed)
    xs = {"gH_coe": GH, "gB_coe": GB, "Sigma_t": ST}
    msgs, ok = [], True
    cnt_sum, nH, M, sw, sw2 = 0, 0, 0, 0.0, 0.0
    bad_n = bad_dir = bad_copy = bad_E = bad_ret = 0
    for i in range(N):
        p = {"z": 3.37, "u": 0.28, "v": -0.41, "w": 0.8677, "E": 0.0253, "wt": 0.6173, "t": 0}
        hist = {"pstack": []}
        n = emit(p, xs, hist, rng)
        k = len(hist["pstack"])
        if n != k:
            bad_ret += 1
        if k not in (1, 2):
            bad_n += 1
        cnt_sum += k
        for g in hist["pstack"]:
            M += 1
            if g["E"] == E_G_H:
                nH += 1
            elif g["E"] != E_G_B:
                bad_E += 1
            if abs(g["u"] ** 2 + g["v"] ** 2 + g["w"] ** 2 - 1.0) > 1e-12:
                bad_dir += 1
            if g["z"] != p["z"] or g["wt"] != p["wt"]:
                bad_copy += 1
            sw += g["w"]
            sw2 += g["w"] ** 2
    mean_n = cnt_sum / N
    z_n = (mean_n - Y) / math.sqrt(0.37 * 0.63 / N)
    fH = nH / M if M else float("nan")
    z_f = (fH - FH) / math.sqrt(FH * (1 - FH) / M) if M else float("inf")
    z_w = (sw / M) / math.sqrt(1 / 3 / M) if M else float("inf")
    z_w2 = (sw2 / M - 1 / 3) / math.sqrt(4 / 45 / M) if M else float("inf")
    checks = [
        ("① 个数只能是 1 或 2", bad_n == 0, "不合格 %d 次" % bad_n),
        ("② 平均个数 vs %.4f" % Y, abs(z_n) < Z_GATE, "均值 %.5f，z = %+.2f" % (mean_n, z_n)),
        ("③ H 线份额 vs %.5f" % FH, abs(z_f) < Z_GATE and bad_E == 0, "份额 %.5f，z = %+.2f；能量不对 %d 个" % (fH, z_f, bad_E)),
        ("④ ⟨w⟩ vs 0、⟨w²⟩ vs 1/3、单位长度", abs(z_w) < Z_GATE and abs(z_w2) < Z_GATE and bad_dir == 0,
         "z = %+.2f / %+.2f；不归一 %d 个" % (z_w, z_w2, bad_dir)),
        ("⑤ z、wt 原样抄；返回值 = 压栈个数", bad_copy == 0 and bad_ret == 0, "抄错 %d 个；返回值不对 %d 次" % (bad_copy, bad_ret)),
    ]
    for name, good, info in checks:
        ok = ok and good
        msgs.append("  %s %s  %s" % ("✅" if good else "❌", name, info))
    return ok, msgs


# ---------- G2：故意写错的版本（Claude 写，只给判据自检用）----------
def _push(hist, p, E, u, v, w, wt=None):
    hist["pstack"].append({"z": p["z"], "u": u, "v": v, "w": w, "E": E, "wt": p["wt"] if wt is None else wt})


def _bad_factory(kind):
    def emit(p, xs, hist, rng):
        y = (xs["gH_coe"] + xs["gB_coe"]) / xs["Sigma_t"]
        xi = rng.random()
        if kind == "不取整":
            n = int(y)                     # 永远 1 个
        elif kind == "四舍五入":
            n = int(round(y))              # 永远 1 个（1.37 → 1）
        else:
            n = math.floor(y + xi)
        fh = xs["gH_coe"] / (xs["gH_coe"] + xs["gB_coe"])
        for _ in range(n):
            xi2 = rng.random()
            if kind == "选线反了":
                E = E_G_B if xi2 < fh else E_G_H
            else:
                E = E_G_H if xi2 < fh else E_G_B
            u, v, w = isotropic_direction(rng)
            if kind == "方向只往前":
                w = abs(w)
            _push(hist, p, E, u, v, w, wt=1.0 if kind == "权重恒1" else None)
        return n
    return emit


if __name__ == "__main__":
    print("G1：你的 emit_photons（y = %.4f，H 份额 = %.5f，N = %d 次）" % (Y, FH, N))
    ok, msgs = run_g1(emit_photons)
    print("\n".join(msgs))
    print("  → G1 %s\n" % ("全过 ✅" if ok else "没过 ❌"))

    print("G2：自检 —— 5 个故意写错的版本，每个都必须被判错")
    all_caught = True
    for kind in ("不取整", "四舍五入", "选线反了", "方向只往前", "权重恒1"):
        good, _ = run_g1(_bad_factory(kind), seed=7351)
        caught = not good
        all_caught = all_caught and caught
        print("  %s 「%s」%s" % ("✅" if caught else "❌", kind, "被抓到" if caught else "没抓到（判据有漏洞，告诉 Claude）"))
    print("  → G2 %s\n" % ("全响 ✅" if all_caught else "有没响的 ❌"))

    print("G3：请另外 Ctrl+F5 跑 check_split.py（gamma 关着，应逐位相同）\n")

    print("G4（只看）：水模体 + KUR 源，开 gamma，2000 条")
    from lesson07_freegas import run_slab
    import lesson07_freegas as L7
    from lesson11_water import kur_table, KEEP_FACTOR, build_water
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    counts = {"n": 0, "H": 0}
    _orig = L7.emit_photons

    def _counting(p, xs, hist, rng):
        before = len(hist["pstack"])
        n = _orig(p, xs, hist, rng)
        counts["n"] += n
        counts["H"] += sum(1 for g in hist["pstack"][before:] if g["E"] == E_G_H)
        return n
    L7.emit_photons = _counting
    try:
        run_slab(2000, 9137, 15.0, nuclides, None, order, z_src=0.0, implicit=True, NL=30,
                 src=(edges, shares), gamma=True)
    finally:
        L7.emit_photons = _orig
    if counts["n"]:
        print("  平均每条历史生 %.3f 个光子；其中 2.2233 MeV 占 %.4f" % (counts["n"] / 2000, counts["H"] / counts["n"]))
    else:
        print("  一个光子都没生 —— gamma 开关没传到 collide？")
