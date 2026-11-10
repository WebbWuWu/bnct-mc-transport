# -*- coding: utf-8 -*-
# ============================================================
#  check_emit.py  --  M4 第 ② 段判据：记产生点（record_site）+ 生光子压栈（emit_photons）
#  判据归 Claude。2026-11-08 写；🆕 2026-11-10 改成新接口（M4 框架第二轮第 3 条选乙：
#  中子碰撞时只记产生点、一个 ξ 都不抽；中子历史结束后再从产生点生光子），加 G6。
#  运行：Windows Ctrl+F5。G1–G2 几秒；G4、G6 要读真库，合起来约 2–3 分钟。
#  被测接口（你写）：photon.record_site(p, xs, hist)          —— 不收 rng
#                    photon.emit_photons(site, hist, rng)     —— 返回压栈个数
#                    hist 里有 "gsites"（产生点列表）和 "pstack"（光子栈）
#
#  G1  假 xs：gH = 0.0411、gB = 0.0093、Σt = 0.036788 → y = 1.3700、H 份额 = 0.81548（都不是特殊值）
#      调 N = 100000 次（每次一个新的 hist）：先 record_site，再把中子 p 的 z、wt 改掉（模拟中子接着飞），
#      再从刚记的产生点 emit_photons
#      ① 每次个数 ∈ {1, 2}（y = 1.37 → ⌊y + ξ⌋ 只能是 1 或 2）
#      ② 平均个数 vs 1.37：z = (均值 − 1.37) / √(0.37·0.63/N)，|z| < 4
#      ③ H 线份额 vs 0.81548：二项 SE，|z| < 4；能量只能是 E_G_H 或 E_G_B
#      ④ ⟨w⟩ vs 0：SE = √(1/3 / M)，|z| < 4 · ⟨w²⟩ vs 1/3：SE = √(4/45 / M)，|z| < 4 · 每个方向 u² + v² + w² = 1（1e-12）
#      ⑤ 产生点：每次只记一条、记的时候一个随机数都没抽、y / share / z / wt 四个值对；
#         光子的 z、wt = 记的时候中子的值（不是中子后来飞到的地方）；返回值 = 压栈个数
#         抓：存了 p 本身（中子接着飞，位置就跟着变了）、record_site 里偷抽了数、share 记反
#  G2  自检：7 个故意写错的版本，每个都必须被 G1 判错
#  G3  不在这里：关着 gamma 跑 check_split.py（框架判据 P5，逐位）
#  G4  只看：水模体 + KUR 源，开 gamma 跑 2000 条，打印平均每条历史生几个光子、2.2 MeV 占比
#      （11/9 量过 0.9 个 / 历史）。🔴 要等 M4 ③ 第 2 步（run_photons 调 emit_photons）写完才有数
#  G6  🆕 开光子以后，中子那边逐位不变：同一个 seed，gamma=True 和 gamma=False 的 run_slab 返回字典，
#      gamma=False 那份里的每一个键都逐位相同（隐式、模拟各跑一次）
#      抓：光子的随机数在中子历史中途被抽了（11/10 第 3 条 / C8 那个「随机数流错位」）
#      自检：包一层 collide，开 gamma 时多抽一个 ξ（= 改之前那种写法），必须判出不同
#      🔴 只做完第 1 步时，这条是「白过」的（还一个光子都没生）；第 2 步写完再跑才算数
# ============================================================
# ---- 路径（2026-11-10 整理文件夹，Claude 加）：本文件在 my_code 的子文件夹里。
#      下面几行把 my_code 根目录和各子文件夹加进搜索路径，后面的 import 才找得到别的文件。----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import math

from rng import lcg, isotropic_direction
from photon import record_site, emit_photons, E_G_H, E_G_B

GH, GB, ST = 0.0411, 0.0093, 0.036788
Y = (GH + GB) / ST
FH = GH / (GH + GB)
N = 100000
Z_GATE = 4.0
Z0, WT0 = 3.37, 0.6173          # 记产生点时中子的位置、权重
Z1, WT1 = 8.81, 0.2297          # 之后中子「接着飞」改成的值（光子不该拿到这两个）


def _close(a, b):
    return abs(a - b) <= 1e-12 * abs(b)


def run_g1(record, emit, seed=48173):
    """对一对 (record, emit) 跑 G1，返回 (是否全过, 说明列表)。"""
    rng = lcg(seed)
    xs = {"gH_coe": GH, "gB_coe": GB, "Sigma_t": ST}
    msgs, ok = [], True
    cnt_sum, nH, M, sw, sw2 = 0, 0, 0, 0.0, 0.0
    bad_n = bad_dir = bad_copy = bad_E = bad_ret = bad_site = 0
    for i in range(N):
        p = {"z": Z0, "u": 0.28, "v": -0.41, "w": 0.8677, "E": 0.0253, "wt": WT0, "t": 0}
        hist = {"pstack": [], "gsites": []}
        st0 = rng.state
        record(p, xs, hist)
        if len(hist["gsites"]) != 1 or rng.state != st0:
            bad_site += 1
            continue
        site = hist["gsites"][-1]
        try:
            if not (_close(site["y"], Y) and _close(site["share"], FH)
                    and site["z"] == Z0 and site["wt"] == WT0):
                bad_site += 1
        except (KeyError, TypeError):
            bad_site += 1
        p["z"], p["wt"] = Z1, WT1                    # 中子接着飞
        n = emit(site, hist, rng)
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
            if g["z"] != Z0 or g["wt"] != WT0:
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
        ("⑤ 产生点记对、没抽数；光子 z、wt 来自产生点；返回值 = 压栈个数",
         bad_site == 0 and bad_copy == 0 and bad_ret == 0,
         "产生点不对 %d 次；光子 z / wt 不对 %d 个；返回值不对 %d 次" % (bad_site, bad_copy, bad_ret)),
    ]
    for name, good, info in checks:
        ok = ok and good
        msgs.append("  %s %s  %s" % ("✅" if good else "❌", name, info))
    return ok, msgs


# ---------- G2：故意写错的版本（Claude 写，只给判据自检用）----------
def _good_record(p, xs, hist):
    cap = xs["gH_coe"] + xs["gB_coe"]
    hist["gsites"].append({"z": p["z"], "wt": p["wt"], "y": cap / xs["Sigma_t"], "share": xs["gH_coe"] / cap})


def _bad_record(kind):
    def record(p, xs, hist):
        cap = xs["gH_coe"] + xs["gB_coe"]
        if kind == "存了 p 本身":
            p["y"], p["share"] = cap / xs["Sigma_t"], xs["gH_coe"] / cap
            hist["gsites"].append(p)                 # 存的是 p 这个字典本身，中子一飞就跟着变
        else:                                        # "share 记反"
            hist["gsites"].append({"z": p["z"], "wt": p["wt"], "y": cap / xs["Sigma_t"], "share": xs["gB_coe"] / cap})
    return record


def _push(hist, site, E, u, v, w, wt=None):
    hist["pstack"].append({"z": site["z"], "u": u, "v": v, "w": w, "E": E, "wt": site["wt"] if wt is None else wt})


def _bad_emit(kind):
    def emit(site, hist, rng):
        y = site["y"]
        xi = rng.random()
        if kind == "不取整":
            n = int(y)                     # 永远 1 个
        elif kind == "四舍五入":
            n = int(round(y))              # 永远 1 个（1.37 → 1）
        else:
            n = math.floor(y + xi)
        for _ in range(n):
            xi2 = rng.random()
            if kind == "选线反了":
                E = E_G_B if xi2 < site["share"] else E_G_H
            else:
                E = E_G_H if xi2 < site["share"] else E_G_B
            u, v, w = isotropic_direction(rng)
            if kind == "方向只往前":
                w = abs(w)
            _push(hist, site, E, u, v, w, wt=1.0 if kind == "权重恒1" else None)
        return n
    return emit


# ---------- G6：比两份 run_slab 返回字典（逐位）----------
def _same(a, b):
    """逐位相同：数用 ==（两边都是 NaN 也算同）；列表 / 数组逐个元素；字典逐个键。"""
    if isinstance(a, dict):
        return isinstance(b, dict) and set(a) == set(b) and all(_same(a[k], b[k]) for k in a)
    if hasattr(a, "tolist"):
        a = a.tolist()
    if hasattr(b, "tolist"):
        b = b.tolist()
    if isinstance(a, (list, tuple)):
        return isinstance(b, (list, tuple)) and len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


def _diff_keys(off, on):
    """gamma=False 那份里的每个键，在 gamma=True 那份里不在或不同的，列出来。"""
    return [k for k in off if k not in on or not _same(off[k], on[k])]


if __name__ == "__main__":
    print("G1：你的 record_site + emit_photons（y = %.4f，H 份额 = %.5f，N = %d 次）" % (Y, FH, N))
    ok, msgs = run_g1(record_site, emit_photons)
    print("\n".join(msgs))
    print("  → G1 %s\n" % ("全过 ✅" if ok else "没过 ❌"))

    print("G2：自检 —— 7 个故意写错的版本，每个都必须被判错")
    all_caught = True
    bads = [(k, _good_record, _bad_emit(k)) for k in ("不取整", "四舍五入", "选线反了", "方向只往前", "权重恒1")]
    bads += [(k, _bad_record(k), _bad_emit("对")) for k in ("存了 p 本身", "share 记反")]
    for kind, rec, emi in bads:
        good, _ = run_g1(rec, emi, seed=7351)
        caught = not good
        all_caught = all_caught and caught
        print("  %s 「%s」%s" % ("✅" if caught else "❌", kind, "被抓到" if caught else "没抓到（判据有漏洞，告诉 Claude）"))
    good, _ = run_g1(_good_record, _bad_emit("对"), seed=7351)
    print("  %s 对照：Claude 写的正确版本 %s" % ("✅" if good else "❌", "过了" if good else "没过（判据本身错了，告诉 Claude）"))
    all_caught = all_caught and good
    print("  → G2 %s\n" % ("全响 ✅" if all_caught else "有问题 ❌"))

    print("G3：请另外 Ctrl+F5 跑 check_split.py（gamma 关着，应逐位相同）\n")

    import transport as T
    import photon as PH
    from transport import run_slab
    from water import kur_table, KEEP_FACTOR, build_water
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])

    print("G4（只看）：水模体 + KUR 源，开 gamma，2000 条")
    counts = {"n": 0, "H": 0}
    _orig = PH.emit_photons

    def _counting(site, hist, rng):
        before = len(hist["pstack"])
        n = _orig(site, hist, rng)
        counts["n"] += n
        counts["H"] += sum(1 for g in hist["pstack"][before:] if g["E"] == E_G_H)
        return n
    PH.emit_photons = _counting
    try:
        run_slab(2000, 9137, 15.0, nuclides, None, order, z_src=0.0, implicit=True, NL=30,
                 src=(edges, shares), gamma=True)
    finally:
        PH.emit_photons = _orig
    if counts["n"]:
        print("  平均每条历史生 %.3f 个光子；其中 2.2233 MeV 占 %.4f\n" % (counts["n"] / 2000, counts["H"] / counts["n"]))
    else:
        print("  一个光子都没生 —— 第 2 步（run_photons 调 emit_photons）还没写？或者 gamma 没传下去\n")

    print("G6：开光子以后中子那边逐位不变（同一个 seed，gamma=True vs gamma=False，300 条）")
    if not counts["n"]:
        print("  ⚠️ 上面 G4 一个光子都没生 → 本条下面就算过了也是「白过」，第 2 步写完再跑一次")
    g6_ok = True
    for implicit in (True, False):
        kw = dict(z_src=0.0, implicit=implicit, NL=30, src=(edges, shares))
        off = run_slab(300, 52741, 15.0, nuclides, None, order, gamma=False, **kw)
        on = run_slab(300, 52741, 15.0, nuclides, None, order, gamma=True, **kw)
        bad = _diff_keys(off, on)
        g6_ok = g6_ok and not bad
        print("  %s %s：gamma=False 的 %d 个键%s" % ("✅" if not bad else "❌", "隐式" if implicit else "模拟",
              len(off), "全部逐位相同" if not bad else "里有 %d 个不同：%s" % (len(bad), ", ".join(bad[:8]))))
    # 自检：包一层 collide —— 开 gamma 时在碰撞当下多抽一个 ξ（就是改之前那种写法），必须判出不同
    _orig_collide = T.collide

    def _draw_in_collide(*args, **kwargs):
        gamma = args[11] if len(args) > 11 else kwargs.get("gamma", False)
        if gamma:
            args[6].random()                         # 第 7 个位置参数是 rng
        return _orig_collide(*args, **kwargs)
    T.collide = _draw_in_collide
    try:
        kw = dict(z_src=0.0, implicit=True, NL=30, src=(edges, shares))
        off = run_slab(300, 52741, 15.0, nuclides, None, order, gamma=False, **kw)
        on = run_slab(300, 52741, 15.0, nuclides, None, order, gamma=True, **kw)
    finally:
        T.collide = _orig_collide
    caught = bool(_diff_keys(off, on))
    print("  %s 自检：碰撞当下多抽一个 ξ 的版本 %s" % ("✅" if caught else "❌", "被判出不同" if caught else "没判出来（判据有漏洞，告诉 Claude）"))
    print("  → G6 %s" % ("全过 ✅" if g6_ok and caught else "没过 ❌"))
