# -*- coding: utf-8 -*-
# ============================================================
#  check_p7.py  --  M4 判据 P7：中子—光子耦合对拍。KUR 中子源 → 15 cm 水（S(α,β) + 角分布 + 生光子），
#                   每个源中子在每层的光子沉积 vs OpenMC
#  判据归 Claude。2026-11-25 写。运行：Windows Ctrl+F5（openmc (3.10.20) 环境）。
#  主体只读 json，几秒；自检要现跑 2 个故意写错的小版本，预计约 10 分钟（Claude 没在真库上量过速度）。
#  ⚠️ Claude 写的时候没有 h5py 和真库，只在假数据上跑通。第一次跑的完整输出发给 Claude；
#     它自己报错（不是判据响）多半是 Claude 的锅。
#
#  前提：
#    ① Windows：explore\m4_coupled_seeds.py 跑完 → outputs\coupled_SA_N20000_s4001 … s4020.json
#    ② WSL：openmc_compare/openmc_coupled.py 跑完 → outputs\openmc_coupled_led.json、openmc_coupled_ttb.json
#
#  比什么（两条独立的路；先确认比的是同一个量）：
#    路 A（你的程序）  run_slab(gamma=True)：Dg_t（径迹长度 × μ_tr × E）、Dg_c（碰撞处记交出的能量），每层「每个源中子沉积多少 eV」
#    路 B（OpenMC）   heating tally + ParticleFilter(['photon', 'electron', 'positron'])：三份加起来 = 光子带进来的能量
#                     （中子那份不要），每层「每个源中子 eV」。为什么要三种粒子见 load_openmc 里的注释。
#    两边一样：15 cm、30 层、x / y 无限大（OpenMC 反射面）、z 两头真空；KUR 源、z = 0 起、沿 +z；
#              水 + 15 ppm B10，H 用 c_H_in_H2O；同一份中子库、光子库；光子截断 1 keV；电子 LED（= 你的 #40）。
#    已知的差别（都在「生光子」和「光子怎么散」这两处）：
#      · 近似 #32 自由电子康普顿：P6 实测对总量 < 0.16%，单层最大 +0.64%
#      · 你只生 H(n,γ) 2.2233 MeV 和 B10(n,α)→Li7* 0.478 MeV 两条线；OpenMC 按库里所有反应生。
#        KUR 源最高 2 MeV，O16 非弹的 γ（阈值 > 6 MeV）出不来；剩 O16(n,γ)：
#        截面约 0.19 mb vs H 332 mb，N_O = N_H / 2，放出约 4.1 MeV → 估计对光子能量约 +0.05%（OpenMC 那边多）
#      · P6 没检验、这里才第一次检验的：中子在哪被俘获、俘获多少（M3 已对过 gH、B10 (n,α) 的反应率）
#
#  门槛（2026-11-25 14:40 写进文件。这时你这边的 20 个种子和 OpenMC 耦合都还没跑，没有任何耦合的数。）
#    规则和 P6 一样（理由一样，容差来自 #32，O16(n,γ) 那 0.05% 比它小一个量级）：
#    一个量「过」= 下面两条满足一条：
#      ① |z| ≤ 3，z = (你 − OpenMC) / √(SE_你² + SE_OpenMC²)        —— 统计上分辨不出差别
#      ② |相对差| < 容差：总量 0.5%、逐层 1%                         —— 分辨出了，但在已知近似的范围里（照实打印大小）
#    [P7-1] 整块板总量（Dg_t、Dg_c 各判一次）：必须过。
#           总量是 30 层合成的一个数 → 误差用种子间散布（纪律 10）：K = 20 个种子的总量求均值和 SE
#           （t 分布 19 个自由度，|t| > 3 的误报约 0.7%）
#    [P7-2] 逐层（Dg_t、Dg_c 各判一次）：R = SE / 均值 > 0.2 的层不判（纪律 9）；不过的层 ≤ 1 层算过。
#           为什么允许 1 层：全对时 30 层至少一层 |z| > 3 的概率 1 − 0.997³⁰ ≈ 9%；允许 1 层后误报约 0.4%。
#  只看、不判：你 vs OpenMC TTB；OpenMC TTB vs LED（= 耦合情况下近似 #40 的大小）。
#
#  抓什么：生光子的个数错（y = Σ_γ / Σ_t 写错、⌊y + ξ⌋ 写错）、光子权重没跟中子权重、
#          两条线的能量或比例错、生光子的位置错（逐层形状偏）、隐式俘获下漏记某一类碰撞的生光子
#  抓不到什么（照实写）：
#    · 两边共用的东西一起错：中子库、光子库、S(α,β) 表、LED 近似本身
#    · 比分辨力小的偏差（每次打印 3 × SE ÷ 均值）
#    · 中子那边的错 —— 会一起带偏，但那归 M2 / M3 的判据；这里偏了要先回头看 M3 有没有被改坏
#  自检（纪律 4；每个 K = 4 × 5000 条，只看总量 Dg_t 那一条必须不过）：
#    S1  生出的光子权重一律用 1，不跟中子权重（隐式俘获下中子权重一路往下掉 → 光子能量偏高）
#    S2  B10 的 γ 能量写成 H 的 2.2233 MeV（复制没改全那一类；11/24 的 c_sq 同类）
# ============================================================
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", "explore", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json
import math
import time

D_SLAB, NL = 15.0, 30                  # 和 OpenMC 那边一样（dose.D、dose.NL）
N_MAIN = 20000
SEEDS_MAIN = list(range(4001, 4021))   # 和 explore\m4_coupled_seeds.py 一样
K_SELF, N_SELF, SEED_SELF = 4, 5000, 62983
Z_GATE = 3.0
REL_TOL_TOT = 0.005
REL_TOL_LAYER = 0.01
R_MAX = 0.2
MAX_BAD_LAYERS = 1


def out_dir():
    from dose import OUT_DIR           # 放在函数里：假数据试跑时不用导入 dose（Claude 那边没有 matplotlib）
    return OUT_DIR


# ---------------- 读 OpenMC ----------------
def load_openmc(e, folder=None):
    folder = folder or out_dir()
    path = os.path.join(folder, "openmc_coupled_%s.json" % e)
    assert os.path.exists(path), "找不到 %s —— 先在 WSL 跑 openmc_compare/openmc_coupled.py" % path
    with open(path, "r", encoding="utf-8") as f:
        r = json.load(f)
    # [Claude assert] 抓：读错文件（比如读成 P6 的光子单独那份）、层数对不上
    assert r["electron"] == e, "%s 里记的电子处理是 %r，不是 %r" % (path, r["electron"], e)
    # ParticleFilter(['photon', 'electron', 'positron']) 放在 mesh 后面 → 每层 3 个数，排法「层 0 的 光子 / 电子 / 正电子，层 1 的 …」
    #   （OpenMC 展平时最后一个过滤器变得最快）。三份加起来才是「光子带进来的能量」：
    #   OpenMC 的 heating 按「进来的能量 + Q − 出去的粒子」记在当时那个粒子头上，康普顿 / 光电 / 电子对交给电子的那份
    #   记在电子名下（LED 下电子一出生就就地沉积），所以只过滤 photon 会漏掉绝大部分。
    #   SE：三份来自同一批历史、正相关，√(平方和) 会略低估；你这边 SE 大约是 OpenMC 的 5 倍（40 万 vs 1000 万条），影响很小。
    raw_m = r["coupled photon heating"]["mean"]
    raw_s = r["coupled photon heating"]["se"]
    assert len(raw_m) == 3 * NL and len(raw_s) == 3 * NL, \
        "OpenMC 分层那个 tally 有 %d 个数，应该 30 层 × 3 种粒子 = %d 个（ParticleFilter 是不是只写了 'photon'？）" % (len(raw_m), 3 * NL)
    m = [sum(raw_m[3 * j: 3 * j + 3]) for j in range(NL)]
    s = [math.sqrt(sum(x * x for x in raw_s[3 * j: 3 * j + 3])) for j in range(NL)]
    share = [sum(raw_m[i::3]) / sum(raw_m) for i in range(3)]       # 只打印看：光子 / 电子 / 正电子各占多少
    raw_t = r["coupled photon total heating"]["mean"]
    raw_ts = r["coupled photon total heating"]["se"]
    assert len(raw_t) == 3, "total 那个 tally 有 %d 个数，应该只有 3 个（只加 ParticleFilter、不加 mesh）" % len(raw_t)
    tot_v = sum(raw_t)
    tot_s = math.sqrt(sum(x * x for x in raw_ts))
    # [Claude assert] 抓：total 漏了 ParticleFilter（混进中子自己的 heating）→ total 会比 30 层加起来大很多
    assert abs(tot_v / sum(m) - 1.0) < 1e-6, \
        "OpenMC：total %.6e ≠ 30 层加起来 %.6e（差 %+.3f%%）—— total 那个 tally 是不是没加 ParticleFilter？" % (tot_v, sum(m), 100 * (tot_v / sum(m) - 1))
    return {"m": m, "s": s, "tot": tot_v, "tot_s": tot_s, "n": r["batches"] * r["particles"], "share": share}


# ---------------- 你这边：sum / sq → 每层均值、SE ----------------
def layer_stats(s_list, q_list, N):
    """每层：均值 = sum / N，SE = √((sq/N − 均值²) / (N − 1))（第 4 课公式）。"""
    m = [s / N for s in s_list]
    se = [math.sqrt(max(q / N - x * x, 0.0) / (N - 1)) for q, x in zip(q_list, m)]
    return m, se


def combine(results):
    """results = 每个种子 run_slab 返回的字典（含 "N"）。返回 {"t": {...}, "c": {...}, "n": 总条数}。"""
    K = len(results)
    assert K >= 2, "至少要 2 个种子才能用种子间散布"
    Ns = [r["N"] for r in results]
    assert len(set(Ns)) == 1, "各种子条数不一样（%s）—— 下面「直接平均」的前提不成立" % sorted(set(Ns))
    N = Ns[0]
    out = {"n": K * N}
    for key, name in (("t", "Dg_t"), ("c", "Dg_c")):
        per = [layer_stats(r[name + "_sum"], r[name + "_sq"], N) for r in results]
        for k, (m, _) in enumerate(per):
            assert len(m) == NL, "种子 %d 的 %s 有 %d 层，应该 %d 层" % (k, name, len(m), NL)
        m = [sum(per[k][0][j] for k in range(K)) / K for j in range(NL)]                    # 条数相同 → 直接平均
        s = [math.sqrt(sum(per[k][1][j] ** 2 for k in range(K))) / K for j in range(NL)]   # 独立种子：SE 平方相加
        tots = [sum(per[k][0]) for k in range(K)]                                           # 每个种子一个总量
        tm = sum(tots) / K
        ts = math.sqrt(sum((x - tm) ** 2 for x in tots) / (K - 1) / K)                      # 种子间散布 → 均值的 SE
        out[key] = {"m": m, "s": s, "tot": tm, "tot_s": ts}
    return out


def load_mine(folder=None):
    folder = folder or out_dir()
    results, missing = [], []
    for seed in SEEDS_MAIN:
        path = os.path.join(folder, "coupled_SA_N%d_s%d.json" % (N_MAIN, seed))
        if not os.path.exists(path):
            missing.append(seed)
            continue
        with open(path, "r", encoding="utf-8") as f:
            r = json.load(f)
        # [Claude assert] 抓：读到没开 gamma 的结果（Dg 全是 0）
        assert sum(r["Dg_t_sum"]) > 0, "%s 里 Dg_t 全是 0 —— 跑的时候 gamma=True 开了吗？" % path
        results.append(r)
    assert not missing, "这些种子还没跑：%s —— 先跑完 explore\\m4_coupled_seeds.py" % missing
    return combine(results)


# ---------------- 判（和 P6 同一套） ----------------
def judge(x, sx, y, sy, tol=REL_TOL_TOT):
    d = x - y
    se = math.sqrt(sx * sx + sy * sy)
    z = d / se if se > 0 else float("inf")
    rel = d / y if y != 0 else float("inf")
    if abs(z) <= Z_GATE:
        return True, z, rel, "✅"
    if abs(rel) < tol:
        return True, z, rel, "✅（分辨出了，< %.1f%%，照实记）" % (100 * tol)
    return False, z, rel, "❌"


def compare(mine, ref, label, verbose=True):
    ok_t, z, rel, tag = judge(mine["tot"], mine["tot_s"], ref["tot"], ref["tot_s"])
    if verbose:
        print("  [P7-1] %s 总量：你 %.5e ± %.2e · OpenMC %.5e ± %.2e · 差 %+.3f%% · z = %+.2f  %s"
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
        print("  [P7-2] %s 逐层（层 · OpenMC · 你 · 差 · z）：" % label)
        print("\n".join(lines))
        print("         不过的层：%s（允许 ≤ %d 层）；不判的层：%s  %s"
              % (bad or "无", MAX_BAD_LAYERS, skipped or "无", "✅" if ok_l else "❌"))
    return ok_t, ok_l


def look_only(mine, ref, label):
    rel_tot = mine["tot"] / ref["tot"] - 1.0
    print("  %s：总量差 %+.3f%%；逐层差（%%）：" % (label, 100 * rel_tot))
    print("    " + " ".join("%+.2f" % (100 * (a / b - 1.0)) for a, b in zip(mine["m"], ref["m"])))


# ---------------- 自检：现跑 2 个故意写错的小版本 ----------------
def run_small(seeds, N):
    from transport import run_slab
    from water import build_water
    from dose import kur_table, KEEP_FACTOR
    from thermal import load_thermal, SAB_PATH, T_TAG
    from elastic_angle import load_elastic_angle, LIB, NUCS
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    thermal = load_thermal(SAB_PATH, T_TAG)
    angle = load_elastic_angle(LIB, NUCS)
    results = []
    for seed in seeds:
        r = run_slab(N, seed, D_SLAB, nuclides, None, order, z_src=0.0, implicit=True, NL=NL,
                     src=(edges, shares), thermal=thermal, angle=angle, gamma=True)
        r["N"] = N
        results.append(r)
    return combine(results)


def selftest(led):
    import photon as PH
    orig_emit, orig_EB = PH.emit_photons, PH.E_G_B

    def emit_w1(site, hist, rng):                 # S1：光子权重一律 1（随机数调用一个不变，只事后改权重）
        n0 = len(hist["pstack"])
        n = orig_emit(site, hist, rng)
        for g in hist["pstack"][n0:]:
            g["wt"] = 1.0
        return n

    seeds = [SEED_SELF + 7919 * k for k in range(K_SELF)]
    all_ok = True
    for name, patch in (("S1 光子权重一律 1", lambda: setattr(PH, "emit_photons", emit_w1)),
                        ("S2 B10 的 γ 用 2.2233 MeV", lambda: setattr(PH, "E_G_B", PH.E_G_H))):
        t0 = time.time()
        patch()
        try:
            # run_photons 在 photon 模块里按名字找 emit_photons、emit_photons 按名字找 E_G_B → 改 PH 上的就生效
            mine = run_small(seeds, N_SELF)
        finally:
            PH.emit_photons, PH.E_G_B = orig_emit, orig_EB
        ok_t, z, rel, _ = judge(mine["t"]["tot"], mine["t"]["tot_s"], led["tot"], led["tot_s"])
        fired = not ok_t
        all_ok = all_ok and fired
        print("  %-24s Dg_t 总量差 %+.2f%%，z = %+.1f  %s（%.0f 秒）"
              % (name, 100 * rel, z, "✅ 响了" if fired else "❌ 没响（判据不够灵，告诉 Claude）", time.time() - t0), flush=True)
    return all_ok


if __name__ == "__main__":
    t_start = time.time()
    led = load_openmc("led")
    ttb = load_openmc("ttb")
    mine = load_mine()

    print("P7：KUR 中子 → 15 cm 水，每个源中子在每层的光子沉积。OpenMC LED %d 个中子；你这边 %d 个种子 × %d 条"
          % (led["n"], len(SEEDS_MAIN), N_MAIN))
    print("  OpenMC LED 的沉积记在谁名下：光子 %.2f%% · 电子 %.2f%% · 正电子 %.2f%%（只看）" % tuple(100 * x for x in led["share"]))
    print()
    ok = True
    for key, label in (("t", "Dg_t（径迹）"), ("c", "Dg_c（碰撞）")):
        ok_t, ok_l = compare(mine[key], led, label)
        ok = ok and ok_t and ok_l
        print()

    print("只看、不判：")
    look_only(mine["t"], ttb, "你的 Dg_t vs OpenMC TTB")
    look_only(ttb, led, "OpenMC TTB vs LED（= 耦合时的 #40）")
    print()

    print("自检：2 个故意写错的版本必须响（现跑，约 10 分钟）")
    ok_s = selftest(led)

    ok = ok and ok_s
    print("\n  → %s    用时 %.0f 秒" % ("✅ P7 全过" if ok else "❌ P7 没全过", time.time() - t_start))
    assert ok, "P7 没全过：看上面哪一行是 ❌（总量 ❌ = 整体偏；逐层 ❌ = 形状偏；自检 ❌ = 判据不够灵）"
