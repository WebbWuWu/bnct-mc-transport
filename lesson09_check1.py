# -*- coding: utf-8 -*-
# ============================================================
#  lesson09_check1.py  --  第 9 课判据 1：没碰过的中子，OpenMC 每层通量 vs 闭式
#  判据归 Claude，2026-10-14 写。被测：lesson09_openmc.py（用户写）跑出来的 json
#  运行：Windows 里 Ctrl+F5（要 h5py 读截面，所以不在 WSL 跑）。几秒。
#        全过最后打印「判据 1 全过」。
# ============================================================
#
#  这一条在比什么（框架第二轮第 3 问，10/13 定「乙」）
#    真软组织、单能 E = 37.2 eV、垂直入射，OpenMC 只数「一次都没碰过」的中子
#    （CollisionFilter 只留 0 次那一箱）在每层里走的径迹长度。
#    没碰过的中子沿 z 直线走，走过 z 的概率 = e^(−Σt·z)，所以第 k 层（[z_k, z_(k+1)]）：
#        φ_k = ( e^(−Σt·(z_k − z0)) − e^(−Σt·(z_(k+1) − z0)) ) / Σt        单位 cm / 源中子
#    z0 = 1e-6 cm 是 OpenMC 源的位置（第 1 层从 z0 算起，不是从 0）。
#
#  记号
#    Σt        37.2 eV 的宏观总截面，/cm。= Σ_i N_i · σt_i，N_i 是核素 i 的数密度（/cm³）
#    「路 A / 路 B」 一条判据两边的两种算法，彼此不共用代码
#    SE        均值的标准误。OpenMC 的 std_dev 就是 SE（批均值离散算的），不是 σ₁
#    z         (OpenMC − 闭式) / SE
#    「有数的层」 相对误差 SE/均值 < 10% 的层。更深的层没碰过的中子太少（第 15 层约 1.5 个），不比
#
#  路 A：OpenMC（C++ 内核，自己读截面、自己追粒子）
#  路 B：纸上公式 + 你的 xs_lookup.sigma_at 在同一个库、同一个温度上查 σt（不追粒子）
#  两边比的是同一个量：「每源中子、在这一层里、还没碰过时走的径迹长度」，cm。
#    OpenMC 网格 tally 的 flux 是整格的径迹长度之和、没有除体积 —— 讲义 §6 那条，判据 1 就是在验它。
#
#  ⚠️ 层与层之间相关性很强（同一个中子穿过好几层；相邻层 z 的相关系数约 0.7，
#     Claude 10/14 在容器里模拟 600 次量的），所以这里不用 χ²（χ² 要求各层独立）。
#
#  判据（4 条，每条旁边写抓什么）
#    1  结构：30 层；均值、SE 都 ≥ 0；有数的层 ≥ 8 层，而且一层比一层小
#       抓：键取错 / tally 取错、层的顺序反了、网格切在 x 上（30 格几乎一样大）、粒子数远少于 10 万
#    2  有数的层逐层 |z| < 4
#       抓：材料（'ao'、密度没乘 Σw）、源不在 z0、flux 被除了体积、库或温度不一样
#       灵敏度（Claude 模拟 600 次）：Σt 差 3% → 600 次全响；差 1% → 约 20% 会响；对的时候 600 次一次没响
#    3  第 1 层的 SE vs 闭式算出的 SE，比值在 [0.6, 1.5]
#       抓：SE 存错（存成 σ₁、多乘或多除了 √批数）。这时第 2 条会过得「太容易」或「太难」，第 3 条直接看出来
#       闭式：第 1 层每个中子贡献 X = min(自由程, 层厚 − z0)，
#             E[X] = (1 − e^(−ΣL))/Σ，E[X²] = 2/Σ² · (1 − e^(−ΣL)(1 + ΣL))，L = 层厚 − z0
#             SE = sqrt( (E[X²] − E[X]²) / (批数 × 每批粒子数) )
#       模拟 600 次的比值落在 0.80–1.27，所以 [0.6, 1.5] 不会误报；×√10 或 ÷√10 都会响
#    4  试响（纪律 4）：把路 B 的 Σt 故意乘 1.03，第 2 条必须响
#       抓：第 2 条太松（SE 太大、有数的层太少）—— 证明它真分得出 3% 的差别
# ============================================================

import os
import json
import math
from xs_lookup import load_nuclide, sigma_at
from lesson08_dose import NUCS, MASS_FRAC, RHO, M_N, N_A, T, D, NL, OUT_DIR

# ---- 和 lesson09_openmc.py 对齐的设置（你跑的时候改了哪个，这里跟着改）----
E_SRC = 37.2                  # eV，判据 1 的单能源（10/13 Claude 定：KUR [10, 50] 箱内、截面平、非特殊值）
Z_SRC = 1e-6                  # cm，源的 z 位置（lesson09_openmc.py 里 Point((0, 0, 1e-6))）
#  批数、每批粒子数从 json 的 "batches" / "particles" 读（框架 10/14 第 2、6 问：两个都是参数，存进 json）
BARN = 1.0e-24                # cm² / b

# ---- 你存的 json（10/14 第 6 问：照默认）----
J1_JSON = os.path.join(OUT_DIR, "openmc_j1.json")
KEY_MEAN = "flux0_mean"       # 30 个数：没碰过的中子每层通量（tally 原样），从浅到深
KEY_SE = "flux0_se"           # 30 个数：对应的 std_dev（= SE）

REL_SE_MAX = 0.10             # 相对误差 < 10% 才算「有数的层」
MIN_LAYERS = 8                # 有数的层至少这么多（模拟里是 9–10 层）
Z_MAX = 4.0
SE_RATIO_LO, SE_RATIO_HI = 0.6, 1.5
WRONG_FACTOR = 1.03           # 试响用：故意错 3%


def sigma_t_macro(e):
    """路 B 的 Σt：每个核素 N_i × σt_i(e)，相加。返回 (Σt, 每个核素一行的明细)。"""
    total_macro = 0.0
    rows = []
    for name in NUCS:
        E, mts, XS, awr, Q = load_nuclide(name, T)
        sig, tot_b = sigma_at(E, XS, e)
        dens = MASS_FRAC[name] * RHO * N_A / (awr * M_N)      # 和 build_tissue 同一条式子
        part = dens * tot_b * BARN
        rows.append((name, tot_b, dens, part))
        total_macro += part
    return total_macro, rows


def closed_form(sig_t):
    """每层没碰过的中子的通量（cm / 源中子），30 个。"""
    dz = D / NL
    phi = []
    for k in range(NL):
        a = max(k * dz, Z_SRC)
        b = (k + 1) * dz
        phi.append((math.exp(-sig_t * (a - Z_SRC)) - math.exp(-sig_t * (b - Z_SRC))) / sig_t)
    return phi


def first_layer_se(sig_t, batches, particles):
    """第 1 层通量的 SE 的闭式（见文件头第 3 条）。"""
    L = D / NL - Z_SRC
    ex = (1.0 - math.exp(-sig_t * L)) / sig_t
    ex2 = 2.0 / sig_t ** 2 * (1.0 - math.exp(-sig_t * L) * (1.0 + sig_t * L))
    return math.sqrt((ex2 - ex ** 2) / (batches * particles))


def z_scores(mean, se, phi, use):
    return [(mean[k] - phi[k]) / se[k] for k in use]


if __name__ == "__main__":
    # ---------- 路 B：Σt ----------
    sig_t, rows = sigma_t_macro(E_SRC)
    print("E = %.4g eV 的 Σt（路 B，xs_lookup，%s）" % (E_SRC, T))
    for name, tot_b, dens, part in rows:
        print("  %-4s  σt = %10.5f b   N = %.6e /cm³   N·σt = %.6e /cm   占 %.2f%%"
              % (name, tot_b, dens, part, 100.0 * part / sig_t))
    print("  Σt = %.6f /cm   （10/13 凭记忆粗估 ≈ 1.48）" % sig_t)

    phi = closed_form(sig_t)

    # ---------- 路 A：OpenMC 的 json ----------
    assert os.path.exists(J1_JSON), "找不到 OpenMC 的 json：%s（文件名不一样就告诉 Claude）" % J1_JSON
    with open(J1_JSON, "r", encoding="utf-8") as f:
        r = json.load(f)
    assert KEY_MEAN in r and KEY_SE in r, \
        "json 里没有键 %r / %r；实际的键是 %r" % (KEY_MEAN, KEY_SE, sorted(r.keys()))
    mean = r[KEY_MEAN]
    se = r[KEY_SE]
    assert "batches" in r and "particles" in r, "json 里没有 \"batches\" / \"particles\"；实际的键是 %r" % sorted(r.keys())
    batches = r["batches"]
    particles = r["particles"]
    print("\njson：%d 批 × 每批 %d 个" % (batches, particles))

    # ---------- 判据 1：结构 ----------
    assert len(mean) == NL and len(se) == NL, \
        "应该各 %d 层，实际 mean %d 个、se %d 个" % (NL, len(mean), len(se))
    for k in range(NL):
        assert mean[k] >= 0 and se[k] >= 0, "第 %d 层出现负数：mean=%r se=%r" % (k + 1, mean[k], se[k])
    use = [k for k in range(NL) if mean[k] > 0 and se[k] > 0 and se[k] / mean[k] < REL_SE_MAX]
    assert len(use) >= MIN_LAYERS, \
        "有数的层（相对误差 < %.0f%%）只有 %d 层，至少要 %d 层：粒子数够吗？网格切对方向了吗？SE 存大了吗（σ₁、多乘 √批数）？" \
        % (100 * REL_SE_MAX, len(use), MIN_LAYERS)
    for i in range(len(use) - 1):
        k1, k2 = use[i], use[i + 1]
        assert mean[k2] < mean[k1], \
            "有数的层应该一层比一层小，第 %d 层 %.6e ≥ 第 %d 层 %.6e（层序反了？网格切在 x 上？）" \
            % (k2 + 1, mean[k2], k1 + 1, mean[k1])
    print("\n判据 1-1 结构 ✓（有数的层：第 %d–%d 层，共 %d 层）" % (use[0] + 1, use[-1] + 1, len(use)))

    # ---------- 判据 1-2：逐层 z ----------
    z = z_scores(mean, se, phi, use)
    dz = D / NL
    print("\n  层   z 范围 (cm)     OpenMC               闭式           z")
    for i, k in enumerate(use):
        print("  %2d  %5.2f–%5.2f   %.6e ± %.2e   %.6e   %+.2f"
              % (k + 1, k * dz, (k + 1) * dz, mean[k], se[k], phi[k], z[i]))
    worst = max(range(len(z)), key=lambda i: abs(z[i]))
    assert abs(z[worst]) < Z_MAX, \
        "第 %d 层 |z| = %.2f ≥ %.0f：OpenMC %.6e vs 闭式 %.6e（先查对齐清单第 1、3、7、8 格）" \
        % (use[worst] + 1, abs(z[worst]), Z_MAX, mean[use[worst]], phi[use[worst]])
    n1 = sum(1 for v in z if abs(v) < 1)
    print("判据 1-2 逐层 |z| < 4 ✓（最大 |z| = %.2f，在第 %d 层；1σ 内 %d / %d 层，只看不判）"
          % (abs(z[worst]), use[worst] + 1, n1, len(z)))

    # ---------- 判据 1-3：第 1 层 SE vs 闭式 SE ----------
    se_th = first_layer_se(sig_t, batches, particles)
    ratio = se[0] / se_th
    assert SE_RATIO_LO < ratio < SE_RATIO_HI, \
        "第 1 层 SE 比值 = %.3f，不在 (%.1f, %.1f)：OpenMC %.4e vs 闭式 %.4e（SE 存成 σ₁ 了？多乘/除了 √批数？json 里的批数、粒子数和真跑的对得上吗？）" \
        % (ratio, SE_RATIO_LO, SE_RATIO_HI, se[0], se_th)
    print("判据 1-3 第 1 层 SE ✓（OpenMC %.4e / 闭式 %.4e = %.3f）" % (se[0], se_th, ratio))

    # ---------- 判据 1-4：试响 ----------
    phi_wrong = closed_form(sig_t * WRONG_FACTOR)
    z_wrong = z_scores(mean, se, phi_wrong, use)
    max_wrong = max(abs(v) for v in z_wrong)
    assert max_wrong >= Z_MAX, \
        "Σt 故意乘 %.2f 后最大 |z| 只有 %.2f < %.0f：第 2 条太松，分不出 3%% 的差别" \
        % (WRONG_FACTOR, max_wrong, Z_MAX)
    print("判据 1-4 试响 ✓（Σt × %.2f 时最大 |z| = %.1f，第 2 条会响）" % (WRONG_FACTOR, max_wrong))

    print("\n判据 1 全过")
