# -*- coding: utf-8 -*-
# ============================================================
#  lesson09_check2.py  --  第 9 课判据 2：OpenMC 的源能谱 vs KUR 原文数字
#  判据归 Claude，2026-10-15 写。被测：lesson09_openmc.py 里你写的 kur()
#  运行：🔴 只在 WSL 里跑（要 import openmc）：
#        先输 wsl → conda activate /home/siyiwu/miniconda3/envs/openmc
#        → cd /mnt/e/BNCT_GPU/my_code → python lesson09_check2.py
#        不跑输运、不读截面，几秒。全过最后打印「判据 2 全过」。
# ============================================================
#
#  这一条在比什么
#    OpenMC 那边的源能量分布 = 你的 kur() 返回的 Tabular（10 个大箱 × M 个对数子箱，histogram）。
#    它必须和 KUR 原文（Fujiwara 2013 表 1）+ 10/4、10/10 定的箱边界、「箱内 ln E 均匀」是同一个谱。
#    源对不上，判据 3、4 的一切对比都没意义（讲义 §5 最后一段）。
#
#  记号
#    「大箱」   KUR 表的 10 个箱：[0.01, 0.5, 1, 3, 10, 50, 200, 1000, 3000, 1e4, 2e6] eV
#    「子箱」   你在每个大箱里再切的 M 个对数等宽小箱（kur() 里的 edges1）
#    p         Tabular 的第二个参数：概率密度（/eV），子箱概率 = p × 子箱宽（讲义 §5 的 10/14 更正）
#    「路 A / 路 B」 一条判据两边的两种算法，彼此不共用代码
#    SE        二项分布的份额 q 抽 n 次：SE = sqrt(q(1−q)/n)
#    z         (抽到的份额 − 应为) / SE
#
#  路 A：你的 kur() → OpenMC 的 Tabular → 它自己的 .sample() 抽能量
#  路 B：本文件里直接写的原文数字（三群通量 3.0e7 / 7.3e8 / 4.7e7）+ 大箱边界 + ln E 均匀的闭式。
#        🔴 路 B 故意不 import kur_table / KUR_FLUX / EPI_EDGES —— 那些是路 A 的上游，共用就不独立了。
#  两边比的是同一个量：「一个源中子的能量落在某个区间里的概率」。
#
#  ⚠️ 路 A 的抽样是 OpenMC Python 层的 Tabular.sample()（numpy 随机数）。输运时用的是 C++ 那份，
#     公式相同（histogram：先按 p×宽 抽子箱，子箱内 E 均匀；10/12 你实测 + 我读 distribution.cpp 确认过）。
#     判据 2-1、2-2 只看 x、p 本身，不依赖抽样，所以这两条对 C++ 也成立。
#
#  判据（4 条，每条旁边写抓什么）
#    2-1 结构：是 Tabular、interpolation 是 'histogram'；x 严格递增，两端 0.01 eV / 2e6 eV；
#        11 个大箱边界都在 x 里；p ≥ 0；p 的个数 = x 的个数或少一个
#        抓：忘写 interpolation（Tabular 默认是 'linear-linear'，不报错）、边界拼错、
#            子箱跨过大箱边界（切子箱的下标错位，10/13「只在第一个箱对」那种）
#    2-2 每个大箱的概率（子箱 p × 宽 加起来）= 原文份额，相对差 < 1e-9（确定性，不抽样）
#        抓：把份额当密度传（dens 换成 shares1 → 宽箱被放大）、子箱份额没除 M、
#            shares[int(j/M)] 下标写错（份额挪到隔壁箱）
#    2-3 抽 20 万个能量，7 项份额 vs 路 B，每项 |z| < 4
#        前 5 项 = 第 8 课 5b 那五项（三群份额、超热里 < 10 eV、快群里 < 100 keV）；
#        后 2 项 = 两个宽箱「箱内」形状（热群里 < 0.1 eV、[1000, 3000] 里 < 2000 eV —— 讲义 §5 的算例）
#        抓：M = 1（箱内 E 均匀，不是 ln E 均匀 → 快群 < 100 keV 从 43.5% 掉到 4.5%）、
#            子箱不是对数等分（线性等分 → 箱内形状偏）
#        M = 50 的近似误差（子箱内是 E 均匀）：7 项里最大 2.1e-4（快群那项，Claude 10/15 算的），
#        那项的 SE 是 4.6e-3，偏差只占 0.05 SE，不会误报（换 40 个种子跑过，|z| 最大 2.6）
#    2-4 试响（纪律 4）：本文件自己按路 B 造一个「M = 1」的 Tabular（每个大箱只有一个子箱），
#        用同一套 2-3 去比，必须响
#        抓：2-3 太松（样本太少、容差太大）—— 证明它分得出「箱内 E 均匀 vs ln E 均匀」
#        （注意：M = 1 的版本 2-2 是过的 —— 大箱概率没错，错的是箱内形状。所以 2-3 不能省）
#
#  你还要手动试两次（v4 10/15 第 8 项，试完改回）：
#    ① kur() 里 Tabular 的 dens 换成 shares1 → 应该 2-2 响
#    ② kur() 里 M = 50 改成 M = 1     → 应该 2-3 响（2-1、2-2 照样过）
# ============================================================

import math
import openmc
from lesson09_openmc import kur

# ---- 路 B 的数字（直接从原文 / 10/4、10/10 的决定抄，不从你的代码 import）----
FLUX_TH, FLUX_EPI, FLUX_FAST = 3.0e7, 7.3e8, 4.7e7          # KUR 三群通量 cm⁻²s⁻¹，Fujiwara 2013 表 1
BIG_EDGES = [0.01, 0.5, 1.0, 3.0, 10.0, 50.0, 200.0, 1000.0, 3000.0, 1.0e4, 2.0e6]   # eV，10 个大箱

N_SAMPLE = 200000             # 和第 8 课 5b 一样多
SEED = 20261015               # 不取 1、0 这类特殊值
Z_MAX = 4.0
REL_TOL = 1e-9                # 2-2 的容差：几百个乘加，舍入 ~1e-15，留足余量
EDGE_TOL = 1e-9               # 判断「大箱边界在不在 x 里」的相对容差


def big_bin_shares():
    """路 B：10 个大箱各自的概率。热、快各一箱；超热 8 箱按 ln E 分超热份额。"""
    tot = FLUX_TH + FLUX_EPI + FLUX_FAST
    q = [FLUX_TH / tot]
    lo, hi = BIG_EDGES[1], BIG_EDGES[9]                       # 超热段 [0.5, 1e4]
    for i in range(1, 9):
        q.append(FLUX_EPI / tot * math.log(BIG_EDGES[i + 1] / BIG_EDGES[i]) / math.log(hi / lo))
    q.append(FLUX_FAST / tot)
    return q


def expected_items():
    """路 B：2-3 的 7 项。每项 = (名字, 判断落进分子的条件, 判断落进分母的条件, 应为的份额)。"""
    tot = FLUX_TH + FLUX_EPI + FLUX_FAST
    allE = lambda e: True
    th = lambda e: e < 0.5
    epi = lambda e: 0.5 <= e < 1.0e4
    fast = lambda e: e >= 1.0e4
    return [
        ("热群份额 (<0.5 eV)", th, allE, FLUX_TH / tot),
        ("超热群份额", epi, allE, FLUX_EPI / tot),
        ("快群份额 (>10 keV)", fast, allE, FLUX_FAST / tot),
        ("超热群里 <10 eV", lambda e: 0.5 <= e < 10.0, epi, math.log(10 / 0.5) / math.log(1e4 / 0.5)),
        ("快群里 <100 keV", lambda e: 1.0e4 <= e < 1.0e5, fast, math.log(1e5 / 1e4) / math.log(2e6 / 1e4)),
        ("热群里 <0.1 eV", lambda e: e < 0.1, th, math.log(0.1 / 0.01) / math.log(0.5 / 0.01)),
        ("[1e3,3e3] 里 <2e3 eV", lambda e: 1000.0 <= e < 2000.0, lambda e: 1000.0 <= e < 3000.0,
         math.log(2.0) / math.log(3.0)),
    ]


def find_edge(x, e):
    """x 里和 e 相对差 < EDGE_TOL 的那个下标；没有返回 None。"""
    for j in range(len(x)):
        if abs(x[j] / e - 1.0) < EDGE_TOL:
            return j
    return None


def check_structure(dsp):
    """2-1。返回 11 个大箱边界在 x 里的下标。"""
    # [Claude assert] 抓：kur() 返回的不是 Tabular（比如返回了列表、或者 Discrete）
    assert isinstance(dsp, openmc.stats.Tabular), "kur() 应该返回 openmc.stats.Tabular，实际是 %r" % type(dsp)
    # [Claude assert] 抓：忘写 interpolation —— Tabular 默认 'linear-linear'，不报错，子箱内形状全变
    assert dsp.interpolation == "histogram", \
        "interpolation 是 %r，应该是 'histogram'（Tabular 的默认值是 'linear-linear'）" % dsp.interpolation
    x = list(dsp.x)
    p = list(dsp.p)
    # [Claude assert] 抓：p 的个数和 x 对不上（histogram 下 p 可以和 x 一样长，也可以少一个）
    assert len(p) == len(x) or len(p) == len(x) - 1, \
        "x 有 %d 个、p 有 %d 个；histogram 下 p 应该和 x 一样多或少一个" % (len(x), len(p))
    for j in range(len(x) - 1):
        # [Claude assert] 抓：子箱边界不递增（切子箱的下标错位、r 的指数写错）
        assert x[j] < x[j + 1], "x 第 %d、%d 个不递增：%r ≥ %r" % (j, j + 1, x[j], x[j + 1])
    for j in range(len(x) - 1):
        # [Claude assert] 抓：份额 / 密度出现负数
        assert p[j] >= 0.0, "p 第 %d 个是负数：%r" % (j, p[j])
    idx = []
    for e in BIG_EDGES:
        j = find_edge(x, e)
        # [Claude assert] 抓：两端或某个大箱边界拼错、子箱跨过大箱边界
        assert j is not None, "大箱边界 %g eV 不在 x 里（x 两端是 %g、%g，共 %d 个）" % (e, x[0], x[-1], len(x))
        idx.append(j)
    # [Claude assert] 抓：x 两端外面还多出东西（比如多 append 了一个边界）
    assert idx[0] == 0 and idx[-1] == len(x) - 1, \
        "0.01 eV 在 x 第 %d 个、2e6 eV 在第 %d 个；应该分别是第一个和最后一个（共 %d 个）" % (idx[0], idx[-1], len(x))
    return idx


def check_bin_probs(dsp, idx):
    """2-2：每个大箱的概率（子箱 p × 宽 之和）vs 路 B。返回最大相对差。"""
    x = list(dsp.x)
    p = list(dsp.p)
    q_b = big_bin_shares()
    worst = 0.0
    print("   大箱             子箱数   OpenMC 概率      原文份额        相对差")
    for i in range(10):
        s = 0.0
        for j in range(idx[i], idx[i + 1]):
            s += p[j] * (x[j + 1] - x[j])
        rel = abs(s / q_b[i] - 1.0)
        worst = max(worst, rel)
        print("   [%8.3g, %8.3g]  %4d   %.10f   %.10f   %.1e"
              % (BIG_EDGES[i], BIG_EDGES[i + 1], idx[i + 1] - idx[i], s, q_b[i], rel))
        # [Claude assert] 抓：份额当密度传、子箱份额没除 M、份额挪到隔壁箱 —— 都是倍数级的差
        assert rel < REL_TOL, \
            "大箱 [%g, %g] eV：OpenMC 的概率 %.10f，原文份额 %.10f（相对差 %.2e）—— p 是密度吗（份额 / 子箱宽）？" \
            % (BIG_EDGES[i], BIG_EDGES[i + 1], s, q_b[i], rel)
    return worst


def sample_items(dsp, label):
    """2-3：抽 N_SAMPLE 个，7 项份额 vs 路 B。返回最大 |z|（不在这里 assert，2-4 要用）。"""
    E = list(dsp.sample(N_SAMPLE, seed=SEED))
    print("   %s：抽 %d 个（seed %d）" % (label, N_SAMPLE, SEED))
    print("   %-22s %8s   %8s   %8s    z" % ("项", "抽到", "应为", "SE"))
    worst = 0.0
    rows = []
    for name, num, den, q in expected_items():
        n_den = 0
        n_num = 0
        for e in E:
            if den(e):
                n_den += 1
                if num(e):
                    n_num += 1
        got = n_num / n_den
        se = math.sqrt(q * (1.0 - q) / n_den)
        z = (got - q) / se
        rows.append((name, got, q, z))
        worst = max(worst, abs(z))
        print("   %-22s %.5f   %.5f   %.5f  %+7.2f" % (name, got, q, se, z))
    return worst, rows


def m1_tabular():
    """2-4 用：按路 B 造一个「每个大箱只有一个子箱」的 Tabular（= 你的 kur() 里 M = 1）。"""
    q_b = big_bin_shares()
    dens = []
    for i in range(10):
        dens.append(q_b[i] / (BIG_EDGES[i + 1] - BIG_EDGES[i]))
    dens.append(0.0)
    return openmc.stats.Tabular(BIG_EDGES, dens, interpolation="histogram")


if __name__ == "__main__":
    dsp = kur()

    # ---------- 2-1 结构 ----------
    idx = check_structure(dsp)
    n_sub = len(dsp.x) - 1
    print("判据 2-1 结构 ✓（Tabular / histogram；%d 个子箱 = 10 个大箱 × %s；p %d 个、x %d 个）"
          % (n_sub, ("%d" % (n_sub // 10)) if n_sub % 10 == 0 else "不等分", len(dsp.p), len(dsp.x)))

    # ---------- 2-2 大箱概率 ----------
    print()
    worst = check_bin_probs(dsp, idx)
    print("判据 2-2 每个大箱的概率 = 原文份额 ✓（最大相对差 %.1e < %.0e）" % (worst, REL_TOL))

    # ---------- 2-3 抽样 7 项 ----------
    print()
    worst_z, rows = sample_items(dsp, "你的 kur()")
    bad = [r for r in rows if abs(r[3]) >= Z_MAX]
    # [Claude assert] 抓：箱内形状不对（M = 1、子箱线性等分）—— 快群 < 100 keV 那项会偏几十到上百 SE
    assert not bad, \
        "判据 2-3 有 %d 项 |z| ≥ %.0f：%s（2-2 过了说明大箱概率对，错在箱内形状：M 是多少？子箱是对数等分吗？）" \
        % (len(bad), Z_MAX, "；".join("%s 抽到 %.5f 应为 %.5f（%+.1f SE）" % r for r in bad))
    print("判据 2-3 抽样 7 项 ✓（最大 |z| = %.2f < %.0f）" % (worst_z, Z_MAX))

    # ---------- 2-4 试响：M = 1 ----------
    print()
    wrong = m1_tabular()
    wrong_idx = check_structure(wrong)
    check_bin_probs(wrong, wrong_idx)          # M = 1 的大箱概率是对的 —— 这一步应该过
    worst_wrong, rows_wrong = sample_items(wrong, "试响：路 B 自己造的 M = 1")
    # [Claude assert] 抓：2-3 太松，分不出「箱内 E 均匀」和「箱内 ln E 均匀」
    assert worst_wrong >= Z_MAX, \
        "M = 1 的版本 2-3 最大 |z| 只有 %.2f < %.0f：2-3 太松" % (worst_wrong, Z_MAX)
    print("判据 2-4 试响 ✓（M = 1 时 2-2 照样过、2-3 最大 |z| = %.1f，会响）" % worst_wrong)

    print("\n判据 2 全过")
