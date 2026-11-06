# -*- coding: utf-8 -*-
# ============================================================
#  check_lookup.py  --  M3 判据 T5：lookup_xs 接入 S(α,β) 截面对不对
#  Claude 2026-11-05 写（判据归 Claude，成品直接跑，不用默写）
#
#  用法：在 VS Code 里打开本文件，Ctrl+F5（Windows 的 openmc (3.10.20) 环境）。
#        它自己搭水模体（lesson11_water.build_water）、读 c_H_in_H2O（thermal.load_thermal，约 629 MB，要等一会儿），
#        然后先跑 check_T5（正式判），再跑 check_T5_selftest（4 种故意写错的必须响、响对的那一条，纪律 4）。
#
#  两条路：
#    你的路   = lookup_xs(nuclides, idx, E, thermal) 返回的字典
#    判据的路 = ① lookup_xs(nuclides, idx, E)（不传 thermal = M2 的老路）· ② 判据自己用 h5py 读 c_H_in_H2O 的
#               inelastic/xs 再插值（不经过 load_thermal）· ③ 用你返回的 tot_ls 重算 Σ = N·σ_tot·BARN
#  3 条（按顺序跑）：
#    [T5-1] 10 eV 以上（10.5 eV、1 MeV）：返回的字典和不传 thermal 的逐位相同          抓：< 写成 >、条件少了能量那一半
#    [T5-2] 10 eV 以下（0.0253 eV、3.7 eV）：除了 H1 自己的 MT=2、tot、Σ 和 Sigma_t，     抓：条件挂错核素（改了 O16 / B10）、
#           别的一个不动（含 MT=102 系数、kerma 系数）；调完再不传 thermal 查一次，            改坏了共用的截面库（下一步还在用）
#           和第一次逐位相同
#    [T5-3] 10 eV 以下：H1 的 MT=2 = 判据自己读库插值（0.0253 eV 还要和 10/31 的 51.69 b 对上）；   抓：先用后改（Σ 没跟着变）、
#           tot 变化量 = 新 − 旧；H1 的 Σ = N·tot·BARN；Sigma_t = 各核素 Σ 之和                        只加新不减旧、插值插错数组
#    🆕 [T5-4]（11/6 加）S(α,β) 表的最低能量以下（取表下限的一半）：和不传 thermal 的逐位相同   抓：表外也换（11/6 用户定乙：
#           表下限以下截面、运动学都走自由 H —— 起因：种子 1011 在表下限以下「修正 E′」抽出负能量，sigma_at_ext 的断言响了）
#  门槛在看数据之前定：「逐位相同」用 ==；算出来的量相对差 < 1e-12（只差浮点舍入）；51.69 容差 0.01 b
#  盲区：只查 4 个能量点；10 eV 两边截面接不接得上（物理上差 < 1%）不在这里查 —— 小跑通看通量谱有没有台阶
# ============================================================

import os
import copy
import h5py
import numpy as np

SAB_PATH = os.path.join(r"E:\BNCT_GPU\lib80x_hdf5", "c_H_in_H2O.h5")
T_TAG = "294K"
BARN = 1.0e-24
E_ABOVE = [10.5, 1.0e6]        # eV，T5-1
E_BELOW = [0.0253, 3.7]        # eV，T5-2、T5-3
SIG_0253_REF = 51.69           # b，10/31 读出（两位小数）
REL = 1e-12


def t5_read_xs(path=SAB_PATH, T=T_TAG):
    """判据自己的路：直接读 inelastic/xs（两行：能量、截面），不经过 load_thermal。"""
    with h5py.File(path, "r") as f:
        root = f["c_H_in_H2O"]
        xs = root[T]["inelastic/xs"][:]
        E_max = float(root.attrs["energy_max"])
    return np.asarray(xs[0], float), np.asarray(xs[1], float), E_max


def _same(a, b):
    """逐位相同：列表 / 数组逐个 ==，字典逐键，数字 ==。"""
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple, np.ndarray)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return a == b


def _rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def _diff_keys(new, old, nuclides):
    """列出 new 和 old 不一样的地方（给 assert 消息用）。"""
    out = []
    for k in old:
        if k in ("sig_ds", "sig_al", "tot_ls"):
            for j in range(len(old[k])):
                if not _same(new[k][j], old[k][j]):
                    out.append("%s[%s]" % (k, nuclides[j]["name"]))
        elif not _same(new[k], old[k]):
            out.append(k)
    return out


def check_T5(lookup_xs, nuclides, idx, thermal, quiet=False):
    xs_E, xs_sig, E_max = t5_read_xs()
    H = idx["ind_H"]
    assert H is not None, "[T5-0] nuclides 里没有 H1 —— 判据要在水模体上跑"
    pl = nuclides[H]["mts"].index(2)
    dens_H = nuclides[H]["dens"]

    for E in E_ABOVE:
        new = lookup_xs(nuclides, idx, E, thermal)
        old = lookup_xs(nuclides, idx, E)
        # [T5-1] 10 eV 以上不该换：和不传 thermal 的逐位相同
        assert _same(new, old), "[T5-1] E = %g eV（> energy_max %g）不该换，可这些变了：%s" % (E, E_max, _diff_keys(new, old, nuclides))
    if not quiet:
        print("T5-1 ✅ 10.5 eV、1 MeV：和不传 thermal 逐位相同（没换）")

    allowed = {"sig_ds[H1]", "sig_al[H1]", "tot_ls[H1]", "Sigma_t"}
    res = {}
    for E in E_BELOW:
        old = lookup_xs(nuclides, idx, E)
        new = lookup_xs(nuclides, idx, E, thermal)
        old2 = lookup_xs(nuclides, idx, E)
        changed = set(_diff_keys(new, old, nuclides))
        # [T5-2] 只许 H1 自己的那几样变；别的核素、MT=102 / kerma 系数一个不动
        assert changed <= allowed, "[T5-2] E = %g eV：不该变的变了：%s（只许变 %s）" % (E, sorted(changed - allowed), sorted(allowed))
        # [T5-2] H1 的 MT=2 以外的道也不许变
        for m in range(len(old["sig_ds"][H])):
            if m != pl:
                assert new["sig_ds"][H][m] == old["sig_ds"][H][m], \
                    "[T5-2] E = %g eV：H1 的第 %d 个道（MT=%d）变了，只许改 MT=2" % (E, m, nuclides[H]["mts"][m])
        # [T5-2] 调完再走老路，结果和第一次一样 → 没有改坏共用的截面库
        assert _same(old2, old), "[T5-2] E = %g eV：传过 thermal 以后再不传，结果变了：%s —— 改到了共用的列表（截面库被改坏）" \
            % (E, _diff_keys(old2, old, nuclides))
        res[E] = (new, old)
    if not quiet:
        print("T5-2 ✅ 0.0253 eV、3.7 eV：只有 H1 的 MT=2、tot、Σ、Sigma_t 变了；再走老路结果不变（截面库没被改坏）")

    for E in E_BELOW:
        new, old = res[E]
        want = float(np.interp(E, xs_E, xs_sig))
        got = new["sig_ds"][H][pl]
        # [T5-3] H1 的 MT=2 = 判据自己读库插值
        assert _rel(got, want) < REL, "[T5-3] E = %g eV：H1 的 MT=2 = %.6g b，判据读库插值 %.6g b" % (E, got, want)
        if E == 0.0253:
            # [T5-3] 物理锚：0.0253 eV 51.69 b（10/31）
            assert abs(got - SIG_0253_REF) < 0.01, "[T5-3] 0.0253 eV：H1 的 MT=2 = %.4f b，10/31 读出来是 %.2f b" % (got, SIG_0253_REF)
        d_tot = new["tot_ls"][H] - old["tot_ls"][H]
        d_sig = got - old["sig_ds"][H][pl]
        # [T5-3] tot 减旧加新：变化量 = 新 MT=2 − 旧 MT=2
        assert abs(d_tot - d_sig) <= REL * max(abs(new["tot_ls"][H]), 1.0), \
            "[T5-3] E = %g eV：tot_ls[H1] 变了 %+.6g b，MT=2 变了 %+.6g b —— 两个应一样（是不是只加新、没减旧？）" % (E, d_tot, d_sig)
        sig_al_H = dens_H * new["tot_ls"][H] * BARN
        # [T5-3] 先改后用：H1 的 Σ 要用新的 tot 算
        assert _rel(new["sig_al"][H], sig_al_H) < REL, \
            "[T5-3] E = %g eV：sig_al[H1] = %.6g /cm，用新 tot 算应是 %.6g /cm —— 替换是不是放在了 Sigma_j 那一行后面？" \
            % (E, new["sig_al"][H], sig_al_H)
        # [T5-3] Sigma_t = 各核素 Σ 之和（用新的）
        assert _rel(new["Sigma_t"], sum(new["sig_al"])) < REL, \
            "[T5-3] E = %g eV：Sigma_t = %.6g /cm，各核素 Σ 加起来 %.6g /cm" % (E, new["Sigma_t"], sum(new["sig_al"]))
        if not quiet:
            print("T5-3 ✅ E = %-6g eV · H1 MT=2 %.4f b（自由 H %.4f b）· tot_ls[H1] %+.4f b · Sigma_t %.5f → %.5f /cm（× %.3f）"
                  % (E, got, old["sig_ds"][H][pl], d_tot, old["Sigma_t"], new["Sigma_t"], new["Sigma_t"] / old["Sigma_t"]))

    # [T5-4] 11/6：表的最低能量以下也不该换（乙方案）—— 用判据自己读的能量网格（xs_E），不用 thermal 里的
    E_low = 0.5 * float(xs_E[0])
    new = lookup_xs(nuclides, idx, E_low, thermal)
    old = lookup_xs(nuclides, idx, E_low)
    assert _same(new, old), "[T5-4] E = %g eV（表下限 %g eV 的一半）不该换，可这些变了：%s" % (E_low, xs_E[0], _diff_keys(new, old, nuclides))
    if not quiet:
        print("T5-4 ✅ E = %g eV（表下限 %g eV 的一半）：和不传 thermal 逐位相同（表外不换）" % (E_low, xs_E[0]))
        print("判据 T5 过")


def _apply(xs, nuclides, j, sig_new, minus_old=True):
    """判据自己做一次「换 MT=2」（只给自检造错用）。"""
    xs = copy.deepcopy(xs)
    pl = nuclides[j]["mts"].index(2)
    old = xs["sig_ds"][j][pl]
    xs["sig_ds"][j][pl] = sig_new
    xs["tot_ls"][j] = xs["tot_ls"][j] + sig_new - (old if minus_old else 0.0)
    xs["sig_al"][j] = nuclides[j]["dens"] * xs["tot_ls"][j] * BARN
    xs["Sigma_t"] = sum(xs["sig_al"])
    return xs


def check_T5_selftest(lookup_xs, nuclides, idx, thermal):
    """纪律 4：5 种故意写错的必须响，而且响对的那一条。都是包在你的 lookup_xs 外面造错，不改你的文件。"""
    xs_E, xs_sig, E_max = t5_read_xs()
    H = idx["ind_H"]
    O = [n["name"] for n in nuclides].index("O16")
    sig = lambda E: float(np.interp(E, xs_E, xs_sig))

    def gt(nuc, ix, E, th=None):                 # < 写成 >：10 eV 以上换、以下不换
        old = lookup_xs(nuc, ix, E)
        return _apply(old, nuc, H, sig(E)) if (th is not None and E > E_max) else old

    def on_O(nuc, ix, E, th=None):               # 条件挂错核素：换到了 O16 头上
        old = lookup_xs(nuc, ix, E)
        return _apply(old, nuc, O, sig(E)) if (th is not None and E < E_max) else old

    def late(nuc, ix, E, th=None):               # 先用后改：MT=2、tot 换了，Σ、Sigma_t 还是旧的
        new = lookup_xs(nuc, ix, E, th)
        if th is None:
            return new
        old = lookup_xs(nuc, ix, E)
        new = copy.deepcopy(new)
        new["sig_al"] = list(old["sig_al"])
        new["Sigma_t"] = old["Sigma_t"]
        return new

    def below(nuc, ix, E, th=None):              # 11/6 前的写法：表下限以下也换（只看 E < E_max）
        old = lookup_xs(nuc, ix, E)
        return _apply(old, nuc, H, sig(E)) if (th is not None and E < E_max) else old

    def no_minus(nuc, ix, E, th=None):           # 只加新不减旧
        old = lookup_xs(nuc, ix, E)
        return _apply(old, nuc, H, sig(E), minus_old=False) if (th is not None and E < E_max) else old

    cases = [("< 写成 >", gt, "[T5-1]"),
             ("挂错核素（O16）", on_O, "[T5-2]"),
             ("先用后改", late, "[T5-3]"),
             ("只加新不减旧", no_minus, "[T5-3]"),
             ("表下限以下也换", below, "[T5-4]")]
    for what, bad, tag_want in cases:
        try:
            check_T5(bad, nuclides, idx, thermal, quiet=True)
        except AssertionError as e:
            tag = str(e)[:6]
            assert tag == tag_want, "自检：「%s」响的是 %s，应是 %s" % (what, tag, tag_want)
            print("自检 ✅ %-14s → %s 响了" % (what, tag))
            continue
        raise AssertionError("自检：「%s」没有任何一条响 —— 判据有漏洞，截图给 Claude" % what)
    print("判据 T5 自检全过：5 种错都响了，而且响的是对的那一条")


if __name__ == "__main__":
    import time
    from lesson07_freegas import lookup_xs, prepare_nuclides
    from lesson11_water import build_water
    from lesson08_dose import kur_table, KEEP_FACTOR
    from thermal import load_thermal

    t0 = time.time()
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    idx = prepare_nuclides(nuclides)
    thermal = load_thermal(SAB_PATH, T_TAG)
    print("水模体 + c_H_in_H2O 读好（%.0f s）" % (time.time() - t0))
    check_T5(lookup_xs, nuclides, idx, thermal)
    check_T5_selftest(lookup_xs, nuclides, idx, thermal)
