# -*- coding: utf-8 -*-
# ============================================================
#  check_split.py  --  run_slab 拆分的逐位回归判据（M2 前置，10/23）
#  判据归 Claude。2026-10-23 写。被测：拆完以后的 run_slab（你写）。
#
#  运行：Windows，Ctrl+F5（要真库的两个情况才会跑；没有 h5py / 真库时只跑假库三个情况）。
#        自动分两种模式，看 outputs\split_ref\ 里有没有基准：
#          没有 → 「存基准」模式：用【现在、拆之前】的 run_slab 跑 5 个情况，存成 json
#          有   → 「检查」模式：用【拆完以后】的 run_slab 再跑一遍，和基准逐位比
#        🔴 所以第一次跑必须在你动 run_slab 之前。之后每拆完一块就跑一次。
#        用时：假库三个约 10 秒；真库两个约 15–30 秒（含读截面）。
# ============================================================
#
#  这一条在比什么
#    拆分只是「把代码搬到几个小函数里」，物理一点不变 —— 所以同一个种子、同一组输入，
#    拆前和拆后的 run_slab 返回的字典必须【逐位相同】（== ，不是「差不多」）。
#    任何一位不同，说明拆的时候改了算法：最常见的是随机数的调用顺序变了
#    （比如把某个 rng.random() 挪到了另一个 rng.random() 前面 —— 分布不变，但逐位全变）。
#
#  两条路
#    路 A：拆之前的 run_slab（存进 json 的基准）
#    路 B：拆之后的 run_slab（现在 import 进来的）
#    比的是同一个量：同一个函数签名、同样 5 组输入、返回字典的每一个键。
#
#  5 个情况（为什么是这 5 个：每个分支至少有一个情况走到）
#    | 情况 | 库   | 模式       | 源               | 特别的           | 走到的分支 |
#    | F1   | 假库 | 模拟       | 单能 1.37 MeV     | z_src = 中间     | 抽核素、抽道、吸收、awr=None 各向同性 |
#    | F2   | 假库 | 隐式       | 三箱表格          | z_src = 0        | 权重、轮盘赌、按散射份额抽核素 |
#    | F3   | 假库 | 模拟       | 单能 0.73 MeV     | t_cut = 0.29     | 时间截断、E_cut |
#    | R1   | 真库 | 隐式       | KUR 表           | 软组织，30 层    | 和第 8 课交付同一套输入（N 小） |
#    | R2   | 真库 | 模拟       | KUR 表           | 软组织           | 真库 + 抽道（order） |
#    假库四个核素：H1（自由气体）、O16（awr = None → 实验室系各向同性）、B10（800 / 801）、N14（600）
#    —— 这样 run_slab 开头找 B、N、102 下标的那段也被走到。
#
#  验证（4 条；每条一行「抓什么」）
#    S1 存基准时，lesson07_freegas.py 必须还是拆之前那一版（md5 和 10/22 收工时一样）
#       抓：先动了 run_slab 才存基准 —— 那就是拿新版和新版比，永远过，判据零分辨率
#    S2 存基准时，同一个种子跑两遍逐位相同（假库三个）
#       抓：程序里有不受种子控制的东西（那样「逐位回归」本身没意义）
#    S3 试响（纪律 4）：把模块常数 Q_val 改大 3e-12（相对），再跑 F2，
#       必须【恰好】Bc_sum、Bc_sq 两个键报差、其余全同
#       抓：比较函数不够严（用了 isclose、漏比某些键、只比了长度）
#    C1 检查时，每个情况：基准里的每个键都在、都逐位相同（列表逐个元素，字典逐个键）
#       10/23 晚改：拆分做完以后，返回字典里【新加的键】只列出来、不判（M2 起要加 kerma 等计数）；少了老键照样报错
#       抓：拆的时候随机数调用顺序变了、某个计数漏加 / 加了两次、层号错、返回值少一项 / 多一项
#
#  Claude 10/23 在故意改坏的版本上试过（假库三个情况）：
#    | 改坏                                         | 结果 |
#    | 抽核素、抽道两个 rng.random() 顺序对调          | ✗ F1、F3 报差（F2 隐式不走这里，对）|
#    | 拆层循环里漏加 a_g（(n,γ) 计数）                | ✗ 三个情况都报 ag_sum、ag_sq |
#    | 轮盘赌存活后 wt_hi 乘 1.0000001                 | ✗ F2 报 16 个键 |
#    | 返回字典多一项                                  | ✗ 三个情况都报「新版多了这个键」|
#    | 🔴 删掉 w == 0 分支里的 N_ar 那一行              | ✓ 没响 —— 盲区，见下 |
#  盲区：w == 0.0（方向正好和板面平行）那个分支，随机情况下几乎走不到（要 rng.random() 正好等于 0.5）。
#        本判据看不见它 —— 拆的时候那一段【照搬，一个字别改】。
#
#  🔴 拆分期间的约定（这条判据成立的前提）
#    · run_slab 的【参数和返回字典】先不变。改签名、改返回值是拆完、回归过了以后的下一步。
#    · 新的 run_slab 放在哪个文件你定：改下面 NEW_MODULE 那一行就行（默认还在 lesson07_freegas）。
# ============================================================

import os
import sys
import json
import math
import hashlib
import importlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

NEW_MODULE = "lesson07_freegas"              # 拆完的 run_slab 在哪个模块（你定，改这一行）
REF_DIR = os.path.join(HERE, "outputs", "split_ref")
BASE_FILE = os.path.join(HERE, "lesson07_freegas.py")
BASE_MD5 = "8acbeec8c42250f3b720346164929e2e"   # 10/22 收工时 lesson07_freegas.py 的 md5（Claude 10/23 15:10 读的）


# ------------------------------------------------------------
#  假库：4 个点的能量网格，截面按「1/v + 平台」随手给的非特殊值
# ------------------------------------------------------------
E_FAKE = np.array([1.0e-5, 0.0253, 1.0e2, 2.0e7])


def fake_nuclides():
    def nuc(name, mts, xs_rows, awr, Q, dens, kT, heat):
        return {"E": E_FAKE, "mts": mts, "XS": [np.array(r) for r in xs_rows],
                "awr": awr, "Q_value": Q, "name": name, "dens": dens, "kT": kT,
                "heat": np.array(heat)}     # 10/23 晚加：M2 起 run_slab 要读 "heat"（MT=301，eV·b）；假库的数按「俘获反冲 + 弹性反冲」量级随手给
    nuclides = [
        nuc("H1", [2, 102],
            [[30.1, 29.7, 20.4, 0.49], [16.75, 0.3326, 5.29e-3, 1.18e-5]],
            0.999167, [0.0, 2224566.0], 6.672e22, 0.025301, [22013.0, 438.03, 1021.0, 4.9e6]),
        nuc("O16", [2, 102],
            [[3.97, 3.91, 3.85, 1.47], [3.0e-3, 1.9e-4, 3.0e-6, 5.0e-5]],
            None, [0.0, 4143080.0], 3.336e22, 0.025301, [0.27, 0.0508, 42.6, 1.6e6]),
        nuc("B10", [2, 102, 800, 801],
            [[2.41, 2.31, 2.13, 0.93], [25.1, 0.50, 7.9e-3, 1.8e-5],
             [12170.0, 242.1, 3.85, 8.6e-3], [181100.0, 3602.3, 57.3, 0.128]],
            9.926921, [0.0, 11454000.0, 2789680.0, 2311890.0], 1.13e18, 0.025301, [4.53e11, 9.003e9, 1.43e8, 1.2e6]),
        nuc("N14", [2, 102, 600],
            [[11.6, 10.05, 10.0, 1.21], [3.77, 0.0750, 1.19e-3, 2.7e-6],
             [92.0, 1.83, 0.029, 0.11]],
            13.882814, [0.0, 10833300.0, 625870.0], 1.47e21, 0.025301, [5.8e7, 1.15e6, 1.9e4, 2.3e6]),
    ]
    order = []
    for n in nuclides:            # 和 build_tissue 一样：按 0.0253 eV（第 1 个格点）截面从大到小排下标
        order.append(sorted(range(len(n["mts"])), key=lambda t: n["XS"][t][1], reverse=True))
    return nuclides, order


FAKE_SRC = ([0.0213, 0.47, 913.0, 1.73e6], [0.07, 0.61, 0.32])

# 情况名 → run_slab 的参数（N, seed, d, E_src 之外的都按名字传）
FAKE_CASES = {
    "F1": dict(N=400, seed=2027, d=7.3, E_src=1.37e6, kw=dict(NL=11)),
    "F2": dict(N=400, seed=331, d=9.1, E_src=None, kw=dict(z_src=0.0, implicit=True, NL=15, src=FAKE_SRC)),
    "F3": dict(N=300, seed=4421, d=7.3, E_src=0.73e6, kw=dict(t_cut=0.29, NL=11)),
}
REAL_CASES = {
    "R1": dict(N=1000, seed=1010, d=15.0, kw=dict(z_src=0.0, implicit=True, NL=30)),
    "R2": dict(N=300, seed=77, d=15.0, kw=dict(z_src=0.0, implicit=False, NL=30)),
}


def real_available():
    try:
        import h5py  # noqa: F401
    except ImportError:
        return False
    from xs_lookup import LIB
    return os.path.exists(os.path.join(LIB, "H1.h5"))


def run_case(mod, name, inputs):
    """跑一个情况，返回 run_slab 的字典。inputs = (nuclides, order, src 或 None)。"""
    nuclides, order, real_src = inputs
    if name in FAKE_CASES:
        c = FAKE_CASES[name]
        return mod.run_slab(c["N"], c["seed"], c["d"], nuclides, c["E_src"], order, **c["kw"])
    c = REAL_CASES[name]
    return mod.run_slab(c["N"], c["seed"], c["d"], nuclides, None, order, src=real_src, **c["kw"])


def prepare_inputs():
    """每个情况的输入只准备一次（真库读截面慢）。"""
    fake = fake_nuclides()
    inputs = {k: (fake[0], fake[1], None) for k in FAKE_CASES}
    if real_available():
        from lesson08_dose import kur_table, build_tissue, KEEP_FACTOR
        edges, shares = kur_table()
        nuclides, order = build_tissue(KEEP_FACTOR * edges[-1])
        # 10/23 晚加：build_tissue 不带 "heat"，这里给每个核素补上（和 build_water 同一个读法）
        from lesson08_dose import T
        from lesson11_mt301 import load_heating
        for n in nuclides:
            E_h, heat = load_heating(n["name"], T)
            assert len(E_h) == len(n["E"]) and np.array_equal(E_h, n["E"]), "%s：heat 的能量网格和截面的不是同一张" % n["name"]
            n["heat"] = heat
        for k in REAL_CASES:
            inputs[k] = (nuclides, order, (edges, shares))
    else:
        print("  （没有 h5py 或真库 → 只跑假库 F1–F3，R1、R2 跳过；Windows 上 Ctrl+F5 会全跑）")
    return inputs


def to_plain(x):
    """把 numpy 数换成 Python float / list，好存 json（float 的 repr 存进去再读出来是逐位一样的）。"""
    if isinstance(x, dict):
        return {k: to_plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [to_plain(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    return x


def diffs(ref, new, path=""):
    """逐位比，返回所有不同之处（列表；空 = 完全相同）。只用 ==，不用 isclose。"""
    out = []
    if isinstance(ref, dict):
        if not isinstance(new, dict):
            return [(path, "类型不同", type(ref).__name__, type(new).__name__)]
        for k in sorted(set(ref) | set(new)):
            if k not in new:
                out.append((path + "/" + k, "新版少了这个键", "", ""))
            elif k not in ref:
                out.append((path + "/" + k, "新版多了这个键", "", ""))
            else:
                out += diffs(ref[k], new[k], path + "/" + k)
        return out
    if isinstance(ref, list):
        if len(ref) != len(new):
            return [(path, "长度不同", len(ref), len(new))]
        for i in range(len(ref)):
            out += diffs(ref[i], new[i], "%s[%d]" % (path, i))
        return out
    if not (ref == new):
        out.append((path, "值不同", repr(ref), repr(new)))
    return out


def report(name, ds):
    if not ds:
        return
    keys = sorted(set(p.split("/")[1].split("[")[0] for p, *_ in ds))
    print("  ✗ %s：%d 处不同，涉及的键：%s" % (name, len(ds), ", ".join(keys)))
    for p, what, a, b in ds[:6]:
        print("      %-28s %-8s 基准 %s  新版 %s" % (p, what, a, b))
    if len(ds) > 6:
        print("      ……（只列前 6 处）")


def md5_of(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def ref_path(name):
    return os.path.join(REF_DIR, "split_ref_%s.json" % name)


def try_fire(mod, inputs, ref_F2):
    """S3 试响：Q_val 相对改大 3e-12，只应改动 Bc 那两个键。"""
    if not hasattr(mod, "Q_val"):
        print("  S3 跳过：%s 里没有模块常数 Q_val（被挪走了？告诉 Claude 换一个试响点）" % NEW_MODULE)
        return
    keep = mod.Q_val
    try:
        mod.Q_val = keep * (1.0 + 3e-12)
        new = to_plain(run_case(mod, "F2", inputs["F2"]))
        ds = diffs(ref_F2, {k: new[k] for k in new if k in ref_F2})   # 新加的键不判（同 C1）
    finally:
        mod.Q_val = keep
    keys = sorted(set(p.split("/")[1].split("[")[0] for p, *_ in ds))
    assert keys == ["Bc_sq", "Bc_sum"], \
        "S3 试响失败：Q_val 只改了 3e-12，应恰好 Bc_sum、Bc_sq 两个键报差，实际报差的键是 %r —— 比较函数本身有问题，先别信 C1" % (keys,)
    print("  S3 试响过：Q_val × (1 + 3e-12) → 恰好 Bc_sum、Bc_sq 报差（%d 处），其余 %d 个键逐位同" % (len(ds), len(ref_F2) - 2))


if __name__ == "__main__":
    mod = importlib.import_module(NEW_MODULE)
    print("被测：%s.run_slab（%s）" % (NEW_MODULE, os.path.abspath(mod.__file__)))
    inputs = prepare_inputs()
    names = [k for k in list(FAKE_CASES) + list(REAL_CASES) if k in inputs]

    if not os.path.exists(ref_path("F1")):
        # ================= 存基准 =================
        print("\n【存基准】outputs\\split_ref\\ 里还没有基准 → 用现在的 run_slab 存")
        got = md5_of(BASE_FILE)
        assert got == BASE_MD5, \
            "S1：lesson07_freegas.py 的 md5 = %s，不是 10/22 收工那版 %s —— 你已经改过它了。" \
            "基准必须用拆之前的版本存：先告诉 Claude，别继续" % (got, BASE_MD5)
        print("  S1 过：lesson07_freegas.py 还是拆之前那版（md5 %s…）" % got[:8])
        os.makedirs(REF_DIR, exist_ok=True)
        refs = {}
        for name in names:
            refs[name] = to_plain(run_case(mod, name, inputs[name]))
            with open(ref_path(name), "w", encoding="utf-8") as f:
                json.dump(refs[name], f)
            print("  存好 %s：%d 个键（frac_T=%.4f  frac_R=%.4f  frac_A=%.4f  nt_cut=%d  n_ki=%d）"
                  % (name, len(refs[name]), refs[name]["frac_T"], refs[name]["frac_R"],
                     refs[name]["frac_A"], refs[name]["nt_cut"], refs[name]["n_ki"]))
        for name in FAKE_CASES:
            with open(ref_path(name), "r", encoding="utf-8") as f:
                back = json.load(f)
            ds = diffs(back, to_plain(run_case(mod, name, inputs[name])))
            report(name, ds)
            assert not ds, "S2：%s 同一个种子跑第二遍（并且经过存盘、读回）和第一遍不逐位相同 —— 先告诉 Claude" % name
        print("  S2 过：假库三个情况同种子再跑一遍，存盘读回后逐位相同")
        try_fire(mod, inputs, refs["F2"])
        print("\n基准存好了（%d 个情况）。现在可以开始拆 run_slab；每拆完一块，再 Ctrl+F5 跑这个文件。" % len(names))
    else:
        # ================= 检查 =================
        print("\n【检查】和 outputs\\split_ref\\ 里的基准逐位比")
        bad = []
        for name in names:
            p = ref_path(name)
            if not os.path.exists(p):
                print("  %s 没有基准（存基准那次没跑到它），跳过 —— 不会当成过" % name)
                bad.append(name + "（无基准）")
                continue
            with open(p, "r", encoding="utf-8") as f:
                ref = json.load(f)
            new = to_plain(run_case(mod, name, inputs[name]))
            extra = sorted(set(new) - set(ref))
            if extra:                                   # 10/23 晚：拆分做完，M2 起返回字典会加新键（kerma 等）—— 只比老的 28 个，新键列出、不判
                print("  %s 新加的键（不在基准里，不判）：%s" % (name, ", ".join(extra)))
            ds = diffs(ref, {k: new[k] for k in new if k in ref})
            if ds:
                report(name, ds)
                bad.append(name)
            else:
                print("  %s 逐位相同（%d 个键）" % (name, len(ref)))
        assert not bad, "C1 没过：%s。先看上面打印的是哪个键、第几个元素不同" % "、".join(bad)
        with open(ref_path("F2"), "r", encoding="utf-8") as f:
            try_fire(mod, inputs, json.load(f))
        print("\nC1 过：%d 个情况全部逐位相同 —— 这次拆分没改算法。" % len(names))
