# -*- coding: utf-8 -*-
# ============================================================
#  check_p2.py  --  M4 判据 P2：Klein–Nishina 抽 μ（photon.sample_kn_mu）对不对
#  Claude 2026-11-21 写（判据归 Claude，成品直接跑）
#
#  用法：VS Code 打开本文件，Ctrl+F5（openmc (3.10.20) 环境）。不读截面库，约十几秒。
#        先跑 check_P2（正式判），再跑 check_P2_selftest（6 种故意写错的必须响，纪律 4）。
#
#  被测：photon.sample_kn_mu(E, rng) —— E 单位 eV，返回 μ = cos θ
#
#  三个能量：2.2233 MeV（H 俘获 γ）、0.4778 MeV（B10 的 γ）、10 keV（接近 Thomson 极限）
#  每个能量抽 N = 100000 次，4 条：
#    [P2-1] 每个 μ 都在 [−1, 1] 里、不是 None                         抓：μ 映射错（只抽 [0,1)）、忘了 return
#    [P2-2] ⟨μ⟩ vs 参考值，z = 差 / √(Var μ / N)                       抓：KN 公式写错（漏项、k 写反）、单位错
#    [P2-3] 往前散（μ ≥ 0）的比例 vs 参考值，二项 SE                   抓：同上；⟨μ⟩ 碰巧对上时还能抓形状
#    [P2-4] 平均每次试几轮 vs 1/ε，ε = σ_KN / (4π rₑ²)                 抓：舍选上界写错（用 f 而不是 f/2）、
#           轮数 = 这 N 次一共抽的随机数个数 ÷ 2（rng.n 计数器）          接受概率 > 1（公式漏项时 f 会超过 2）
#           ⚠️ 默认每试一轮正好抽 2 个数（一个抽 μ、一个舍选）
#
#  两条路（纪律 2）：
#    你的路   = 舍选抽样（μ 均匀 + 以 f/2 接受）
#    判据的路 = ⟨μ⟩、往前比例：KN 公式数值积分（Claude 11/21 用 200 万点梯形积分算，
#               2.2233、0.4778 MeV 两行和第 13 课 §3.2 的表一致：0.435 / 77%、0.284 / 69%）
#               ε：Heitler 的 KN 总截面闭式解（另一种写法，不经过 f(μ)），本文件里现算
#    先确认比的是同一个量：σ_KN = π rₑ² ∫ f dμ，所以舍选效率 = ∫ f dμ / (2 × 2) = σ_KN / (4π rₑ²)
#  门槛在看数据之前定（11/21 23:58）：|z| < 4（和 check_emit 一样；3 个能量 × 3 条 z = 9 次检验，
#  全对时误报的概率约 0.06%）
#  盲区：只看 ⟨μ⟩、前后比例、效率三个数，不看完整的直方图；E′ 的公式不在这里（在 track_photon 里）
# ============================================================

# ---- 路径：本文件在 my_code\checks\，下面几行把 my_code 根目录加进搜索路径 ----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import math
from rng import lcg
from photon import sample_kn_mu, ME_C2

N = 100000
N_SELF = 5000        # 自检用：写错的都错得很明显，5000 次足够响（11/21 试过）
Z_GATE = 4.0

# 能量（eV）→ (⟨μ⟩, ⟨μ²⟩, 往前比例)：KN 数值积分（见文件头）
REF = {
    2.2233e6: (0.43522, 0.50449, 0.77367),
    0.4778e6: (0.28414, 0.44144, 0.68834),
    1.0e4:    (0.01517, 0.40011, 0.51066),
}


def kn_efficiency(E):
    """舍选效率 ε = σ_KN / (4π rₑ²)，σ_KN 用 Heitler 闭式解（不经过 f(μ)，第二条路）。"""
    a = E / ME_C2
    L = math.log(1 + 2 * a)
    s = 2 * ((1 + a) / a**2 * (2 * (1 + a) / (1 + 2 * a) - L / a)
             + L / (2 * a) - (1 + 3 * a) / (1 + 2 * a)**2)       # σ_KN / (π rₑ²)
    return s / 4


def run_one(fn, E, n, seed):
    """对一个能量跑 n 次，返回 (说明列表, 是否全过)。"""
    m1, m2, fw = REF[E]
    eps = kn_efficiency(E)
    rng = lcg(seed)
    n0 = rng.n
    s1 = 0.0
    nf = 0
    bad = 0
    for _ in range(n):
        mu = fn(E, rng)
        if mu is None or not (-1.0 <= mu <= 1.0):
            bad += 1
            if bad >= 20:                      # 已经判定不过，不用再跑完（写错的版本可能很慢）
                return ["    [P2-1] μ 越界或 None 已有 %d 个，提前停  ❌" % bad], False
            continue
        s1 += mu
        if mu >= 0.0:
            nf += 1
    good = n - bad
    lines = []
    ok_all = True

    ok1 = (bad == 0)
    lines.append("    [P2-1] μ 都在 [−1, 1]：越界或 None %d 个  %s" % (bad, "✅" if ok1 else "❌"))
    ok_all &= ok1
    if good == 0:
        return lines, False

    mean = s1 / good
    se_m = math.sqrt((m2 - m1 * m1) / good)
    z2 = (mean - m1) / se_m
    ok2 = abs(z2) < Z_GATE
    lines.append("    [P2-2] ⟨μ⟩ = %.5f vs %.5f  z = %+.2f  %s" % (mean, m1, z2, "✅" if ok2 else "❌"))
    ok_all &= ok2

    pf = nf / good
    z3 = (pf - fw) / math.sqrt(fw * (1 - fw) / good)
    ok3 = abs(z3) < Z_GATE
    lines.append("    [P2-3] 往前比例 = %.5f vs %.5f  z = %+.2f  %s" % (pf, fw, z3, "✅" if ok3 else "❌"))
    ok_all &= ok3

    tries = (rng.n - n0) / 2.0 / n
    se_t = math.sqrt(1 - eps) / eps / math.sqrt(n)
    z4 = (tries - 1 / eps) / se_t
    ok4 = abs(z4) < Z_GATE
    lines.append("    [P2-4] 平均每次试 %.4f 轮 vs 1/ε = %.4f（ε = %.4f）  z = %+.2f  %s"
                 % (tries, 1 / eps, eps, z4, "✅" if ok4 else "❌"))
    ok_all &= ok4
    return lines, ok_all


def check_P2(fn=sample_kn_mu, n=N, seed=730241, verbose=True):
    ok_all = True
    for i, E in enumerate(REF):
        lines, ok = run_one(fn, E, n, seed + 1000 * i)
        if verbose:
            print("  E = %.4g eV（α = %.4f），N = %d" % (E, E / ME_C2, n))
            for ln in lines:
                print(ln)
        ok_all &= ok
    return ok_all


# ---------------- 自检：故意写错的版本（纪律 4），每个都必须响 ----------------
MAX_TRIES = 1000      # 写错的版本可能几乎永远不接受（W6 的效率约 10⁻⁶）：试满 1000 轮还没接受就返回 None，P2-1 会响


def _kn_template(E, rng, f_of, accept_scale=0.5, map_mu=lambda x: 2 * x - 1, me=ME_C2, k_inv=False):
    alpha = E / me
    for _ in range(MAX_TRIES):
        mu = map_mu(rng.random())
        k = 1 + alpha * (1 - mu) if k_inv else 1 / (1 + alpha * (1 - mu))
        f = f_of(k, mu)
        if rng.random() < f * accept_scale:
            return mu
    return None


_F_OK = lambda k, mu: k**2 * (k + 1 / k - 1 + mu**2)
WRONG = [
    ("W1 漏了「−1」（11/21 你犯的）", lambda E, r: _kn_template(E, r, lambda k, mu: k**2 * (k + 1 / k + mu**2))),
    ("W2 接受概率用 f 而不是 f/2",     lambda E, r: _kn_template(E, r, _F_OK, accept_scale=1.0)),
    ("W3 k 写反成 1 + α(1 − μ)",       lambda E, r: _kn_template(E, r, _F_OK, k_inv=True)),
    ("W4 只用 Thomson (1+μ²)/2",       lambda E, r: _kn_template(E, r, lambda k, mu: 1 + mu**2)),
    ("W5 μ 只在 [0, 1) 里抽",          lambda E, r: _kn_template(E, r, _F_OK, map_mu=lambda x: x)),
    ("W6 mₑc² 用了 MeV（0.511）",      lambda E, r: _kn_template(E, r, _F_OK, me=0.511)),
]


def check_P2_selftest():
    all_fired = True
    for name, fn in WRONG:
        fired = not check_P2(fn, n=N_SELF, seed=55117, verbose=False)
        print("  %-28s %s" % (name, "✅ 响了" if fired else "❌ 没响"))
        all_fired &= fired
    # 对照：模板配上正确的 f，必须不响（证明模板本身没写错）
    quiet = check_P2(lambda E, r: _kn_template(E, r, _F_OK), n=N_SELF, seed=55117, verbose=False)
    print("  %-28s %s" % ("对照：模板 + 正确公式", "✅ 不响" if quiet else "❌ 误响"))
    return all_fired and quiet


if __name__ == "__main__":
    import time
    t0 = time.time()
    print("P2 正式判：photon.sample_kn_mu")
    ok = check_P2()
    print("  → %s\n" % ("✅ P2 通过" if ok else "❌ P2 没过"))
    print("P2 自检：6 个故意写错的版本必须响")
    ok_s = check_P2_selftest()
    print("  → %s" % ("✅ 自检通过" if ok_s else "❌ 自检没过"))
    print("\n用时 %.1f 秒" % (time.time() - t0))
