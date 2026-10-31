# -*- coding: utf-8 -*-
# ============================================================
#  check_rng.py  --  随机数换法乙（模 2^63 + 每条历史跳到第 i*S 个数）的判据
#  判据归 Claude。2026-10-31 写。被测：rng.py 的 lcg / jump（你写）、
#  lesson07_freegas.run_slab 历史循环开头那一行起点（你写）。
#
#  运行：Windows，Ctrl+F5。只用 check_split 的假库，不读真截面，约 10–20 秒。
#  全过的最后一行是「5 条全过」；任何一条不过，assert 停在那一条，消息说哪里不对。
# ============================================================
#
#  验证（5 条；每条一行「抓什么」）
#    R1 大发生器上，jump(x0, n) 和「一步一步走 n 步」逐位相同（n = 0、7、12345、3*S+17）
#       抓：拼接 / 自拼公式写错、自拼两行顺序反、忘了 n = n>>1、% 写成 /
#       试响（纪律 4）：本文件里放了一个「漏乘 A」的错版 jump（B_ans + B_cur），R1 必须对它响。
#       🔴 要在大模数、n = 12345 上试：小发生器 (5,3,16) 里走 4 步那档 A = 1，漏乘也算对（10/30 讲义订正）
#    R2 lcg(9137) 调 777 次 uint() 后的状态 == jump(9137, 777)
#       抓：lcg 类里的常数和 jump 的默认参数不是同一组（两条独立的路：类一步步走 vs 跳跃）
#    R3 random()：20 万个数全在 [0,1)、均值离 0.5 不超过 4 个 SE；构造「下一状态 = m−1」时仍 < 1
#       抓：取高 53 位那一行写错（移位位数、除数、运算顺序）—— ξ 跑出 [0,1) 或挤到一端
#    R4 run_slab（假库 F1 模拟、F2 隐式）里，每条历史用掉的 uint 次数 < S = 152917
#       抓：一条历史用的数挤进下一条历史的段 —— 10/28「序列重叠」同一种病
#    R5 run_slab 每条历史开头设起点：第 i 条 = jump(seed, i*S)，共 N 次；第一次设起点之前没用过随机数；
#       同一种子跑两遍返回字典逐位相同
#       抓：起点用了 rng.state（接着上一条跳）、漏乘 S、设在循环外 / start_history 之后、有不受种子控制的东西
# ============================================================

import math
import importlib

import rng as rng_mod
from rng import lcg, jump, S, a, c, m
from check_split import fake_nuclides, FAKE_CASES, FAKE_SRC   # 只借假库，不跑它的主程序

slab = importlib.import_module("lesson07_freegas")


# ------------------------------------------------------------
#  R1  跳跃 = 一步一步走
# ------------------------------------------------------------
def walk(x0, n):
    x = x0
    for _ in range(n):
        x = (x * a + c) % m
    return x


def jump_bad(x0, n, a=a, c=c, m=m):
    """试响用的错版：拼进答案时漏乘 A_cur（B_ans + B_cur）。"""
    A_ans, B_ans, A_cur, B_cur = 1, 0, a, c
    while n > 0:
        if n & 1:
            A_ans = (A_cur * A_ans) % m
            B_ans = (B_ans + B_cur) % m          # ← 故意错
        B_cur = ((A_cur + 1) * B_cur) % m
        A_cur = (A_cur * A_cur) % m
        n = n >> 1
    return (A_ans * x0 + B_ans) % m


X0 = 20261031
for n in [0, 7, 12345, 3 * S + 17]:
    want = walk(X0, n)
    got = jump(X0, n)
    assert got == want, "R1 ✗ n=%d：jump 得 %d，一步一步走得 %d —— 跳跃函数和逐步走不一致" % (n, got, want)
assert jump_bad(X0, 12345) != walk(X0, 12345), \
    "R1 试响 ✗：漏乘 A 的错版在 n=12345 上居然也对上了 —— R1 抓不住这种错，判据没分辨率"
print("R1 ✓ 跳跃 = 一步一步走（n = 0、7、12345、3S+17）；漏乘 A 的错版会响")


# ------------------------------------------------------------
#  R2  lcg 类和 jump 用同一组常数
# ------------------------------------------------------------
g = lcg(9137)
for _ in range(777):
    g.uint()
assert g.state == jump(9137, 777), \
    "R2 ✗：lcg(9137) 走 777 步的状态 %d ≠ jump(9137,777) = %d —— 类里的常数和 jump 默认参数不是同一组" \
    % (g.state, jump(9137, 777))
print("R2 ✓ lcg 类走 777 步 = jump(9137, 777)")


# ------------------------------------------------------------
#  R3  random() 的区间
# ------------------------------------------------------------
NR = 200000
g = lcg(4471)
lo, hi, tot = 1.0, 0.0, 0.0
for _ in range(NR):
    xi = g.random()
    lo = min(lo, xi)
    hi = max(hi, xi)
    tot += xi
mean = tot / NR
se = math.sqrt(1.0 / 12.0 / NR)                  # 均匀分布 σ² = 1/12，SE = σ/√N
assert 0.0 <= lo and hi < 1.0, "R3 ✗：random() 跑出了 [0,1)：最小 %r，最大 %r" % (lo, hi)
assert abs(mean - 0.5) < 4 * se, "R3 ✗：random() 均值 %.6f 离 0.5 超过 4 个 SE（SE = %.2e）—— 分布挤到一端" % (mean, se)
g = lcg(0)
g.state = ((m - 1 - c) * pow(a, -1, m)) % m      # 构造：下一次 uint() 正好得 m−1（最大的状态）
xi_top = g.random()
assert xi_top < 1.0, "R3 ✗：状态 = m−1 时 random() 返回 %r，取到了 1 —— -log(1-ξ) 会炸" % xi_top
print("R3 ✓ random()：20 万个数在 [%.2e, %.12f]，均值 %.5f（z = %.2f）；最大状态时 = %r"
      % (lo, hi, mean, (mean - 0.5) / se, xi_top))


# ------------------------------------------------------------
#  R4、R5  run_slab 里的用量和起点（从外面替换 uint、jump，你的程序不加计数）
# ------------------------------------------------------------
log = {"calls": [], "used": [], "since": 0, "before_first": None}
orig_uint = lcg.uint
orig_jump = slab.jump


def counting_uint(self):
    log["since"] += 1
    return orig_uint(self)


def recording_jump(x0, n, *args, **kw):
    if not log["calls"]:
        log["before_first"] = log["since"]       # 第一次设起点之前用了几个数（应为 0）
    else:
        log["used"].append(log["since"])         # 上一条历史用了几个数
    log["since"] = 0
    log["calls"].append((x0, n))
    return orig_jump(x0, n, *args, **kw)


def run_fake(name, nuclides, order):
    cs = FAKE_CASES[name]
    return slab.run_slab(cs["N"], cs["seed"], cs["d"], nuclides, cs["E_src"], order, **cs["kw"])


nuclides, order = fake_nuclides()
results = {}
lcg.uint = counting_uint
slab.jump = recording_jump
try:
    for name in ["F1", "F2"]:
        log.update(calls=[], used=[], since=0, before_first=None)
        results[name] = run_fake(name, nuclides, order)
        log["used"].append(log["since"])          # 最后一条历史
        cs = FAKE_CASES[name]
        N, seed = cs["N"], cs["seed"]

        # R5 起点
        assert len(log["calls"]) == N, \
            "R5 ✗ %s：设起点 %d 次，应等于历史数 N = %d —— 那一行不在历史循环开头" % (name, len(log["calls"]), N)
        assert log["before_first"] == 0, \
            "R5 ✗ %s：第一次设起点之前已经用了 %d 个随机数 —— 设起点那行要放在 start_history 之前" % (name, log["before_first"])
        for i, (x0, n) in enumerate(log["calls"]):
            assert x0 == seed and n == i * S, \
                "R5 ✗ %s：第 %d 条历史的起点是 jump(%d, %d)，应为 jump(seed=%d, i*S=%d)" % (name, i, x0, n, seed, i * S)

        # R4 用量
        worst = max(log["used"])
        assert worst < S, \
            "R4 ✗ %s：有一条历史用了 %d 个随机数，超过步长 S = %d —— 会挤进下一条历史的段" % (name, worst, S)
        print("R4 ✓ %s：%d 条历史，每条最多用 %d 个随机数（S = %d，占 %.3f%%）"
              % (name, N, worst, S, 100.0 * worst / S))
finally:
    lcg.uint = orig_uint
    slab.jump = orig_jump

again = run_fake("F1", nuclides, order)
assert again == results["F1"], "R5 ✗：F1 同一种子跑两遍，返回字典不逐位相同 —— 程序里有不受种子控制的东西"
print("R5 ✓ F1、F2：第 i 条历史起点 = jump(seed, i*S)，共 N 次，之前没用过随机数；F1 同种子两遍逐位相同")

print("\n5 条全过")
