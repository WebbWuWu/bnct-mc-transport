# -*- coding: utf-8 -*-
# ============================================================
#  profile_water.py  --  加速的第一步：量。水模体每秒跑多少条历史、时间花在哪几个函数
#  2026-11-11 Claude 写（11/10 定「加速第一步是量，Claude 给计时脚本，用户跑 5 分钟」）。只量，不改程序。
#
#  怎么跑：Windows、openmc 环境，VS Code 里打开本文件 Ctrl+F5（不是 F5：调试器会让程序慢好几倍，量出来的数没法用）。
#          跑的时候别同时跑别的重活。约 3–5 分钟。结果同时存到 outputs\profile_water.txt，发给 Claude。
#
#  🔴 跑之前先在纸上写下两个预测（新规矩 ①）：
#       1. 每秒大约多少条历史？（线索：M2 那次 10 万条用了 519 秒）
#       2. 下面「自己花的时间」那张表，排第一的是哪个函数？它占百分之几？
#
#  量三样：
#    A. 干净计时（不开剖析器）：三种配置各跑 N_TIME 条 → 每秒历史数、每条历史平均碰撞数
#         ① M2 配置（都不开） ② M3 配置（S(α,β) + 角分布） ③ M3 + 光子（gamma=True；光子主循环没写完时只多了生光子那一步）
#    B. 剖析（cProfile）：M3 配置跑 N_PROF 条 → 每个函数「自己花的时间」（tottime）和「连它调用的一起」（cumtime）
#    C. 调用次数：lookup_xs / fly / collide 各被调了多少次、平均一次多少微秒
#
#  怎么读 B：
#    tottime = 这个函数自己那几行花的时间，不算它调用的别的函数 → 找「哪几行最烫」看这一列
#    cumtime = 从进这个函数到出来的总时间 → 找「哪一块最烫」看这一列
#    ⚠️ 开了剖析器以后，每次函数调用都要多记一笔，「调用次数很多的小函数」会显得比实际更贵。所以每秒历史数以 A 为准，B 只看比例和排名。
# ============================================================
import os
import sys
import time
import io
import cProfile
import pstats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                         # my_code 根目录：本文件在 my_code\explore\ 里
sys.path.insert(0, ROOT)

import transport
from transport import run_slab
from water import build_water
from dose import kur_table, KEEP_FACTOR, D, NL, OUT_DIR
from thermal import load_thermal, SAB_PATH, T_TAG
from elastic_angle import load_elastic_angle, LIB, NUCS

N_TIME = 300          # A：每种配置跑几条（约 2 秒一百条的话，三种配置共约半分钟到一分钟）
N_PROF = 300          # B：剖析跑几条（开剖析器会慢 2–3 倍）
SEED = 4711           # 不取 1、0 这种特殊值；和以前的种子都不同，只用来计时
TOP = 18              # 表里列前几名


class Tee:
    """print 的东西同时打到屏幕和一个字符串里（最后存成 txt）。"""
    def __init__(self):
        self.buf = io.StringIO()

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s)
        self.buf.write(s + "\n")


def count_calls(names):
    """给 transport 里的几个函数各套一层「数次数 + 累计时间」的壳。返回 (计数字典, 还原函数)。
    只在 A 段用：壳本身每次调用多花约 0.2 微秒，比这几个函数自己小得多。"""
    stat = {n: [0, 0.0] for n in names}
    orig = {n: getattr(transport, n) for n in names}

    def wrap(n):
        f = orig[n]
        rec = stat[n]

        def g(*a, **k):
            t0 = time.perf_counter()
            out = f(*a, **k)
            rec[1] += time.perf_counter() - t0
            rec[0] += 1
            return out
        return g
    for n in names:
        setattr(transport, n, wrap(n))

    def restore():
        for n in names:
            setattr(transport, n, orig[n])
    return stat, restore


if __name__ == "__main__":
    out = Tee()
    edges, shares = kur_table()
    nuclides, order = build_water(KEEP_FACTOR * edges[-1])
    thermal = load_thermal(SAB_PATH, T_TAG)
    angle = load_elastic_angle(LIB, NUCS)
    base = dict(z_src=0.0, implicit=True, NL=NL, src=(edges, shares))
    configs = [
        ("① M2（都不开）", dict(thermal=None, angle=None, gamma=False)),
        ("② M3（S(α,β) + 角分布）", dict(thermal=thermal, angle=angle, gamma=False)),
        ("③ M3 + 光子（gamma=True）", dict(thermal=thermal, angle=angle, gamma=True)),
    ]

    out("=" * 78)
    out("A. 干净计时：水模体 15 cm、30 层、KUR 源、隐式俘获，每种配置 %d 条，种子 %d" % (N_TIME, SEED))
    out("=" * 78)
    out("%-28s %10s %12s %14s %16s" % ("配置", "用时 s", "历史 / 秒", "碰撞 / 历史", "10 万条要多久"))
    rate = {}
    for name, kw in configs:
        t0 = time.perf_counter()
        r = run_slab(N_TIME, SEED, D, nuclides, None, order, **base, **kw)
        dt = time.perf_counter() - t0
        rate[name] = N_TIME / dt
        out("%-28s %10.2f %12.1f %14.1f %13.1f 分钟" % (name, dt, N_TIME / dt, r["coll_per_hist"], 1e5 / (N_TIME / dt) / 60))
    out("")

    out("C. 三个大函数各被调了多少次（M3 配置，%d 条；套了计数的壳，所以总用时比上面略长）" % N_TIME)
    names = ["lookup_xs", "fly", "collide", "start_history", "end_history"]
    stat, restore = count_calls(names)
    try:
        t0 = time.perf_counter()
        run_slab(N_TIME, SEED, D, nuclides, None, order, **base, **configs[1][1])
        dt_all = time.perf_counter() - t0
    finally:
        restore()
    out("%-16s %12s %14s %14s %10s" % ("函数", "调用次数", "累计时间 s", "平均一次 µs", "占总用时"))
    acc = 0.0
    for n in names:
        c, t = stat[n]
        acc += t
        out("%-16s %12d %14.3f %14.1f %9.1f%%" % (n, c, t, (t / c * 1e6 if c else 0.0), t / dt_all * 100))
    out("%-16s %12s %14.3f %14s %9.1f%%   ← run_slab 自己那个循环、jump、断言等" % ("（其余）", "", dt_all - acc, "", (dt_all - acc) / dt_all * 100))
    out("")

    out("=" * 78)
    out("B. 剖析（cProfile）：M3 配置，%d 条。每秒历史数以 A 为准，这里只看比例和排名" % N_PROF)
    out("=" * 78)
    pr = cProfile.Profile()
    pr.enable()
    run_slab(N_PROF, SEED, D, nuclides, None, order, **base, **configs[1][1])
    pr.disable()
    for key, title in (("tottime", "自己花的时间（tottime）前 %d 名 —— 找「哪几行最烫」" % TOP),
                       ("cumtime", "连调用的一起（cumtime）前 %d 名 —— 找「哪一块最烫」" % TOP)):
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s)
        ps.strip_dirs().sort_stats(key).print_stats(TOP)
        out("--- " + title)
        txt = s.getvalue()
        out(txt[txt.find("ncalls"):] if "ncalls" in txt else txt)

    path = os.path.join(OUT_DIR, "profile_water.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write(out.buf.getvalue())
    print("存好：", path, "→ 把这个文件（或者屏幕输出）发给 Claude；先对一下你开头写的两个预测")
