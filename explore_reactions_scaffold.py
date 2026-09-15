# -*- coding: utf-8 -*-
# ============================================================
#  explore_reactions_scaffold.py  --  第 5 课 · 第 1 步
#  把一个核素的反应道摊开看。2026-09-15
# ============================================================
#
#  【这一步只打印，不计算总截面。】
#
#  为什么先做这一步：9/9 你发现"文件里没有总截面"。它不是数据缺了，
#  是 OpenMC 的 HDF5 格式【故意不存】MT=1，因为它是冗余的。它要你自己加。
#  但"加哪些"这件事，不看一眼数据是定不下来的。
#
#  9/8 那条纪律在这里第二次生效：**数字要打出来，而且要有人看。**
#
#  你要填 5 个空，标着 ① ② ③ ④ ⑤。
#  **其余一行都不要动。**（9/8 lesson03 就是因为没说这句，被顺手删掉 114 行）
#
#  跑法：VS Code 里打开本文件 → 确认焦点在编辑器（不是"输出"面板）→ F5
# ============================================================

import os
import h5py
import numpy as np

LIB    = r"E:\BNCT_GPU\lib80x_hdf5"
NUC    = "B10"          # 先跑 B10（59 个道，坑最多）；跑通了把这里改成 "H1" 再跑一遍
T      = "294K"
E_TEST = 0.0253         # eV，热中子的"标准"能量。B-10 的 (n,alpha) 在这里是 3844 barn

path = os.path.join(LIB, NUC + ".h5")
# assert 挡在它要挡的那个操作前面 —— 路径错了就一行说清楚。
assert os.path.exists(path), "文件不在这里：%s" % path


with h5py.File(path, "r") as f:
    root = f[NUC]

    # --------------------------------------------------------
    #  能量网格
    # --------------------------------------------------------
    E = root["energy"][T][:]
    print("%s 的 %s 能量网格：%d 点，%.3e ~ %.3e eV" % (NUC, T, len(E), E[0], E[-1]))

    # --------------------------------------------------------
    #  定位热能点 —— 顺便认识 numpy 的第一个函数
    # --------------------------------------------------------
    # np.searchsorted(E, x) 回答："x 插到哪个位置，能保持 E 仍然升序？"
    #
    # side="right" 这个参数不是可有可无的：
    #   side="left"  遇到 x 正好等于某个网格点时，插在它【前面】
    #   side="right" 遇到 x 正好等于某个网格点时，插在它【后面】
    # 你自己写的 find 在 e == E[mid] 时 return mid，等价于 side="right" 再减 1。
    # 用错一边，就是你 9/7 反复踩的"差一格"，而且只在 x 恰好落在网格点上时才现形。
    i0 = int(np.searchsorted(E, E_TEST, side="right")) - 1
    print("热能点落在第 %d 格：E[%d]=%.6e <= %.4f < E[%d]=%.6e"
          % (i0, i0, E[i0], E_TEST, i0 + 1, E[i0 + 1]))

    # ---- 判据 B：numpy 的答案必须和你自己写的二分一样 ----
    #   这条抓"我对 searchsorted 的理解错了"，也顺带再测一次你的 find。
    #   两条【独立写出来】的算路走到同一个下标，才说明两边都对。
    from whiteboard_05_rewrite import find
    k_mine = find(list(E), E_TEST)
    assert k_mine == i0, "你的 find 给 %s，numpy 给 %d —— 必须一致" % (k_mine, i0)
    print("判据B：自写二分与 numpy.searchsorted 一致，都是 %d" % i0)
    print()

    # --------------------------------------------------------
    #  ① 取出 reactions 组下面【所有反应道的名字】
    # --------------------------------------------------------
    #  提示：h5py 的 group 用起来像 dict。
    #        root["reactions"] 是一个 group，
    #        list(某个 group.keys()) 给你一个名字列表，
    #        形如 ['reaction_002', 'reaction_102', 'reaction_800', ...]
    #
    names = list(root["reactions"].keys())

    print("%s 一共有 %d 个反应道" % (NUC, len(names)))
    print()
    print("  道名            %-14s  MT  redundant  thr_idx  len(xs)   %.4f eV 处"
      % ("反应", E_TEST))
    print("  " + "-" * 74)

    n_redundant = 0
    n_threshold = 0

    for name in sorted(names):
        rx = root["reactions"][name]

        # ----------------------------------------------------
        #  ② 从反应道的【属性】里取出 mt 和 redundant
        # ----------------------------------------------------
        #  提示：属性【不在】 keys() 里。它挂在 .attrs 上，.attrs 也像 dict：
        #        某个东西.attrs["属性名"]
        #  要取的两个属性名就叫 "mt" 和 "redundant"。
        #
        #  （如果这里报 KeyError，说明这个版本的库没存该属性 —— 停下来告诉我，别猜。）
        mt        = rx.attrs["mt"]
        redundant = rx.attrs["redundant"]
        label = rx.attrs["label"]

        # ----------------------------------------------------
        #  ③ 取出这个道在 294K 下的 xs 数据集，以及它的 threshold_idx
        # ----------------------------------------------------
        #  xs 的路径是   rx[T]["xs"]
        #
        #  ⚠️ 先把【数据集对象本身】拿到手，**不要急着加 [:]**。
        #     因为 threshold_idx 是挂在【数据集】的 .attrs 上的，
        #     而 [:] 一取就变成普通 numpy 数组，属性就丢了。
        #     —— 这就是"取数据"和"取属性"是两件事。
        xs_ds = rx[T]["xs"]
        thr   = int(xs_ds.attrs["threshold_idx"])
        xs    = xs_ds[:]        # 这一行已经写好，别动

        # ---- 判据 A：xs 数组 + 阈值偏移，必须正好铺满能量网格 ----
        #   这条抓"我以为所有道都躺在同一张网格上"。
        #   阈值反应（比如 (n,2n)）的 xs 只从第 thr 个点开始存，前面的零不存，
        #   所以 len(xs) 比 len(E) 短，短掉的正好是 thr。
        #
        #   ④ 把这个等式写出来。
        assert len(E)-len(xs)==thr, \
            "%s: len(xs)=%d, thr=%d, len(E)=%d —— 对不上" % (name, len(xs), thr, len(E))

        # ----------------------------------------------------
        #  ⑤ 取出这个道在 E_TEST 处的截面值
        # ----------------------------------------------------
        #  这是坑一的具体形态，也是今天唯一需要你动脑的一个空：
        #
        #      xs[0] 对应的【不是】 E[0]，而是 E[thr]。
        #
        #  那么，能量网格上的下标 i0，对应 xs 里的哪个下标？
        #  在纸上写两行再填。填错了不会报错，只会给你一个安静的错数字。
        if i0 < thr:
            val = 0.0          # 还没到阈值，这个道在这个能量下根本不发生
        else:
            val = xs[i0-thr]

        if redundant:
            n_redundant += 1
        if thr > 0:
            n_threshold += 1

        print("  %-16s %-14s %4d   %6d   %7d  %7d   %14.6g"
        % (name, label, mt, redundant, thr, len(xs), val))

    print("  " + "-" * 74)
    print("  合计 %d 个道：其中 redundant=1 的 %d 个，有阈值的 %d 个"
          % (len(names), n_redundant, n_threshold))
    print()
    print("跑完了先别写求和。看着这张表回答三个问题（写在纸上）：")
    print("  1. redundant=1 的那几个道，MT 号分别是多少？它们各自是谁的和？")
    print("  2. 如果把它们也加进去，B10 在 0.0253 eV 的总截面会变成多少？")
    print("  3. 有阈值的道，在 0.0253 eV 处贡献是多少？为什么？")
