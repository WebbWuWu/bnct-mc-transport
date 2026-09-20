# -*- coding: utf-8 -*-
"""
第 6 课开工前的一次性测量：插值方式会不会破坏「total == sum(各道)」。

背景（原理讲义 第6课 §4）：
    线性插值对数据是【线性的】，所以它和求和可以交换次序；
    log-log 插值不是线性的，不能交换。
    而 9/16 白板题⑦ 写下的判据是  assert abs(sum(sigma)/total - 1) < 1e-12 。
    如果 total 和各道分别做 log-log 插值，这条判据会在一半以上的能量点炸，
    而代码没有任何错 —— 真因是用两条不同的路算了同一个量。

本脚本只测量，不下结论。跑完自己看数字。
Claude 2026-09-17 在 B10 上跑过一次，结果见讲义 §4.1。
"""
import os
import numpy as np
import h5py

LIB = r"E:\BNCT_GPU\lib80x_hdf5"
NUCLIDE = "B10"
T = "294K"


def lin(Ea, Eb, Sa, Sb, e):
    """线性插值（lin-lin）"""
    return Sa + (Sb - Sa) * (e - Ea) / (Eb - Ea)


def loglog(Ea, Eb, Sa, Sb, e):
    """双对数插值（白板题⑥ 的幂形式）。端点有 0 时退回线性。"""
    if Sa <= 0 or Sb <= 0:
        return lin(Ea, Eb, Sa, Sb, e)
    p = np.log(Sb / Sa) / np.log(Eb / Ea)
    return Sa * (e / Ea) ** p


def load(nuclide=NUCLIDE, T=T):
    """返回 (E, total, leaves)；leaves = [(mt, 对齐到全网格的 xs 数组), ...]"""
    with h5py.File(os.path.join(LIB, nuclide + ".h5"), "r") as f:
        root = f[nuclide]
        E = np.array(root["energy"][T])
        total = np.zeros(len(E))
        leaves = []
        for name in root["reactions"]:
            rx = root["reactions"][name]
            # 9/15 纪律：读 redundant 属性，绝不硬编码 MT 列表
            if int(rx.attrs["redundant"]) == 1:
                continue
            xs = rx[T]["xs"]
            thr = int(xs.attrs["threshold_idx"])   # 阈值反应从阈值点才开始存
            a = np.zeros(len(E))
            a[thr:thr + len(xs)] = np.array(xs)
            total += a
            leaves.append((int(rx.attrs["mt"]), a))
    return E, total, leaves


def main():
    E, total, leaves = load()
    print("%s %s: 网格 %d 点, 非冗余道 %d, E = %.3e ~ %.3e eV"
          % (NUCLIDE, T, len(E), len(leaves), E[0], E[-1]))

    # ---- 测 1: 总截面在格子几何中点上，两种插值差多少 ----
    rel = []
    for k in range(len(E) - 1):
        if E[k + 1] <= E[k] or total[k] <= 0 or total[k + 1] <= 0:
            continue
        e = np.sqrt(E[k] * E[k + 1])
        a = lin(E[k], E[k + 1], total[k], total[k + 1], e)
        b = loglog(E[k], E[k + 1], total[k], total[k + 1], e)
        rel.append(abs(b / a - 1))
    rel = np.array(rel)
    print("\n[1] 总截面格子中点  log-log vs lin-lin 相对差")
    print("    中位 %.3e  平均 %.3e  最大 %.3e   超过 1e-3 的格子 %d/%d"
          % (np.median(rel), rel.mean(), rel.max(), (rel > 1e-3).sum(), len(rel)))

    # ---- 测 2: 先插各道再求和  vs  先求和再插 ----
    worst, worst_e, cnt, big = 0.0, 0.0, 0, 0
    for k in range(len(E) - 1):
        if E[k + 1] <= E[k] or total[k] <= 0 or total[k + 1] <= 0:
            continue
        e = np.sqrt(E[k] * E[k + 1])

        s_log = sum(loglog(E[k], E[k + 1], a[k], a[k + 1], e) for _, a in leaves)
        t_log = loglog(E[k], E[k + 1], total[k], total[k + 1], e)
        r = abs(s_log / t_log - 1)
        cnt += 1
        if r > 1e-6:
            big += 1
        if r > worst:
            worst, worst_e = r, e

        # lin-lin 必须逐格严格相等 —— 这条 assert 就是「线性可交换」的判据
        s_lin = sum(lin(E[k], E[k + 1], a[k], a[k + 1], e) for _, a in leaves)
        t_lin = lin(E[k], E[k + 1], total[k], total[k + 1], e)
        assert abs(s_lin - t_lin) <= 1e-9 * max(1.0, abs(t_lin)), \
            "lin-lin 竟然不可交换，k=%d, e=%.6e, 差 %.3e" % (k, e, abs(s_lin - t_lin))

    print("\n[2] log-log 下  先插各道再求和  vs  先求和再插总截面")
    print("    最大相对差 %.3e  在 E = %.4e eV;  超过 1e-6 的格子 %d/%d"
          % (worst, worst_e, big, cnt))
    print("    lin-lin :  %d 个格子全部通过严格相等的 assert" % cnt)

    print("\n判读提示：第 [2] 项才是会咬人的那个。")
    print("9/16 的判据是  abs(sum(sigma)/total - 1) < 1e-12 。")
    print("对照上面 log-log 那一行，自己判断这条判据能不能活下来。")


if __name__ == "__main__":
    main()
