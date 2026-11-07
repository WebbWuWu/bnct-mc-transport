# explore/diag_angle_interp.py —— 11/7 角分布诊断第 1 步（Claude 写，不默写、不进默写范围）
# 看 O16、H1 弹性（MT=2）角分布 mu 数据集的 interpolation 属性（1 = 直方图，2 = 线性）
# 第二条路：按「线性」「直方图」两种理解分别对 pdf 积分，哪种 = 1 就是哪种
# 运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5
import os
import numpy as np
import h5py

LIB = r"E:\BNCT_GPU\lib80x_hdf5"

for nuc in ["H1", "O16"]:
    with h5py.File(os.path.join(LIB, nuc + ".h5"), "r") as f:
        dist = f[nuc]["reactions"]["reaction_002"]["product_0/distribution_0"]
        ang = dist["angle"]
        E = ang["energy"][:]
        mu_ds = ang["mu"]
        print("=" * 20, nuc, "=" * 20)
        print("distribution_0 属性:", dict(dist.attrs))
        print("angle 组属性:", dict(ang.attrs))
        print("energy 数据集属性:", dict(ang["energy"].attrs))
        print("mu 数据集属性名:", list(mu_ds.attrs.keys()))
        if "interpolation" in mu_ds.attrs:
            itp = np.asarray(mu_ds.attrs["interpolation"])
            vals, cnts = np.unique(itp, return_counts=True)
            print("interpolation 取值: 个数 =", dict(zip(vals.tolist(), cnts.tolist())),
                  "（共 %d 个能量点）" % len(E))
        else:
            itp = None
            print("mu 上没有 interpolation 属性")

        mu_all = mu_ds[:]
        offs = mu_ds.attrs["offsets"]
        ends = np.append(offs[1:], mu_all.shape[1])
        print("%11s %5s %4s | %8s %8s | %8s %8s %8s | %8s"
              % ("E(eV)", "npts", "itp", "归一线性", "归一直方",
                 "<mu>线性", "<mu>直方", "<mu>cdf", "cdf末"))
        for E0 in [1.0e4, 1.0e5, 5.0e5, 1.0e6, 2.0e6, 5.0e6, 1.0e7]:
            j = int(np.argmin(abs(E - E0)))
            m, p, c = mu_all[0:3, offs[j]:ends[j]]
            dm = np.diff(m)
            n_lin = np.sum(dm * (p[:-1] + p[1:]) / 2.0)
            n_his = np.sum(dm * p[:-1])
            m_lin = np.sum(dm * (p[:-1] * (2*m[:-1] + m[1:]) + p[1:] * (m[:-1] + 2*m[1:])) / 6.0)
            m_his = np.sum(dm * p[:-1] * (m[:-1] + m[1:]) / 2.0)
            m_cdf = m[-1]*c[-1] - m[0]*c[0] - np.sum(dm * (c[:-1] + c[1:]) / 2.0)
            it = "-" if itp is None else str(int(itp[j]))
            print("%11.4e %5d %4s | %8.5f %8.5f | %8.5f %8.5f %8.5f | %8.5f"
                  % (E[j], len(m), it, n_lin, n_his, m_lin, m_his, m_cdf, c[-1]))
