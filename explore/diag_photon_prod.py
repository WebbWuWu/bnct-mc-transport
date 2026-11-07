# explore/diag_photon_prod.py —— 11/7 第 3 项：lib80x 的中子文件里有没有「光子产物」（Claude 写，不默写）
# 看 H1、O16、B10 每个反应道下面的 product_*，particle 属性是不是 photon
# 运行：Windows 的 openmc (3.10.20) 环境，Ctrl+F5
import os
import h5py

LIB = r"E:\BNCT_GPU\lib80x_hdf5"

for nuc in ["H1", "O16", "B10"]:
    with h5py.File(os.path.join(LIB, nuc + ".h5"), "r") as f:
        rx = f[nuc]["reactions"]
        print("=" * 20, nuc, "：共 %d 个反应道" % len(rx), "=" * 20)
        n_ph = 0
        for name in rx:
            g = rx[name]
            parts = []
            for key in g:
                if key.startswith("product_"):
                    p = g[key].attrs.get("particle", b"?")
                    p = p.decode() if isinstance(p, bytes) else str(p)
                    parts.append(p)
            if "photon" in parts:
                n_ph += 1
                print("  %-14s MT=%-4s 产物: %s" % (name, g.attrs.get("mt", "?"), parts))
        print("  → 带光子产物的反应道：%d 个" % n_ph)
