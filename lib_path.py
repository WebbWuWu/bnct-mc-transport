# -*- coding: utf-8 -*-
# ============================================================
#  lib_path.py  --  截面库在哪（整个程序只在这一处写）
#  2026-11-10 Claude 建（整理文件夹时）。以前这个路径在 xs_lookup / thermal / elastic_angle / photon 里各写了一遍。
#
#  库：LANL 版 ENDF/B-VIII.0 `lib80x_hdf5`（约 15 GB，不在仓库里）。
#  换电脑时二选一：① 改下面那一行的路径；② 不改文件，设环境变量 BNCT_LIB = 你的库目录。
# ============================================================
import os

LIB_DIR = os.environ.get("BNCT_LIB") or (r"E:\BNCT_GPU\lib80x_hdf5" if os.name == "nt" else "/mnt/e/BNCT_GPU/lib80x_hdf5")   # Windows / WSL 各一个默认值
