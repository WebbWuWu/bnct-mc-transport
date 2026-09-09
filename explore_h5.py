import os
import h5py

PATH = r"E:\BNCT_GPU\lib80x_hdf5\H1.h5"

# assert 挡在它要挡的那个操作前面 —— 路径错了就一行说清楚，
# 不是让 h5py 吐 10 行内部调用栈。
assert os.path.exists(PATH), "文件不在这里：%s" % PATH

with h5py.File(PATH, "r") as f:
    print("顶层：", list(f.keys()))
    print("H1 下面：", list(f["H1"].keys()))
    print("温度档：", list(f["H1"]["energy"].keys()))
    E = f["H1"]["energy"]["294K"][:]
    print("294K 网格 %d 点，%.3e ~ %.3e eV" % (len(E), E[0], E[-1]))