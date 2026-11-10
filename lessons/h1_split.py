# ---- 路径（2026-11-10 整理文件夹，Claude 加）：本文件在 my_code 的子文件夹里。
#      下面几行把 my_code 根目录和各子文件夹加进搜索路径，后面的 import 才找得到别的文件。----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

from xs_lookup import load_nuclide, find_cell, xs_at

e = 0.0246375                                   # H-1 294K 网格第 242 个点（《库结构说明》§7 的 2.464e-02）
E, mts, XS, awr, Q = load_nuclide("H1", "294K")    # 4 个返回值，顺序和 xs_lookup 的 __main__ 一样
lo = find_cell(E, e)                            # 格子只找一次，所有道共用这个 lo

total = 0.0                                     # 累加器：数的是全部道，所以放在 for 外面
for j in range(len(XS)):                        # 走一遍所有叶子道
    s = xs_at(E, XS[j], e, lo)                  # 第 j 个道在 e 处的截面（barn）
    print("MT =", mts[j], "  截面 =", s)
    total += s
print("加起来 =", total)
