# -*- coding: utf-8 -*-
# ============================================================
#  check_load_kT.py —— 默写 load_kT + sigma_at_ext 的对拍（Claude 10/16 写）
#  运行：Windows，Ctrl+F5（要 h5py，读真库）。几秒。
#  比什么：你的 whiteboard\rewrite_load_kT1.py 里的两个函数 vs 原函数（lesson07_freegas.py），逐位相同才算过
#    ① load_kT：H1、O16 的 294K（两个核素，防止「只对 H1」）
#    ② sigma_at_ext：第 6 课假库，3 个能量 —— 低于表下限（外推）、正好在表下限、表里面（不外推）
#    ③ sigma_at_ext：真库 H1 在表下限以下 3.1e-6 eV（真表的道数 > 2，抓「只写对两个道」）
# ============================================================
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))      # my_code 根目录：lesson07_freegas、xs_lookup 在那里
sys.path.insert(0, HERE)                       # whiteboard：默写稿在这里

MOD = "rewrite_load_kT1"                       # 你的默写文件名（不带 .py）
from lesson07_freegas import load_kT as kT_ref, sigma_at_ext as ext_ref
from xs_lookup import load_nuclide
mine = __import__(MOD)

ok = True
print("① load_kT")
for nuc in ["H1", "O16"]:
    a, b = mine.load_kT(nuc, "294K"), kT_ref(nuc, "294K")
    same = (a == b)
    ok = ok and same
    print("  %-4s 你的 %r   原函数 %r   %s" % (nuc, a, b, "逐位相同" if same else "✗ 不同"))

print("② sigma_at_ext，假库 E = [1e-5, 2e7]，两个道 0.8 / 0.2")
E_f = [1e-5, 2e7]
XS_f = [[0.8, 0.8], [0.2, 0.2]]
for e in [3.1e-6, 1e-5, 0.37]:
    a, b = mine.sigma_at_ext(E_f, XS_f, e), ext_ref(E_f, XS_f, e)
    same = (list(a[0]) == list(b[0]) and a[1] == b[1])
    ok = ok and same
    print("  e = %-8g 你的 %r\n             原函数 %r   %s" % (e, a, b, "逐位相同" if same else "✗ 不同"))

print("③ sigma_at_ext，真库 H1 294K，e = 3.1e-6 eV（表下限以下）")
E, mts, XS, awr, Q = load_nuclide("H1", "294K")
a, b = mine.sigma_at_ext(E, XS, 3.1e-6), ext_ref(E, XS, 3.1e-6)
same = (list(a[0]) == list(b[0]) and a[1] == b[1])
ok = ok and same
print("  表下限 %g eV，%d 个道；总截面 你的 %r  原函数 %r   %s" % (E[0], len(XS), a[1], b[1], "逐位相同" if same else "✗ 不同"))

print("=" * 60)
print("对拍全过：load_kT + sigma_at_ext 默写逐位相同" if ok else "有 ✗ —— 对着上面那行看是哪个函数、哪个能量")
