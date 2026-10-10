# -*- coding: utf-8 -*-
# ============================================================
#  whiteboard/check_fgk.py  --  free_gas_kinematics 默写稿 vs 原文：同种子逐位对拍
#  Claude 2026-10-10 写。用法：把下面 MOD 改成你默写稿的文件名（不带 .py），Ctrl+F5。
# ============================================================
#
#  两条路：lesson07_freegas.free_gas_kinematics（原文，10/3 判据全过的那份）
#          vs  你的默写稿里同名函数
#  同一组输入、同一个种子 → 四个返回值必须【逐位相同】（==，不是「差不多」），
#  而且两边用掉的随机数个数也要一样（看发生器跑完后的 state）。
#  输入取 8 组：H-1 和 O-16 两种 A，热区 / 超热 / keV 三档能量，方向都不取 (0,0,1)（纪律 5）。
#  抓：任何一步写错（点乘写成各乘各、A 和 A+1、少一步归一化、u/v/w 抄错一个分量、
#      多吃 / 少吃一个随机数）—— 一个比特都对不上就响。
# ============================================================

import os
import sys
import math

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))      # my_code 根目录：rng、lesson07_freegas 在那里
sys.path.insert(0, HERE)                       # whiteboard：默写稿在这里

MOD = "rewrite_fgk1"                           # ← 改成你的默写稿文件名

from rng import lcg
from lesson07_freegas import free_gas_kinematics as fgk_ref, sample_target
mine = __import__(MOD)
fgk_mine = mine.free_gas_kinematics

KT = 0.025301
CASES = []
for A in (0.999167, 15.857510):
    for E in (0.0371, 3.73, 2.9e3):
        CASES.append((E, A))
CASES.append((0.0123, 0.999167))
CASES.append((47.3, 15.857510))

worst = 0
for n, (E, A) in enumerate(CASES):
    s = 1.0 / math.sqrt(0.37**2 + 0.48**2 + 0.79**2)
    u, v, w = 0.37 * s, -0.48 * s, 0.79 * s
    # 先用同一个发生器抽靶核（两边共用这一步，不是被测的部分）
    r0 = lcg(9100 + n)
    y = math.sqrt(E * A / KT)
    x, mu_T, _ = sample_target(y, r0)
    seed_state = r0.state
    ra = lcg(0); ra.state = seed_state
    rb = lcg(0); rb.state = seed_state
    a = fgk_ref(E, u, v, w, A, KT, x, mu_T, ra)
    b = fgk_mine(E, u, v, w, A, KT, x, mu_T, rb)
    same = all(p == q for p, q in zip(a, b)) and len(a) == len(b)
    print("第 %d 组  E=%-8g A=%-9g  原文 E'=%.15g  默写 E'=%.15g  %s"
          % (n, E, A, a[0], b[0], "逐位相同" if same else "✗"))
    # [Claude assert] 抓：默写稿任何一步和原文不同 —— 四个返回值不能有一个比特不同
    assert same, "第 %d 组（E=%r, A=%r）对不上：\n  原文 %r\n  默写 %r" % (n, E, A, a, b)
    # [Claude assert] 抓：多吃 / 少吃了随机数（比如重抽 r 的循环写错）—— 两边发生器最后停在同一个位置
    assert ra.state == rb.state, "第 %d 组：原文用完随机数后 state=%d，默写 state=%d —— 吃的随机数个数不同" \
        % (n, ra.state, rb.state)

print("8 组全部逐位相同，随机数用量一致 ✓")
