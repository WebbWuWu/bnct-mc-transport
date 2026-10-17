# -*- coding: utf-8 -*-
# ============================================================
#  whiteboard/copy_fgk_v2.py  --  free_gas_kinematics 精简版（抄写稿）
#  Claude 2026-10-17 写（用户 17:26：「重复内容太多，用 (vnx,vny,vnz) 这种表示，帮我改精简」）
#  和 lesson07_freegas.py 里的原文【逐位相同】（同输入、同种子，四个返回值 ==，随机数用量一样）：
#      Claude 在容器里用 2 万组随机输入对拍过；你可以用 whiteboard\check_fgk.py（MOD = "copy_fgk_v2"）自己再跑一次
#  lesson07_freegas.py 本身还没动 —— 换不换、谁来换，你定（见对话）
# ============================================================
#
#  精简的办法只有一个：三个分量不再各写一行，而是装进一个长度 3 的列表，一行算完三个
#     原文                                    精简版
#     vnx=y*u / vny=y*v / vnz=y*w      →      vn = [y*d[i] for i in range(3)]
#     sqrt(rx**2+ry**2+rz**2)（写了 4 次）→    norm(r)（写一次函数，到处用）
#     uu*u+vv*v+ww*w                   →      dot(p, d)
#
#  新名字 ↔ 原文名字（意思一个没变，讲义 §4.3 的四步照旧）
#     d    = (u, v, w)            中子方向
#     vn   = (vnx, vny, vnz)      中子速度（以 v_T 为尺）
#     p    = (uu, vv, ww)         各向同性抽的辅助方向
#     r    = (rx, ry, rz)         p 去掉沿 d 的部分
#     e1   = (e1x, e1y, e1z)      r 归一化
#     rb   = (rbx, rby, rbz)      靶核方向
#     vb   = (vbx, vby, vbz)      靶核速度
#     vcm  = (vcmx, vcmy, vcmz)   质心速度
#     vc   = (vncmx, vncmy, vncmz) CM 里的中子速度；s = vncm（它的速率）
#     q    = (u1, v1, w1)         CM 里抽的新方向
#     vl   = (vlx, vly, vlz)      LAB 新速度；L = vl（它的速率）
#
#  你没学过的用法（规则 13，例子和项目无关）
#   ① 列表推导：  [表达式 for i in range(3)]  = 把 i 取 0、1、2 各算一次，排成列表
#        price = [12.0, 8.5, 3.0]
#        with_tax = [price[i] * 1.1 for i in range(3)]      # [13.2, 9.35, 3.3]（第 i 个只和第 i 个算）
#   ② 函数返回的三个数可以只用一个名字接，再用 [0] [1] [2] 取
#        def size(): return 3, 4, 5
#        s = size();  s[1]                                  # 4
#      （a, b, c = size() 也行 —— 两种写法拿到的是同样三个数）
#   ③ 一个小函数到处用，比同一个式子抄四遍好记，也不会某一遍抄错
#
#  抄写规矩（10/14）：assert、print 不抄，看一眼就行
# ============================================================

import math
from rng import isotropic_direction


def dot(a, b):
    """两个三维向量的点乘：对应分量相乘再加起来（不是各乘各）。"""
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]


def norm(a):
    """三维向量的长度。"""
    return math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)


def free_gas_kinematics(E, u, v, w, A, kT, x, mu_T, rng):
    """四步矢量变换（讲义 §4.3）。输入输出和原文完全一样：返回 E', u', v', w'。"""
    d = (u, v, w)                                         # 中子方向
    y = math.sqrt(E*A/kT)
    # 1.  中子速度 = y × 中子方向
    vn = [y*d[i] for i in range(3)]
    # 2a. 抽一个方向 p，减掉它沿 d 的部分 → r；剩下太短就重抽（最多 100 轮）
    for k in range(100):
        p = isotropic_direction(rng)
        al = dot(p, d)
        r = [p[i] - d[i]*al for i in range(3)]
        if norm(r) >= 1e-6:
            break
    else:
        # [Claude assert] for 的 else —— 抓：剩余长度永远不够（正常撞满概率约 1e-1230）
        assert False, "重抽 r 满 100 轮（正常不可能，必是 bug）：E=%r A=%r kT=%r y=%r" % (E, A, kT, y)
    # 2b. r 归一化 → e1（和 d 垂直的单位向量）
    nr = norm(r)
    e1 = [r[i]/nr for i in range(3)]
    # 2c. 靶核方向 = mu_T × d + sqrt(1 − mu_T²) × e1
    st = math.sqrt(1 - mu_T**2)
    rb = [mu_T*d[i] + st*e1[i] for i in range(3)]
    # [Claude assert] 2c 之后 —— 抓：e1 没和中子方向垂直 / 没归一化、某个分量算错
    assert abs(dot(rb, rb) - 1.0) < 1e-9 and abs(dot(rb, d) - mu_T) < 1e-9, \
        "靶核方向不对：长度²=%.12g（应为 1），与中子方向夹角余弦=%.12g（应为 mu_T=%.12g）" % (dot(rb, rb), dot(rb, d), mu_T)
    # 2d. 靶核速度 = x × 靶核方向
    vb = [x*rb[i] for i in range(3)]
    # 3.  质心速度 = (中子速度 + A × 靶核速度) / (1 + A)
    vcm = [(vn[i] + A*vb[i])/(1 + A) for i in range(3)]
    # 4.  CM 里中子速度 = 中子速度 − 质心速度，记下速率 s
    vc = [vn[i] - vcm[i] for i in range(3)]
    s = norm(vc)
    # 5+6. CM 里各向同性抽新方向 q，乘同一个速率 s，加回质心速度 → LAB 新速度 vl
    q = isotropic_direction(rng)
    vl = [vcm[i] + q[i]*s for i in range(3)]
    L = norm(vl)
    # [Claude assert] 除以 L 之前 —— 抓：LAB 新速度为 0，方向无定义（④续 甲）
    assert L > 0.0, "LAB 新速度为 0，方向无定义：E=%r A=%r kT=%r x=%r mu_T=%r" % (E, A, kT, x, mu_T)
    # E' = (kT/A) × |vl|²；新方向 = vl / |vl|
    return (kT/A)*L**2, vl[0]/L, vl[1]/L, vl[2]/L
