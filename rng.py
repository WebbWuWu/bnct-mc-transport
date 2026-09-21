# -*- coding: utf-8 -*-
# ============================================================
#  rng.py  --  线性同余随机数发生器 lcg + 各向同性方向抽样
#  2026-09-19 抽成独立文件，lesson03_tally 和 lesson06_continuous 共用
# ============================================================
#
# 🔴 模块名 rng 和代码里的局部变量 rng 同名（rng = lcg(seed)）。
#    一律这样导入：  from rng import lcg, isotropic_direction
#    （现有的 lesson03 / lesson06 都是这么写的，所以现在不会出事。）
#
#    如果写成  import rng  再写  rng = rng.lcg(seed)：
#      · 赋值那一行不报错 —— 模块名 rng 被悄悄换成了一个 lcg 对象；
#      · 之后再写 rng.isotropic_direction(...) 才报 AttributeError，
#        报错的位置离真正出错的那一行很远。
# ============================================================

import math

a=1664525
c=1013904223
m=2**32

class lcg:
    def __init__(self,seed):
        self.state=seed
    def uint(self):
        self.state=(self.state*a+c)%m
        return self.state
    def random(self):
        return self.uint()/m


def isotropic_direction(rng):
    mu=rng.random()*2.0-1.0
    phi=2*rng.random()*math.pi
    sin_theta=math.sqrt(1-mu**2)
    u=sin_theta*math.cos(phi)
    v=sin_theta*math.sin(phi)
    w=mu
    return u,v,w