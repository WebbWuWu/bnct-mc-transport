# -*- coding: utf-8 -*-
# ============================================================
#  lesson01_flight.py  --  第1课 任务B：抽样中子飞行距离
# ============================================================
#
# 【物理背景】
#   中子在均匀介质中，两次碰撞之间的飞行距离 s 服从指数分布：
#       p(s) = Sigma_t * exp(-Sigma_t * s)
#
#   用"逆变换法"抽样：
#       s = -ln(1 - xi) / Sigma_t        其中 xi 是 [0,1) 均匀随机数
#
#   Sigma_t（宏观总截面，单位 cm^-1）的物理意义：
#       中子每走 1 cm 发生一次反应的概率
#   平均自由程 lambda = 1 / Sigma_t
#       两次碰撞之间平均走多远

# 【语法点 16】import 把外部的工具库拿进来
#   math 是 Python 自带的数学库，之后可以用：
#       math.log(x)   自然对数 ln(x)
#       math.exp(x)   e 的 x 次方
#       math.pi       圆周率
import math


A = 1664525
C = 1013904223
M = 2**32


# 和任务A完全一样，自己再敲一遍，别复制
class LCG:
    def __init__(self, seed):
        self.state = seed

    def next_uint(self):
        self.state = (A*self.state+C)%M               # TODO 同任务A
        return self.state

    def random(self):
        return self.next_uint()/M                    # TODO 同任务A


# 【语法点 17】写在 class 外面的 def 就是普通函数
#   它需要什么就通过参数传进来（这里要传一个 rng 对象和一个 sigma_t 数值）
def sample_flight(rng, sigma_t):
    """抽一次飞行距离，单位 cm"""
    xi = rng.random()
    # TODO 填这里：写出逆变换公式 s = -ln(1 - xi) / sigma_t
    #   提示：自然对数用 math.log(...)
    return -math.log(1-xi)/sigma_t


# ============================================================
#  主程序
# ============================================================
SIGMA_T = 1.6658                       # 水的宏观总截面，单位 cm^-1
N = 1000000

rng = LCG(seed=2026)

# 【语法点 18】空列表 [] 和 append
#   samples = [] 造一个空列表
#   samples.append(x) 在末尾追加一个元素
samples = []
for i in range(N):
    samples.append(sample_flight(rng, SIGMA_T))


# ---- 检验 1：样本均值应等于 1/Sigma_t ----
# 【语法点 19】sum() 和 len() 是 Python 内置函数
#   sum(列表) 求和，len(列表) 求元素个数
mean = sum(samples) / N
print("样本均值    =", mean, "cm")
print("理论 1/St   =", 1.0 / SIGMA_T, "cm")
print("相对偏差    =", (mean - 1.0/SIGMA_T) / (1.0/SIGMA_T) * 100, "%")


# ---- 检验 2：中位数应等于 ln2/Sigma_t ----
# 【语法点 20】列表方法 .sort() 原地从小到大排序
#   注意它没有返回值，是直接把 samples 本身改掉
samples.sort()
median = samples[N // 2]               # 排好序后取正中间那个
print("\n样本中位数  =", median, "cm")
print("理论 ln2/St =", math.log(2.0) / SIGMA_T, "cm")


# ---- 检验 3：未碰撞穿透率，应等于 exp(-Sigma_t * d) ----
# 【语法点 21】for 也可以直接遍历一个写死的列表
#   d 依次取 0.5, 1.0, 2.0, 3.0
for d in [0.5, 1.0, 2.0, 3.0]:
    count = 0
    for s in samples:
        # 【语法点 22】if 判断
        #   条件成立才执行缩进的那一行
        #   比较符：>  <  >=  <=  ==（相等）  !=（不等）
        #   注意两个等号 == 才是判断相等，一个等号 = 是赋值
        if s > d:                      # 飞行距离大于 d，说明它穿过了 d 还没碰撞
            count += 1
    print("d =", d, "cm   MC =", count / N,
          "   理论 =", math.exp(-SIGMA_T * d))
