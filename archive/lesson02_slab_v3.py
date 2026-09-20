# -*- coding: utf-8 -*-
# ============================================================
#  lesson02_slab_v3.py  --  无提示重写（考试模式）
# ============================================================
#
#  开写之前先做这件事：**关掉所有自动补全**
#    1. 状态栏右下角的 Copilot 图标 → Disable Completions
#    2. Ctrl+Shift+P → "Preferences: Open User Settings (JSON)" → 加上：
#         "editor.inlineSuggest.enabled": false,
#         "editor.quickSuggestions": false,
#         "editor.suggestOnTriggerCharacters": false,
#         "editor.parameterHints.enabled": false
#    写完 v3 再改回 true —— 这是考试模式，不是长期设置。
#
#  规则：从这里往下全部自己写。
#  不看 lesson02_slab.py，不看 lesson02_slab_v2.py，不看第 1 课。
#  卡住了就停下来想；想不出来就在下面记一行「卡在哪」，不要去翻旧文件。
#
# ------------------------------------------------------------
#  要实现的东西
# ------------------------------------------------------------
#  一块厚度 d 的均匀平板，z 从 0 到 d，x/y 方向无限大。
#  中子从 z=0 处垂直射入（初始方向 +z）。
#  统计三个互斥且穷尽的结局：透射 T、反射 R、吸收 A。
#
# ------------------------------------------------------------
#  你必须自己写出来的部件
# ------------------------------------------------------------
#  1. LCG 随机数发生器
#     乘数 1664525，增量 1013904223，模 2**32。
#     要有一个返回 [0,1) 浮点数的方法。
#
#  2. 各向同性方向抽样，返回单位向量 (u, v, w)
#
#  3. 单条历史的循环
#     抽飞行距离 → 推进 → 出界判断 → 碰撞后定散射/吸收 → 散射则换方向继续
#     加一个最大事件数上限防死循环
#
#  4. 平均碰撞次数统计
#     定义：每条历史发生的碰撞总数（**包括最后那次导致吸收的碰撞**）除以 N。
#     想清楚两件事：
#       - 计数器在哪一层循环里初始化、在哪一行自增
#       - 历史有几种结束方式？每一种都把这条历史的账记进总数了吗？
#
#  5. assert 自检（v2 完全没写，这次必须有）
#     把下面「自检标准」里能写成 assert 的都写成 assert，不要只 print。
#
# ------------------------------------------------------------
#  三个 case（N = 1000000，seed = 2026）
# ------------------------------------------------------------
#     case 1:  d=2.0   Sigma_t=1.0   Sigma_s=0.0
#     case 2:  d=2.0   Sigma_t=1.0   Sigma_s=0.8
#     case 3:  d=5.0   Sigma_t=1.0   Sigma_s=0.9
#
#  三个 case 都要跑。只跑 case 1 就下结论是 v2 踩过的坑。
#
# ------------------------------------------------------------
#  自检标准
# ------------------------------------------------------------
#  - 三个 case 的 T+R+A 都必须等于 1（用 assert，容差 1e-12）
#  - case 1 的 R 必须严格等于 0
#  - case 1 的 T 对上 exp(-2) = 0.135335，差值在 3 倍标准误内
#  - case 1 的平均碰撞次数 = 1 - exp(-2) = 0.8647
#  - case 2 平均碰撞约 2.24，case 3 约 5.17
#  - 触顶最大事件数的历史数必须是 0
#
# ------------------------------------------------------------
#  写完之后
# ------------------------------------------------------------
#  和 v2 做 diff。差异的地方就是昨天 Copilot 替你想的部分——
#  逐条看清楚你自己写的和它写的差在哪，这份 diff 是今天最有价值的产物。
#
# ------------------------------------------------------------
#  卡住的地方记在这里（不许翻旧文件，只许记）：
#
#
#
# ============================================================
#  以下开始写你的代码
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

abb=lcg(2026)
for i in range(10):
    print("i=%.6f"%abb.random())

def isotropic_direction(rng):
    mu=rng.random()*2.0-1.0
    phi=2*rng.random()*math.pi
    sin_theta=math.sqrt(1-mu**2)
    u=sin_theta*math.cos(phi)
    v=sin_theta*math.sin(phi)
    w=mu
    return u,v,w

Max=10000

def run_slab(N,seed,d,sigma_t,sigma_s):
    rng=lcg(seed)
    n_T=0
    n_R=0
    n_A=0
    n_stuck=0
    n_i=0
    for i in range(N):
        z=0.0
        w=1.0
        for step in range(Max):
            # 为什么是 1-rng.random() 而不是 rng.random()：
            #   random() 返回 [0,1) —— 0 取得到，1 取不到。若写 -log(xi)，xi=0 时 log 收到 0 → -inf，程序炸。
            #   取 1-xi 后 log 的参数落在 (0,1]，永远取不到 0。两种写法同分布，但安全性完全不同。
            #   命中概率约 2^-32：跑一亿条历史大概撞一次 —— 调试时永远不出现，演示时出现。
            #   这个选择配套的是本文件里 LCG 的区间，换 RNG 要重新核对。见 whiteboard_02_spec.md。
            s=-math.log(1-rng.random())/sigma_t
            z=z+s*w
            if z<0:
                n_R+=1
                break
            if z>d:
                n_T+=1
                break
            if rng.random()<sigma_s/sigma_t:
                u,v,w=isotropic_direction(rng)
                n_i+=1
            else:
                n_A+=1
                n_i+=1
                break
        else:
            n_stuck+=1
    if n_stuck>0:
        print("warning!有%d条历史触顶max_event"% n_stuck)
    return n_T/N,n_R/N,n_A/N,n_i/N

N=100000
print("case1:纯吸收：d=2.0,sigma_t=1.0,sgma_s=0.0")
T,R,Ab,TA=run_slab(N,2026,2.0,1.0,0.0)
print("T=%.6f,R=%.6f,A=%.6f,T+R+A=%.6f,Ta=%.6f"%(T,R,Ab,T+R+Ab,TA))
assert abs(T + R + Ab - 1.0) < 1e-12, \
    "守恒破了：T+R+A=%.15f" % (T + R + Ab)

assert R == 0.0, \
    "纯吸收不可能反射，实际 R=%.6f" % R

p = math.exp(-2.0)
se = (p * (1 - p) / N) ** 0.5
assert abs(T - p) < 3 * se, \
    "T=%.6f 偏离解析解 %.6f 达 %.2f 个标准误" % (T, p, (T - p) / se)

print("  case 1 自检通过")

print("case2 各向同性散射：d=2.0,sigma_t=1.0,sigma_s=0.8")
Ta,Ra,Ac,TB=run_slab(N,2026,2.0,1.0,0.8)
print("T=%.6f,R=%.6f,A=%.6f,T+R+A=%.6f,Ta=%.6f"%(Ta,Ra,Ac,Ta+Ra+Ac,TB))
assert abs(Ta+Ra+Ac - 1.0) < 1e-12, \
    "守恒破了：T+R+A=%.15f" % (Ta+Ra+Ac)
    
assert Ra != 0.0, \
    "各向异性散射一定有反射，实际 R=%.6f" %Ra
    

        
print("case3 各向同性散射：d=5.0,sigma_t=1.0,sigma_s=0.9")
Tb,Rb,Ad,TC=run_slab(N,2026,5.0,1.0,0.9)
print("T=%.6f,R=%.6f,A=%.6f,T+R+A=%.6f,Ta=%.6f"%(Tb,Rb,Ad,Tb+Rb+Ad,TC))
assert abs(Tb+Rb+Ad - 1.0) < 1e-12, \
    "守恒破了：T+R+A=%.15f" % (Tb+Rb+Ad)
    
assert Rb != 0.0, \
    "各向异性散射一定有反射，实际 R=%.6f" %Rb