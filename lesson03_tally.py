import mcstat
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

Max=10000

NL=20

def run_slab(N,seed,d,sigma_t,sigma_s):
    rng=lcg(seed)
    n_T=0
    n_R=0
    n_A=0
    n_stuck=0
    n_i=0
    dz=d/NL
    track_sum=[0.0]*NL
    track_sq=[0.0]*NL
    v_sum=v_sq=0.0
    for i in range(N):
        z=0.0
        w=1.0
        track_h=[0.0]*NL
        n_hist=0
        for step in range(Max):
            
            # 为什么是 1-rng.random() 而不是 rng.random()：
            #   random() 返回 [0,1) —— 0 取得到，1 取不到。若写 -log(xi)，xi=0 时 log 收到 0 → -inf，程序炸。
            #   取 1-xi 后 log 的参数落在 (0,1]，永远取不到 0。两种写法同分布，但安全性完全不同。
            #   命中概率约 2^-32：跑一亿条历史大概撞一次 —— 调试时永远不出现，演示时出现。
            #   这个选择配套的是本文件里 LCG 的区间，换 RNG 要重新核对。见 whiteboard_02_spec.md。
            s=-math.log(1-rng.random())/sigma_t
            z0=z
            z1=z0+s*w
            if z1>d:
                z_end=d
            elif z1<0:
                z_end=0.0
            else:
                z_end=z1
            assert 0.0 <= z0 <= d, "z0 跑出板外了: %.15f" % z0
            seg = 0.0                                  # 本段实际记了多少
            if w == 0.0:
                k0 = int(z0 / dz)
                assert 0 <= k0 < NL, "层号越界: k0=%d z0=%.9f" % (k0, z0)
                track_h[k0] += s
                seg    += s
                expect  = s

            else:
                k0 = int(z0 / dz)                      # 起点层
                k1 = int(z_end / dz)                   # 终点层
                if k1 >= NL:                           # z_end 恰好 == d 时夹一下
                    k1 = NL - 1

                cur = z0                               # 笔尖当前位置
                k   = k0                               # 笔尖当前层
                guard=0
                while True:
                    guard+=1
                    assert guard<=NL+2,\
                        "拆层循环未退出"
                    # ① 这一步走到哪：当前层的出口边界 vs z_end，取近的那个
                    #    w>0 出口是 (k+1)*dz ；w<0 出口是 k*dz
                    if w>0:
                        k_o=(k+1)*dz
                        z_step=min(k_o,z_end)
                    else:
                        k_o=k*dz
                        z_step=max(k_o,z_end)
                    # ② 记账：dz_step = abs(这一步终点 - cur)
                    #         track_h[k] += dz_step / abs(w)
                    #         seg        += 同一个值
                    #         记账前先 assert 0 <= k < NL
                    assert 0<=k<NL,\
                        "k越界"
                    dz_step=abs(z_step-cur)
                    track_h[k] += dz_step / abs(w)
                    seg+=dz_step/abs(w)
                    # ③ cur = 这一步的终点
                    cur=z_step
                    # ④ 到 k1 了就 break，否则 k += 1（w>0）或 k -= 1（w<0）
                    if k==k1:
                        break
                    if w > 0:
                        k += 1
                    else:
                        k -= 1

                expect = abs(z_end - z0) / abs(w)

            assert math.isclose(seg, expect, rel_tol=1e-9, abs_tol=1e-12), \
                "段记账不守恒: seg=%.15f expect=%.15f z0=%.6f z_end=%.6f w=%.9f" \
                % (seg, expect, z0, z_end, w)

            z = z1                                     
            if z<0:
                n_R+=1
                break
            if z>d:
                n_T+=1
                break
            if rng.random()<sigma_s/sigma_t:
                u,v,w=isotropic_direction(rng)
                n_i+=1
                n_hist+=1
            else:
                n_A+=1
                n_i+=1
                n_hist+=1
                break
        else:
            n_stuck+=1
        for j in range(NL):
            track_sum[j]+=track_h[j]
            track_sq[j]+=track_h[j]**2
        x=sigma_t*sum(track_h)-n_hist
        v_sum+=x
        v_sq+=x**2
    if n_stuck>0:
        print("warning!有%d条历史触顶max_event"% n_stuck)
    assert n_T+n_R+n_A+n_stuck == N
    scale = abs(sigma_t*sum(track_sum)) + abs(n_i)
    assert math.isclose(sigma_t*sum(track_sum) - n_i, v_sum,
                    rel_tol=0.0, abs_tol=1e-9*scale), \
                        "验证B总账不一致" 
    return n_T/N,n_R/N,n_A/N,n_i/N,track_sum,track_sq,dz,v_sum,v_sq,sigma_t

def print_flux_table(track_sum, track_sq, N, dz, NL, title):
    print()
    print(title)
    print(f"{'k':>3} {'z_lo':>7} {'z_hi':>7} {'phi':>12} {'abs_err':>11} {'R':>9}")
    print("-" * 54)
    for k in range(NL):
        z_lo = k * dz
        z_hi = z_lo + dz
        Phi,Var,Err=mcstat.mean_se(track_sum[k],track_sq[k],N)
        phi = Phi/dz          # 平均径迹 / 体积，你有公式
        err = Err/dz         # 标准误 / 体积
        R   = err / phi if phi > 0 else 0.0
        print(f"{k+1:>3} {z_lo:>7.3f} {z_hi:>7.3f} {phi:>12.6f} {err:>11.2e} {R:>9.5f}")
        
#test
N = 1000000

print("case1:纯吸收：d=2.0,sigma_t=1.0,sgma_s=0.0")
T,R,Ab,TA,sum1,sq1,dz1,v_sum1,v_sq1,sigma_t1=run_slab(N,2026,2.0,1.0,0.0)
print("T=%.6f,R=%.6f,A=%.6f,T+R+A=%.6f,Ta=%.6f"%(T,R,Ab,T+R+Ab,TA))
print_flux_table(sum1,sq1,N,dz1,NL,"case1:纯吸收：d=2.0,sigma_t=1.0,sgma_s=0.0")
assert abs(T + R + Ab - 1.0) < 1e-12, \
    "守恒破了：T+R+A=%.15f" % (T + R + Ab)
for j in range(NL):
    j1,j2,Err1=mcstat.mean_se(sum1[j],sq1[j],N)
    err1=Err1/dz1
    assert abs(sum1[j]/N/dz1-((math.exp(-sigma_t1*j*dz1)-math.exp(-sigma_t1*(dz1+j*dz1)))/(dz1*sigma_t1)))<3*err1,\
        "3sigma检验不合格"
mean1,var1,err_t1=mcstat.mean_se(v_sum1,v_sq1,N)
assert abs(mean1)<3*err_t1,\
    "case1碰撞次数检验不合格"

assert R == 0.0, \
    "纯吸收不可能反射，实际 R=%.6f" % R

p = math.exp(-2.0)
se = (p * (1 - p) / N) ** 0.5
assert abs(T - p) < 3 * se, \
    "T=%.6f 偏离解析解 %.6f 达 %.2f 个标准误" % (T, p, (T - p) / se)

print("  case 1 自检通过")

print("case2 各向同性散射：d=2.0,sigma_t=1.0,sigma_s=0.8")
Ta,Ra,Ac,TB,sum2,sq2,dz2,v_sum2,v_sq2,sigma_t2=run_slab(N,2026,2.0,1.0,0.8)
print("T=%.6f,R=%.6f,A=%.6f,T+R+A=%.6f,Ta=%.6f"%(Ta,Ra,Ac,Ta+Ra+Ac,TB))
print_flux_table(sum2,sq2,N,dz2,NL,"case2 各向同性散射：d=2.0,sigma_t=1.0,sigma_s=0.8")
assert abs(Ta+Ra+Ac - 1.0) < 1e-12, \
    "守恒破了：T+R+A=%.15f" % (Ta+Ra+Ac)
    
assert Ra != 0.0, \
    "各向异性散射一定有反射，实际 R=%.6f" %Ra

mean2,var2,err_t2=mcstat.mean_se(v_sum2,v_sq2,N)   
assert abs(mean2)<3*err_t2,\
    "case2碰撞次数检验不合格"
    

        
print("case3 各向同性散射：d=5.0,sigma_t=1.0,sigma_s=0.9")
Tb,Rb,Ad,TC,sum3,sq3,dz3,v_sum3,v_sq3,sigma_t3=run_slab(N,2026,5.0,1.0,0.9)
print("T=%.6f,R=%.6f,A=%.6f,T+R+A=%.6f,Ta=%.6f"%(Tb,Rb,Ad,Tb+Rb+Ad,TC))
print_flux_table(sum3,sq3,N,dz3,NL,"case3 各向同性散射：d=5.0,sigma_t=1.0,sigma_s=0.9")
assert abs(Tb+Rb+Ad - 1.0) < 1e-12, \
    "守恒破了：T+R+A=%.15f" % (Tb+Rb+Ad)
    
assert Rb != 0.0, \
    "各向异性散射一定有反射，实际 R=%.6f" %Rb

mean3,var3,err_t3=mcstat.mean_se(v_sum3,v_sq3,N)
assert abs(mean3)<3*err_t3,\
    "case3碰撞次数检验不合格"

#通量随深度先升后降，因为通量随深度的增加而减少，入射面表面没有反散射回流的中子，在表面附近反应的中子无法积累，
# 当到达一定深度后，反散射中子加未反应中子的衰减大于一开始的衰减。
#相对误差随深度单调增加，因为中子数量随深度越来越低，当到底最深处，相对误差最大，而这是最需要算准的区域。