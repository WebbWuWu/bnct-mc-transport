# -*- coding: utf-8 -*-
# ============================================================
#  lesson07_freegas.py  --  第 7 课：自由气体靶核热运动 + 上散射
#  2026-09-28 开文件（Claude：只有这段头注释、import、两个函数的
#  签名行和步骤注释；函数体、run_slab 的改动都由用户写）
# ============================================================
#
#  框架（9/28 用户定的九个决定）见 进度日志 9/28 条目。
#  lesson06_continuous.py 不动 —— 它是假库逐位对拍的基准。
#
#  讲义：项目文档\原理讲义_第7课_热化与自由气体模型.md
#        §4.2 抽靶核  §4.3 四步矢量变换  §8 判据
# ============================================================

import math
from rng import lcg, isotropic_direction
import os
import h5py
import numpy as np
from xs_lookup import LIB,sigma_at,find_cell

def path_A_prob(y):
    return y*math.sqrt(math.pi)/(y*math.sqrt(math.pi)+2)

def load_kT(nuc, T):
    path=os.path.join(LIB,nuc+".h5")
    assert os.path.exists(path),"文件不在这里:%s"%path
    with h5py.File(path,"r") as f:
        root=f[nuc]
        kT=float(root["kTs"][T][()])
        # [Claude assert] 抓：读错了温度那一项 / 单位不对。
        #   两条独立的路：库里存的 kT  vs  玻尔兹曼常数 × 标签上的温度。
        #   容差 0.5% 的来历：标签按整数 K 四舍五入，最多差 0.5 K / 250 K = 0.2%；
        #   294K 那一项实测差 0.13%（库按 293.6 K 做的）。读错成 250K 或 600K 那一项，差 15% 以上。
        T_label = float(T[:-1])
        assert abs(kT / (8.617333262e-5 * T_label) - 1.0) < 5e-3, \
            "kT=%r eV 与温度标签 %s 对不上：按标签应约为 %.6f eV，差超过 0.5%%" % (kT, T, 8.617333262e-5 * T_label)
        return kT

def sample_target(y, rng):
    """抽靶核（讲义 §4.2）。
    输入  y   = 中子速度 / 靶核最可几速度 = sqrt(E*A/kT)
          rng = 调用方传进来的发生器（第 4 问：③）
    输出  x, mu_T, 第几轮接受
    """
    # 1. 入口检查 y（Claude 给的 assert 1）
    assert y > 0 and math.isfinite(y), "y=%r 不是正的有限数（E 或 kT 出了问题）" % (y,)
    # 2. 最多 1000 轮：
    for k in range(1000):
        xi2=rng.random()
        xi3=rng.random()
        xi4=rng.random()
        if rng.random()<path_A_prob(y):
            x2=-math.log(1-xi2)-math.log(1-xi3)*math.cos((1-xi4)/2*math.pi)**2
             #A路配方
        else:
            x2=-math.log((1-xi2)*(1-xi3))
        x=math.sqrt(x2)
        mu_T=2*rng.random()-1
        arg = x2 + y**2 - 2*x*y*mu_T
        assert arg > -1e-12, "根号里为负且不是舍入误差：arg=%.3g x=%.6g y=%.6g mu_T=%.6g" % (arg, x, y, mu_T)
        v_rel = math.sqrt(max(arg, 0.0))
        assert v_rel <= x + y + 1e-12, "上界 x+y 失效：x=%.6g y=%.6g mu_T=%.6g v_rel=%.6g" % (x, y, mu_T, v_rel)
        pac=v_rel/(x+y)
        if rng.random()<pac:
            return x,mu_T,k+1
    else:
        assert False, "抽靶核 1000 轮全被拒（正常情况不可能，必是 bug）：y=%r" % (y,)


def free_gas_kinematics(E, u, v, w, A, kT, x, mu_T, rng):
    """四步矢量变换（讲义 §4.3）。
    输出  E', 新方向 u', v', w'
    """
    # 单位：全程以靶核最可几速度 v_T 为尺（9/29 ③ 丙）
    # 0.  入口：由 E、A、kT 算 y —— 和 sample_target 外面那次同一个公式
    y=math.sqrt(E*A/kT)
    # 1.  中子速度向量 = y × 中子方向 (u, v, w)
    vnx=y*u
    vny=y*v
    vnz=y*w
    # 2a. 最多 100 轮：各向同性抽方向 r → 减掉它沿中子方向的部分
    for k in range(100):
        uu,vv,ww=isotropic_direction(rng)
        al=uu*u+vv*v+ww*w
        rx=uu-u*al
        ry=vv-v*al
        rz=ww-w*al
        if math.sqrt(rx**2+ry**2+rz**2)>=1e-6:
            break
    else:
        # [Claude assert] for 的 else —— 抓：剩余长度永远不够（正常撞满概率约 1e-1230）
        assert False, "重抽 r 满 100 轮（正常不可能，必是 bug）：E=%r A=%r kT=%r y=%r" % (E, A, kT, y)
    #     → 剩下的长度 ≥ 1e-6 就跳出，否则重抽（② 丙 + 续2 甲）
    # 2b. 剩下的部分除以自己的长度 → e1
    e1x=rx/math.sqrt(rx**2+ry**2+rz**2)
    e1y=ry/math.sqrt(rx**2+ry**2+rz**2)
    e1z=rz/math.sqrt(rx**2+ry**2+rz**2)
    # 2c. 靶核方向 = mu_T × 中子方向 + sqrt(1 − mu_T²) × e1（不抽 phi、不造 e2：② 续 乙）
    rbx=mu_T*u+math.sqrt(1-mu_T**2)*e1x
    rby=mu_T*v+math.sqrt(1-mu_T**2)*e1y
    rbz=mu_T*w+math.sqrt(1-mu_T**2)*e1z
    # [Claude assert] 2c 之后 —— 抓：e1 没和中子方向垂直 / 没归一化、2c 某个分量抄错（判据 1 看不见这类错）
    nb2 = rbx**2 + rby**2 + rbz**2
    cb = rbx*u + rby*v + rbz*w
    assert abs(nb2 - 1.0) < 1e-9 and abs(cb - mu_T) < 1e-9, \
        "靶核方向不对：长度²=%.12g（应为 1），与中子方向夹角余弦=%.12g（应为 mu_T=%.12g）" % (nb2, cb, mu_T)
    # 2d. 靶核速度向量 = x × 靶核方向
    vbx=x*rbx
    vby=x*rby
    vbz=x*rbz
    # 3.  质心速度 = (中子速度 + A × 靶核速度) / (1 + A)
    vcmx=(vnx+A*vbx)/(1+A)
    vcmy=(vny+A*vby)/(1+A)
    vcmz=(vnz+A*vbz)/(1+A)
    # 4.  CM 里中子速度 = 中子速度 − 质心速度，记下速率
    vncmx=vnx-vcmx
    vncmy=vny-vcmy
    vncmz=vnz-vcmz
    vncm=math.sqrt(vncmx**2+vncmy**2+vncmz**2)
    # 5.  CM 里各向同性抽新方向，乘同一个速率
    u1,v1,w1=isotropic_direction(rng)
    v1x=u1*vncm
    v1y=v1*vncm
    v1z=w1*vncm
    # 6.  加回质心速度 → LAB 新速度 v' → E' = (kT/A) × |v'|² → 新方向 = v' / |v'|
    vlx=vcmx+v1x
    vly=vcmy+v1y
    vlz=vcmz+v1z
    vl=math.sqrt(vlx**2+vly**2+vlz**2)
    # [Claude assert] 除以 vl 之前 —— 抓：LAB 新速度为 0，方向无定义（④续 甲）
    assert vl > 0.0, "LAB 新速度为 0，方向无定义：E=%r A=%r kT=%r x=%r mu_T=%r" % (E, A, kT, x, mu_T)
    El=(kT/A)*(vl**2)
    u2=vlx/vl
    v2=vly/vl
    w2=vlz/vl
    #     （|v'| = 0 时由 Claude 给的 assert 当场停：④续 甲）
    return El,u2,v2,w2

def sigma_at_ext(E, XS, e):
    # [Claude assert] 抓：能量被上游算成负数 / nan / inf，却被当成「低能」往下外推。
    #   不加的话：负数 → math.sqrt 报 "math domain error"，报错的位置离真正的错很远；
    #            nan → 外推出 nan 截面，一路传到 run_slab 那条「Sigma_t 为零」，消息是错的。
    assert e > 0.0 and math.isfinite(e), "要查的能量 e=%r 不是正的有限数（上游算错了），不能外推" % (e,)
    if e>=E[0]:
        sig,total=sigma_at(E,XS,e)
    else:
        sig1,total1=sigma_at(E,XS,E[0])
        sig=[0.0]*len(sig1)
        for i in range(len(sig1)):
            sig[i]=sig1[i]*math.sqrt(E[0]/e)
        total=total1*math.sqrt(E[0]/e)
    return sig,total
# ------------------------------------------------------------
#  Claude 给的三条 assert（成品，去掉行首 # 粘到对应位置）
# ------------------------------------------------------------
# 入口 —— 抓：E 变负或 kT=0 导致 y 是 nan/inf/0
# assert y > 0 and math.isfinite(y), "y=%r 不是正的有限数（E 或 kT 出了问题）" % (y,)
#
# 每轮算完 v_rel 后 —— 抓：mu_T 抽的范围错了（phi 和 mu 那对易混）
# assert v_rel <= x + y + 1e-12, "上界 x+y 失效：x=%.6g y=%.6g mu_T=%.6g v_rel=%.6g" % (x, y, mu_T, v_rel)
#
# for 循环的 else 里 —— 抓：接受概率是 nan 等导致永远不接受（正常撞满概率约 1e-513）
# assert False, "抽靶核 1000 轮全被拒（正常情况不可能，必是 bug）：y=%r" % (y,)

def table_1overE(edges):
    m=[0.0]*(len(edges)-1)
    m1=[0.0]*(len(edges)-1)
    for i in range(len(edges)-1):
        m[i]=math.log(edges[i+1]/edges[i])
    tot=sum(m)
    for j in range(len(edges)-1):
        m1[j]=m[j]/tot
    return m1

def make_cdf(shares):
    tot=[0.0]
    for i in range(len(shares)):
        tot.append(sum(shares[0:(i+1)]))
    tot[-1]=1
    return tot

def sample_E_table(edges, cdf, rng):
    x1=rng.random()
    n1=find_cell(cdf,x1)
    le=edges[n1]
    re=edges[n1+1]
    x2=rng.random()
    E=le*(re/le)**x2
    return E
# ============================================================
#  run_slab（第 7 课版）
#  2026-09-29 Claude 从 lesson06_continuous.py 第 16–171 行原样复制（机械搬运，用户同意），
#  一个字未改；下面三行 import / 常量也是从 lesson06 开头原样抄来，函数要用。
#  改动由用户写。lesson06_continuous.py 本身不动。
# ============================================================
from xs_lookup import sigma_at, BARN
from sample_reaction import sample_reaction

Max=10000

def run_slab(N,seed,d,nuclides,E_src,order,z_src=None,t_cut=None,implicit=False, NL=20, src=None):
    if z_src is None:
        z_src=d/2

    assert len(order) == len(nuclides)

    for i in range(len(nuclides)):
        assert len(order[i]) == len(nuclides[i]["XS"])
    
    assert t_cut is None or t_cut > 0.0, "t_cut=%r 必须是 None（不截断）或正数" % (t_cut,)
    E_cut=[]
    len_nuc=len(nuclides)
    rng=lcg(seed)
    if src!=None:
        edges,shares=src
        cdf=make_cdf(shares)
    n_T=0
    n_R=0
    n_A=0
    n_stuck=0
    n_i=0
    nt_cut=0
    n_sd=0
    dz=d/NL
    track_sum=[0.0]*NL
    track_sq=[0.0]*NL
    v_sum=v_sq=0.0
    tol_sigt_all=0.0
    for i in range(N):
        wt=1.0
        t=0
        z=z_src
        if src==None:
            E=E_src
        else:
            E=sample_E_table(edges,cdf,rng)
        tol_sigt_s=0.0
        u=0
        v=0
        w=1.0
        track_h=[0.0]*NL
        n_hist=0
        for step in range(Max):
            sig_ds=[]
            sig_al=[]
            tot_ls=[]
            for j in range(len(nuclides)):
                Sig_j,tot_j=sigma_at_ext(nuclides[j]["E"],nuclides[j]["XS"],E)
                Sigma_j = nuclides[j]["dens"] * tot_j * BARN
                sig_ds.append(Sig_j)
                sig_al.append(Sigma_j)
                tot_ls.append(tot_j) 
            Sigma_t=sum(sig_al)
            assert Sigma_t > 0.0, "Sigma_t 为零，飞行距离会除零：E=%.6g eV" % E
            # 为什么是 1-rng.random() 而不是 rng.random()：
            #   random() 返回 [0,1) —— 0 取得到，1 取不到。若写 -log(xi)，xi=0 时 log 收到 0 → -inf，程序炸。
            #   取 1-xi 后 log 的参数落在 (0,1]，永远取不到 0。两种写法同分布，但安全性完全不同。
            #   命中概率约 2^-32：跑一亿条历史大概撞一次 —— 调试时永远不出现，演示时出现。
            #   这个选择配套的是本文件里 LCG 的区间，换 RNG 要重新核对。见 whiteboard_02_spec.md。
            s=-math.log(1-rng.random())/Sigma_t
            has_sign = False
            
            if t_cut!=None:
                t1=s/math.sqrt(E)
                if t+t1>t_cut:
                    s=(t_cut-t)*math.sqrt(E)
                    has_sign = True
                    t=t_cut
                else:
                    t+=t1
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
            tol_sigt_s += Sigma_t * seg
            z = z1                                     
            if z<0:
                n_R+=1
                break
            if z>d:
                n_T+=1
                break
            if has_sign ==True:
                nt_cut+=1
                E_cut.append(E)
                break
            share2=[0.0]*len(nuclides)
            si=[0.0]*len(nuclides)
            sig_s=[0.0]*len(nuclides)
            for j in range(len(nuclides)):
                si[j]=sig_al[j]/Sigma_t
                pl=nuclides[j]["mts"].index(2)
                sig_s[j]=sig_ds[j][pl]*nuclides[j]["dens"]*BARN
            sigma_s=sum(sig_s)
            for j in range(len(nuclides)):
                share2[j]=sig_s[j]/sigma_s
            if implicit==False:
                if len_nuc==1:
                    jn=0
                else:
                    xi1=rng.random()
                    tot=make_cdf(si)
                    jn=find_cell(tot,xi1)
                xi=rng.random()
                k=sample_reaction(sig_ds[jn],tot_ls[jn],xi,order[jn])
                is_scat=nuclides[jn]["mts"][k]==2
            else:
                wt=wt*sigma_s/Sigma_t
                if len_nuc==1:
                    jn=0
                else:
                    cdf2=make_cdf(share2)
                    xi3=rng.random()
                    jn=find_cell(cdf2,xi3)
                is_scat=True
            if is_scat:
                n_i+=1
                n_hist+=1
                if nuclides[jn]["awr"] is None:
                    u,v,w=isotropic_direction(rng)
                else:
                    A=nuclides[jn]["awr"]
                    kT=nuclides[jn]["kT"]
                    y=math.sqrt(E*A/kT)
                    x1,mu_T,k2=sample_target(y,rng)
                    E,u,v,w=free_gas_kinematics(E,u,v,w,A,kT,x1,mu_T,rng)
            else:
                n_A+=1
                n_i+=1
                n_hist+=1
                break
        else:
            if t_cut==None:
                n_sd+=1
            else:
                n_stuck+=1
        for j in range(NL):
            track_sum[j]+=track_h[j]
            track_sq[j]+=track_h[j]**2
        x = tol_sigt_s - n_hist
        tol_sigt_all += tol_sigt_s
        v_sum+=x
        v_sq+=x**2
    if n_stuck>0:
        print("warning!有%d条历史触顶max_event"% n_stuck)
    assert n_T+n_R+n_A+n_stuck+nt_cut+n_sd == N
    scale = abs(tol_sigt_all) + abs(n_i)
    assert math.isclose(tol_sigt_all - n_i, v_sum,rel_tol=0.0, abs_tol=1e-9*scale), \
                        "验证B总账不一致：累加器=%.6f  n_i=%d  v_sum=%.6f" % (tol_sigt_all, n_i, v_sum)
    # [Claude assert] 抓：E_cut 追加的位置 / 层放错（写到 break 后面 → 一条也记不上；放进每一步 → 记多了）
    assert len(E_cut) == nt_cut, "t_cut 时刻的能量记了 %d 条，但被截断的历史有 %d 条" % (len(E_cut), nt_cut)
    return n_T/N,n_R/N,n_A/N,n_i/N,track_sum,track_sq,dz,v_sum,v_sq,nt_cut,n_sd,n_stuck,E_cut


if __name__ == "__main__":
    E_f = [1e-5, 2e7]; XS_f = [[0.8, 0.8], [0.2, 0.2]]      # 第 6 课那个假库
    print(sigma_at_ext(E_f, XS_f, 3.1e-6))    # 低于表下限，要外推,sig=[1.44,0.36],total=1.8
    print(sigma_at_ext(E_f, XS_f, 0.37))      # 在表里，不外推,sig=[0.8,0.2],total=1.0
