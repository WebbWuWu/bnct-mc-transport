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
from xs_lookup import LIB,sigma_at,find_cell,load_nuclide

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


# ---- 10/24 Claude 换成精简版（用户 13:36「你把这个 fgk 直接换掉」）：和 whiteboard\copy_fgk_v2.py 逐字相同，
#      物理和随机数用法一个没变，只是三个分量装进长度 3 的列表 + dot、norm 两个小函数。
#      判据：whiteboard\check_fgk.py 8 组逐位 + 随机数用量一致；check_split.py F1–F3、R1、R2 逐位 ----
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
wt_lo=0.1
wt_hi=0.2
Q_val=2341569.0#核数脚本 10/9，801 份额 0.937016

def prepare_nuclides(nuclides):
    ind_H=None
    ind_B=None
    loca_800=None
    loca_801=None
    Q_800=None
    Q_801=None
    ind_N=None
    loca_600=None
    Q_600=None
    ind_102=[]
    for j in range(len(nuclides)):
        if nuclides[j]["name"]=="B10":
            ind_B=j
            loca_800=nuclides[j]["mts"].index(800)
            loca_801=nuclides[j]["mts"].index(801)
            Q_800=nuclides[j]["Q_value"][loca_800]
            Q_801=nuclides[j]["Q_value"][loca_801]
        if nuclides[j]["name"]=="H1":
            ind_H=j
        if nuclides[j]["name"]=="N14":
            ind_N=j
            loca_600=nuclides[j]["mts"].index(600)
            Q_600=nuclides[j]["Q_value"][loca_600]
        loca_102=nuclides[j]["mts"].index(102)
        ind_102.append(loca_102)
    return {
        "ind_H":ind_H,
        "ind_B":ind_B,
        "loca_800":loca_800,
        "loca_801":loca_801,
        "Q_800":Q_800,
        "Q_801":Q_801,
        "ind_N":ind_N,
        "loca_600":loca_600,
        "Q_600":Q_600,
        "ind_102":ind_102
    }
def make_tallies(nuclides, NL):
    """整个模拟的计数：10 个计数器、12 个逐层数组、按核素的 sigma_x / sigma_x2、总账 3 个、E_cut。
    返回一个字典，键名 = 原来的变量名（10/23 拆分 ① 框架第 1 问 甲；Claude 从 run_slab 原样搬来，用户要求）。"""
    E_cut=[]
    n_T=0
    n_ki=0
    n_R=0
    n_A=0
    nt_w=0
    nr_w=0
    n_stuck=0
    n_i=0
    nt_cut=0
    n_sd=0
    track_sum=[0.0]*NL
    track_sq=[0.0]*NL
    Br_sum=[0.0]*NL
    Br_sq=[0.0]*NL
    Bc_sum=[0.0]*NL
    Bc_sq=[0.0]*NL
    Nc_sum=[0.0]*NL
    Nc_sq=[0.0]*NL
    ag_sum=[0.0]*NL
    ag_sq=[0.0]*NL
    DH_sum=[0.0]*NL
    DH_sq=[0.0]*NL
    sigma_x={}
    sigma_x2={}
    kerma_sum={}
    kerma_sq={}
    kerma_tot_sum={}
    kerma_tot_sq={}
    for i in range(len(nuclides)):
        name=nuclides[i]["name"]
        sigma_x[name]=[0.0]*NL
        sigma_x2[name]=[0.0]*NL
        kerma_sum[name]=[0.0]*NL
        kerma_sq[name]=[0.0]*NL
        kerma_tot_sum[name]=0.0
        kerma_tot_sq[name]=0.0
    gH_sum=[0.0]*NL
    gH_sq=[0.0]*NL
    gB_sum=[0.0]*NL
    gB_sq=[0.0]*NL
    v_sum=v_sq=0.0
    track_tot_sum=0.0
    track_tot_sq=0.0
    tol_sigt_all=0.0
    return {
        "E_cut":E_cut,
        "n_T":n_T,
        "n_ki":n_ki,
        "n_R":n_R,
        "n_A":n_A,
        "nt_w":nt_w,
        "nr_w":nr_w,
        "n_stuck":n_stuck,
        "n_i":n_i,
        "nt_cut":nt_cut,
        "n_sd":n_sd,
        "track_sum":track_sum,
        "track_sq":track_sq,
        "Br_sum":Br_sum,
        "Br_sq":Br_sq,
        "Bc_sum":Bc_sum,
        "Bc_sq":Bc_sq,
        "Nc_sum":Nc_sum,
        "Nc_sq":Nc_sq,
        "ag_sum":ag_sum,
        "ag_sq":ag_sq,
        "DH_sum":DH_sum,
        "DH_sq":DH_sq,
        "sigma_x":sigma_x,
        "sigma_x2":sigma_x2,
        "v_sum":v_sum,
        "v_sq":v_sq,
        "kerma_sum":kerma_sum,
        "kerma_sq":kerma_sq,
        "gH_sum":gH_sum,
        "gH_sq":gH_sq,
        "gB_sum":gB_sum,
        "gB_sq":gB_sq,
        "tol_sigt_all":tol_sigt_all,
        "track_tot_sum":track_tot_sum,
        "track_tot_sq":track_tot_sq,
        "kerma_tot_sum":kerma_tot_sum,
        "kerma_tot_sq":kerma_tot_sq,
    }

def finish(tal, N, dz):
    """收尾：撞满警告 + 三条总账 assert + 拼 28 个键的返回字典（run_slab 对外的说明书）。
    10/23 拆分 ① 框架第 3 问 甲；Claude 从 run_slab 原样搬来（用户要求）。"""
    # ---- 10/23 拆分 ① C（Claude 搬运，用户要求）：计数字典里的 28 样取回原名，下面的总账 assert、返回字典是 run_slab 原来的结尾，一字未改 ----
    E_cut=tal["E_cut"]
    n_T=tal["n_T"]
    n_ki=tal["n_ki"]
    n_R=tal["n_R"]
    n_A=tal["n_A"]
    nt_w=tal["nt_w"]
    nr_w=tal["nr_w"]
    n_stuck=tal["n_stuck"]
    n_i=tal["n_i"]
    nt_cut=tal["nt_cut"]
    n_sd=tal["n_sd"]
    track_sum=tal["track_sum"]
    track_sq=tal["track_sq"]
    Br_sum=tal["Br_sum"]
    Br_sq=tal["Br_sq"]
    Bc_sum=tal["Bc_sum"]
    Bc_sq=tal["Bc_sq"]
    Nc_sum=tal["Nc_sum"]
    Nc_sq=tal["Nc_sq"]
    ag_sum=tal["ag_sum"]
    ag_sq=tal["ag_sq"]
    DH_sum=tal["DH_sum"]
    DH_sq=tal["DH_sq"]
    sigma_x=tal["sigma_x"]
    sigma_x2=tal["sigma_x2"]
    v_sum=tal["v_sum"]
    v_sq=tal["v_sq"]
    tol_sigt_all=tal["tol_sigt_all"]
    if n_stuck>0:
        print("warning!有%d条历史触顶max_event"% n_stuck)
    scale = abs(tol_sigt_all) + abs(n_i)
    assert math.isclose(tol_sigt_all - n_i, v_sum,rel_tol=0.0, abs_tol=1e-9*scale), \
                        "验证B总账不一致：累加器=%.6f  n_i=%d  v_sum=%.6f" % (tol_sigt_all, n_i, v_sum)
    # [Claude assert] 抓：某种结局的历史没进任何计数器（比如轮盘赌杀掉的），或者一条历史被数了两次。
    #   两条路：各个结局计数器之和  vs  历史循环跑了 N 圈。
    assert n_T + n_R + n_A + n_stuck + nt_cut + n_sd + n_ki == N, \
        "历史去向对不上：漏出 %d+%d、吸收 %d、撞满 %d+%d、截断 %d、轮盘赌 %d，合计 %d ≠ N=%d" \
        % (n_T, n_R, n_A, n_stuck, n_sd, nt_cut, n_ki,
           n_T + n_R + n_A + n_stuck + n_sd + nt_cut + n_ki, N)
    # [Claude assert] 抓：E_cut 追加的位置 / 层放错（写到 break 后面 → 一条也记不上；放进每一步 → 记多了）
    assert len(E_cut) == nt_cut, "t_cut 时刻的能量记了 %d 条，但被截断的历史有 %d 条" % (len(E_cut), nt_cut)
    # 10/9 第 12 条：返回值改字典（框架第 9 问 丙）。这一块由 Claude 代写（用户 20:26 要求），共 26 项。
    #   前 4 项是「除以 N 之后的比例 / 平均数」，键名带 frac_ / per_hist，免得和条数混。
    return {
        # ---- 原来的 13 个（值一个不改，顺序同原来的元组）----
        "frac_T": n_T/N,            # 透射比例（条数 / N）
        "frac_R": n_R/N,            # 反射比例
        "frac_A": n_A/N,            # 真吸收比例（隐式模式下恒为 0）
        "coll_per_hist": n_i/N,     # 平均每条历史碰撞几次
        "track_sum": track_sum,
        "track_sq": track_sq,
        "dz": dz,
        "v_sum": v_sum,
        "v_sq": v_sq,
        "nt_cut": nt_cut,
        "n_sd": n_sd,
        "n_stuck": n_stuck,
        "E_cut": E_cut,
        # ---- 10/9 新加的 13 个 ----
        "Br_sum": Br_sum,           # 硼逐道（eV）
        "Br_sq": Br_sq,
        "Bc_sum": Bc_sum,           # 硼常数对照（eV）
        "Bc_sq": Bc_sq,
        "Nc_sum": Nc_sum,           # 氮（eV）
        "Nc_sq": Nc_sq,
        "ag_sum": ag_sum,           # (n,γ) 反应次数
        "ag_sq": ag_sq,
        "DH_sum": DH_sum,           # D_H 反冲（eV）
        "DH_sq": DH_sq,
        "nt_w": nt_w,               # 透射权重和
        "nr_w": nr_w,               # 反射权重和
        "n_ki": n_ki,               # 轮盘赌杀掉的条数
        "sigma_x":sigma_x,
        "sigma_x2":sigma_x2,
        "kerma_sum":tal["kerma_sum"],
        "kerma_sq":tal["kerma_sq"],
        "gH_sum":tal["gH_sum"],
        "gH_sq":tal["gH_sq"],
        "gB_sum":tal["gB_sum"],
        "gB_sq":tal["gB_sq"],
        "track_tot_sum":tal["track_tot_sum"],
        "track_tot_sq":tal["track_tot_sq"],
        "kerma_tot_sum":tal["kerma_tot_sum"],
        "kerma_tot_sq":tal["kerma_tot_sq"],
    }

# ============================================================
#  10/23 拆分 ②（Claude 搬运，用户要求「直接拆完」）：历史循环拆成 5 个小函数。
#  每一行都是从原 run_slab 历史循环里原样搬来的，只是把名字换成从两个小字典里取：
#    p     粒子状态（这条中子现在在哪、往哪飞、多少能量、多少权重、走了多久）：z u v w E wt t
#    hist  这条历史自己的计数（每条历史开头清零，纪律 6）：track_h B_r B_c N_ar a_g D_H d_his n_hist tol_sigt_s
#  规矩（和拆分 ① 一样）：idx 只读；tal、hist、p 在函数里「改」，不「换」。
#  随机数调用的先后顺序一个都没动 —— check_split.py 逐位回归守着它。
# ============================================================

def start_history(nuclides, NL, z_src, E_src, src, edges, cdf, rng):
    """一条新历史：粒子从源出发（p），这条历史的计数全部清零（hist）。原历史循环第 1–22 行。"""
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
    B_r=[0.0]*NL
    B_c=[0.0]*NL
    N_ar=[0.0]*NL
    a_g=[0.0]*NL
    D_H=[0.0]*NL
    gH_h=[0.0]*NL
    gB_h=[0.0]*NL
    d_his={}
    kerma_h={}
    for m in range(len(nuclides)):           # 原来这里用的是 i，和外层历史循环的 i 同名（10/17 挂着的那条）；挪进函数后改成 m，结果不变
        name=nuclides[m]["name"]
        d_his[name]=[0.0]*NL
        kerma_h[name]=[0.0]*NL
    n_hist=0
    p={"z":z, "u":u, "v":v, "w":w, "E":E, "wt":wt, "t":t}
    hist={"track_h":track_h, "B_r":B_r, "B_c":B_c, "N_ar":N_ar, "a_g":a_g, "D_H":D_H,
          "d_his":d_his, "n_hist":n_hist, "tol_sigt_s":tol_sigt_s,"kerma_h":kerma_h
          ,"gH_h":gH_h,"gB_h":gB_h}
    return p, hist


def lookup_xs(nuclides, idx, E):
    """在能量 E 处查好这一步要用的所有截面和「每 cm 记多少」的系数。飞行中 E 不变，所以每步查一次。
    M2 的 kerma 系数（Σ N_i·H_i(E)）以后就加在这里。原历史循环每步第 1–28 行。"""
    sig_ds=[]
    sig_al=[]
    tot_ls=[]
    kerco={}
    for j in range(len(nuclides)):
        Sig_j,tot_j=sigma_at_ext(nuclides[j]["E"],nuclides[j]["XS"],E)
        Sigma_j = nuclides[j]["dens"] * tot_j * BARN
        sig_ds.append(Sig_j)
        sig_al.append(Sigma_j)
        tot_ls.append(tot_j)
        _,heat=sigma_at_ext(nuclides[j]["E"],[nuclides[j]["heat"]],E)
        ker_coe=nuclides[j]["dens"]*BARN*heat
        kerco[nuclides[j]["name"]]=ker_coe
    Sigma_t=sum(sig_al)
    assert Sigma_t > 0.0, "Sigma_t 为零，飞行距离会除零：E=%.6g eV" % E
    B_reac=0.0
    B_con=0.0
    N_coe=0.0
    a_coe=0.0
    gB_coe=0.0
    gH_coe=0.0
    if idx["ind_B"] !=None:
        sig800=sig_ds[idx["ind_B"]][idx["loca_800"]]
        sig801=sig_ds[idx["ind_B"]][idx["loca_801"]]
        tran_coe=nuclides[idx["ind_B"]]["dens"]*BARN
        B_reac=tran_coe*(sig800*idx["Q_800"]+sig801*idx["Q_801"])
        B_con=tran_coe*(sig800+sig801)*Q_val
        gB_coe=tran_coe*sig801
    if idx["ind_H"] !=None:
        sigH1=sig_ds[idx["ind_H"]][idx["ind_102"][idx["ind_H"]]]
        gH_coe=nuclides[idx["ind_H"]]["dens"]*BARN*sigH1
    if idx["ind_N"] !=None:
        sig600=sig_ds[idx["ind_N"]][idx["loca_600"]]
        N_coe=nuclides[idx["ind_N"]]["dens"]*BARN*sig600*idx["Q_600"]
    for j in range(len(nuclides)):
        sig102=sig_ds[j][idx["ind_102"][j]]
        a_coe+=nuclides[j]["dens"]*BARN*sig102
    return {"sig_ds":sig_ds, "sig_al":sig_al, "tot_ls":tot_ls, "Sigma_t":Sigma_t,
            "B_reac":B_reac, "B_con":B_con, "N_coe":N_coe, "a_coe":a_coe,
            "kerco":kerco,"gH_coe":gH_coe,"gB_coe":gB_coe}


def fly(p, xs, hist, rng, d, dz, NL, t_cut):
    """飞一段：抽飞行距离 → 时间截断 → 这一段按层拆开记账（径迹长度估计）→ 更新位置。
    改 p（z、t）和 hist；返回这一段是不是被时间截断了（原来的 has_sign）。
    M5 三维只改这一个函数。原历史循环每步第 29–111 行。"""
    E=p["E"]
    w=p["w"]
    wt=p["wt"]
    z=p["z"]
    t=p["t"]
    Sigma_t=xs["Sigma_t"]
    B_reac=xs["B_reac"]
    B_con=xs["B_con"]
    N_coe=xs["N_coe"]
    a_coe=xs["a_coe"]
    track_h=hist["track_h"]
    B_r=hist["B_r"]
    B_c=hist["B_c"]
    N_ar=hist["N_ar"]
    a_g=hist["a_g"]
    kerco=xs["kerco"]
    gH_coe=xs["gH_coe"]
    gB_coe=xs["gB_coe"]
    kerma_h=hist["kerma_h"]
    gH_h=hist["gH_h"]
    gB_h=hist["gB_h"]
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
        track_h[k0] += s*wt
        B_r[k0]+=wt*s*B_reac
        B_c[k0]+=wt*s*B_con
        N_ar[k0]+=wt*s*N_coe
        a_g[k0]+=wt*s*a_coe
        for name in sorted(kerco.keys()):
            kerma_h[name][k0]+=wt*s*kerco[name]
        gH_h[k0]+=wt*s*gH_coe
        gB_h[k0]+=wt*s*gB_coe
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
            track_h[k] += wt*dz_step / abs(w)
            B_r[k]+=wt*dz_step / abs(w)*B_reac
            B_c[k]+=wt*dz_step / abs(w)*B_con
            N_ar[k]+=wt*dz_step / abs(w)*N_coe
            a_g[k]+=wt*dz_step / abs(w)*a_coe
            for name in sorted(kerco.keys()):
                kerma_h[name][k]+=wt*dz_step / abs(w)*kerco[name]
            gH_h[k]+=wt*dz_step / abs(w)*gH_coe
            gB_h[k]+=wt*dz_step / abs(w)*gB_coe
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
    hist["tol_sigt_s"] += Sigma_t * seg
    p["z"] = z1
    p["t"] = t
    return has_sign


def collide(p, xs, hist, tal, nuclides, order, rng, implicit, dz):
    """撞一下：抽核素 → 抽道（模拟）或减权重 + 轮盘赌（隐式）→ 散射（改 E、方向，记 D_H）或吸收。
    改 p（E、u、v、w、wt）、hist、tal；返回这条中子还活着没有（False = 吸收或被轮盘赌杀掉，历史结束）。
    M3 的 S(α,β)、真实角分布都改在这里。原历史循环每步第 129–184 行。"""
    E=p["E"]
    u=p["u"]
    v=p["v"]
    w=p["w"]
    wt=p["wt"]
    z=p["z"]
    sig_ds=xs["sig_ds"]
    sig_al=xs["sig_al"]
    tot_ls=xs["tot_ls"]
    Sigma_t=xs["Sigma_t"]
    len_nuc=len(nuclides)
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
        if wt<wt_lo:
            xi4=rng.random()
            if xi4<wt/wt_hi:
                wt=wt_hi
            else:
                tal["n_ki"]+=1
                p["wt"]=wt
                return False
        if len_nuc==1:
            jn=0
        else:
            cdf2=make_cdf(share2)
            xi3=rng.random()
            jn=find_cell(cdf2,xi3)
        is_scat=True
    p["wt"]=wt
    if is_scat:
        tal["n_i"]+=1
        hist["n_hist"]+=1
        E_be=E
        if nuclides[jn]["awr"] is None:
            u,v,w=isotropic_direction(rng)
        else:
            A=nuclides[jn]["awr"]
            kT=nuclides[jn]["kT"]
            y=math.sqrt(E*A/kT)
            x1,mu_T,k2=sample_target(y,rng)
            E,u,v,w=free_gas_kinematics(E,u,v,w,A,kT,x1,mu_T,rng)
        k_hit=int(z/dz)
        hist["D_H"][k_hit]+=wt*(E_be-E)
        nam=nuclides[jn]["name"]
        hist["d_his"][nam][k_hit]+=wt*(E_be-E)
        p["E"]=E
        p["u"]=u
        p["v"]=v
        p["w"]=w
        return True
    else:
        tal["n_A"]+=1
        tal["n_i"]+=1
        hist["n_hist"]+=1
        return False


def end_history(tal, hist, nuclides, NL):
    """一条历史结束：把它自己的计数加进整个模拟的 sum 和 sq（第 4 课：每条历史是一个样本）。原历史循环最后 21 行。"""
    track_h=hist["track_h"]
    B_r=hist["B_r"]
    B_c=hist["B_c"]
    N_ar=hist["N_ar"]
    a_g=hist["a_g"]
    D_H=hist["D_H"]
    d_his=hist["d_his"]
    kerma_h=hist["kerma_h"]
    gH_h=hist["gH_h"]
    gB_h=hist["gB_h"]
    for j in range(NL):
        for m in range(len(nuclides)):
            name=nuclides[m]["name"]
            tal["sigma_x"][name][j]+=d_his[name][j]
            tal["sigma_x2"][name][j]+=d_his[name][j]**2
            tal["kerma_sum"][name][j]+=kerma_h[name][j]
            tal["kerma_sq"][name][j]+=kerma_h[name][j]**2
        tal["gH_sum"][j]+=gH_h[j]
        tal["gH_sq"][j]+=gH_h[j]**2
        tal["gB_sum"][j]+=gB_h[j]
        tal["gB_sq"][j]+=gB_h[j]**2
        tal["track_sum"][j]+=track_h[j]
        tal["track_sq"][j]+=track_h[j]**2
        tal["Br_sum"][j]+=B_r[j]
        tal["Br_sq"][j]+=B_r[j]**2
        tal["Bc_sum"][j]+=B_c[j]
        tal["Bc_sq"][j]+=B_c[j]**2
        tal["Nc_sum"][j]+=N_ar[j]
        tal["Nc_sq"][j]+=N_ar[j]**2
        tal["ag_sum"][j]+=a_g[j]
        tal["ag_sq"][j]+=a_g[j]**2
        tal["DH_sum"][j]+=D_H[j]
        tal["DH_sq"][j]+=D_H[j]**2
    sum_th=sum(track_h)
    tal["track_tot_sum"]+=sum_th
    tal["track_tot_sq"]+=sum_th**2
    for n in range(len(nuclides)):
        name=nuclides[n]["name"]
        sum_kh=sum(kerma_h[name])
        tal["kerma_tot_sum"][name]+=sum_kh
        tal["kerma_tot_sq"][name]+=sum_kh**2
    x = hist["tol_sigt_s"] - hist["n_hist"]
    tal["tol_sigt_all"] += hist["tol_sigt_s"]
    tal["v_sum"]+=x
    tal["v_sq"]+=x**2


def run_slab(N,seed,d,nuclides,E_src,order,z_src=None,t_cut=None,implicit=False, NL=20, src=None):
    if z_src is None:
        z_src=d/2

    assert len(order) == len(nuclides)

    for i in range(len(nuclides)):
        assert len(order[i]) == len(nuclides[i]["XS"])

    assert t_cut is None or t_cut > 0.0, "t_cut=%r 必须是 None（不截断）或正数" % (t_cut,)
    idx=prepare_nuclides(nuclides)
    rng=lcg(seed)
    edges=None
    cdf=None
    if src!=None:
        edges,shares=src
        cdf=make_cdf(shares)
    tal=make_tallies(nuclides, NL)      # 10/23 拆分 ① B：全模拟的计数装进一个字典（框架第 1 问 甲）
    dz=d/NL
    for i in range(N):
        p, hist = start_history(nuclides, NL, z_src, E_src, src, edges, cdf, rng)
        for step in range(Max):
            xs = lookup_xs(nuclides, idx, p["E"])
            has_sign = fly(p, xs, hist, rng, d, dz, NL, t_cut)
            # 这一段飞完：漏出 / 截断 → 历史结束（三种结局的计数留在这里，一眼能看全）
            if p["z"]<0:
                tal["n_R"]+=1
                tal["nr_w"]+=p["wt"]
                break
            if p["z"]>d:
                tal["n_T"]+=1
                tal["nt_w"]+=p["wt"]
                break
            if has_sign ==True:
                tal["nt_cut"]+=1
                tal["E_cut"].append(p["E"])
                break
            if not collide(p, xs, hist, tal, nuclides, order, rng, implicit, dz):
                break
        else:
            if t_cut==None:
                tal["n_sd"]+=1
            else:
                tal["n_stuck"]+=1
        end_history(tal, hist, nuclides, NL)
    return finish(tal, N, dz)          # 10/23 拆分 ① C



if __name__ == "__main__":
    E_f = [1e-5, 2e7]; XS_f = [[0.8, 0.8], [0.2, 0.2]]      # 第 6 课那个假库
    print(sigma_at_ext(E_f, XS_f, 3.1e-6))    # 低于表下限，要外推,sig=[1.44,0.36],total=1.8
    print(sigma_at_ext(E_f, XS_f, 0.37))      # 在表里，不外推,sig=[0.8,0.2],total=1.0
