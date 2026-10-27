# -*- coding: utf-8 -*-
# ============================================================
#  whiteboard/copy_build_water.py  --  build_water 抄写稿（10/27 抄，10/28 默写（1））
#  Claude 2026-10-26 晚开好文件头；从「你从这里往下写」那一行下面开始，照 lesson11_water.py 手敲
# ============================================================
#
#  照着抄的原文：lesson11_water.py
#    · 第 1–11 行：import + 三个常数（rho_w、share_B、nucs_w）
#    · 第 13–81 行：def build_water(E_keep) 整个函数
#
#  要写：
#    · 用到的 import（只抄 build_water 真用到的：load_nuclide、find_cell、load_kT、load_heating、numpy、T、M_N、N_A）
#    · 三个常数
#    · build_water 函数体，分四块：
#        ① 读三个核素：每个核素读一次截面（load_nuclide）、一次 MT=301（load_heating），装进 stodic[核素名]
#        ② 质量份额：H、O 按 awr × M_N 算水分子里的质量占比，再乘 (1 − 硼份额)；B10 = 15e-6 → 三个和为 1
#        ③ 每个核素：筛道（E < E_keep 范围里有非零截面的才留）→ 数密度 dens = 份额 × ρ × N_A /(awr × M_N)
#           → 读 kT → 拼核素字典（9 个键）→ 在 0.0253 eV 处按截面从大到小排 order
#        ④ return nuclides, order
#
#  不写：
#    · 所有 assert 和它们上面那行 [Claude assert] 注释（看一眼就行，10/14 定）
#    · print
#    · if __name__ == "__main__": 整块（10/10 定）
#
#  抄的时候顺手看清楚三件事（10/28 默写要用，也是老师问答会问的）：
#    1. 为什么 H 的份额要乘 2、O 不乘？「(1 − share_B)」乘在哪几个上？
#    2. dens 公式里每一项的单位：份额（无量纲）× ρ（g/cm³）× N_A（/mol）÷ 摩尔质量（g/mol）→ ?
#    3. 筛道条件 np.any(XS[i][E < E_keep] > 0) 是什么意思：E < E_keep 先得到什么，再拿去干什么？
#
#  下面两行是 Claude 加的（不用抄）：这个文件在 whiteboard\ 里，要能 import my_code 根目录的模块
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# ============================ 你从这里往下写 ============================
from xs_lookup import load_nuclide,find_cell
from lesson07_freegas import load_kT,run_slab
from lesson11_mt301 import load_heating
import numpy as np
from lesson08_dose import T,M_N,N_A,kur_table,KEEP_FACTOR
import os,time
from lesson08_dose import save_result,load_result,D,NL,N,SEED,OUT_DIR

rho_w=0.998
share_B=15e-6
nucs_w=["H1","O16","B10"]

def build_water(E_keep):
    stodic={}
    for nuc in nucs_w:
        E1,mts1,XS1,awr1,Q1=load_nuclide(nuc,T)
        E2,heat2=load_heating(nuc,T)
        stodic[nuc]={
            "E":E1,
            "mts":mts1,
            "XS":XS1,
            "awr":awr1,
            "Q":Q1,
            "heat":heat2
        }
    M_H=M_N*stodic["H1"]["awr"]
    M_O=M_N*stodic["O16"]["awr"]
    M_water=2*M_H+M_O
    share_H=M_H*2/M_water*(1-share_B)
    share_O=M_O/M_water*(1-share_B)
    frac={
        "H1":share_H,
        "O16":share_O,
        "B10":share_B
    }
    nuclides=[]
    order=[]
    for nuc in nucs_w:
        keys=["E","mts","XS","awr","Q","heat"]
        E,mts,XS,awr,Q,heat=[stodic[nuc][k] for k in keys]
        mts1=[]
        XS1=[]
        Q1=[]
        for i in range(len(XS)):
            if np.any(XS[i][E<E_keep]>0):
                mts1.append(mts[i])
                XS1.append(XS[i])
                Q1.append(Q[i])
        dens=frac[nuc]*rho_w*N_A/(awr*M_N)
        kT=load_kT(nuc,T)
        dirc={
            "E":E,
            "mts":mts1,
            "XS":XS1,
            "awr":awr,
            "Q_value":Q1,
            "name":nuc,
            "dens":dens,
            "kT":kT,
            "heat":heat
        }
        loa_e=find_cell(E,0.0253)
        order.append(sorted(range(len(mts1)),key=lambda t:XS1[t][loa_e],reverse=True))
        nuclides.append(dirc)
    return nuclides,order
        
