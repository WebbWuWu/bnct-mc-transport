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
    '''返回 (nuclides, order)，和 run_slab 的两个参数一一对应'''
    stodic={}
    for nuc in nucs_w:
        E1,mts1,XS1,awr1,Q1=load_nuclide(nuc, T)
        E2,heat2=load_heating(nuc,T)
                # [Claude assert] 抓：两次读到的不是同一张能量网格（温度标签不同 / 读错核素）—— heat 会和截面错位，查出来全是隔壁格点的值
        assert len(E2) == len(E1) and np.array_equal(E2, E1), \
            "%s：load_heating 的能量网格（%d 点）和 load_nuclide 的（%d 点）不是同一张" % (nuc, len(E2), len(E1))
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
        # [Claude assert] 抓：漏乘 (1 − 硼份额)、硼份额加错地方 —— OpenMC 用 'wo' 会把总和悄悄归一化成 1，两边的数密度就对不上
    assert abs(sum(frac.values()) - 1.0) < 1e-14, \
        "三个质量份额加起来 = %.17g，不是 1（H %.9g、O %.9g、B10 %.3g）" % (sum(frac.values()), frac["H1"], frac["O16"], frac["B10"])
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
        # [Claude assert] 抓：份额配错核素、awr 用错、dens 公式里 ρ 或 N_A 漏乘
    #   两条路：① 水的化学式 —— 每个水分子 2 个 H、1 个 O，N_H / N_O 必须正好是 2
    #           ② 你 4-1、4-2 的手算（天然水 18.015、M(B10) 10.013；和库里质量差 ≤ 0.03%，所以容差取 1e-3）
    dens_of = {n["name"]: n["dens"] for n in nuclides}
    assert abs(dens_of["H1"] / dens_of["O16"] - 2.0) < 1e-12, \
        "N_H / N_O = %.15g，水里应正好是 2" % (dens_of["H1"] / dens_of["O16"])
    assert abs(dens_of["H1"] / 6.672e22 - 1) < 1e-3 and abs(dens_of["B10"] / 9.003e17 - 1) < 1e-3, \
        "数密度和 4-1、4-2 手算差超过 0.1%%：N_H = %.5e（手算 6.672e22），N_B10 = %.5e（手算 9.003e17）" % (dens_of["H1"], dens_of["B10"])
    return nuclides,order

if __name__ == "__main__":
    edges,shares=kur_table()
    E_keep=KEEP_FACTOR*edges[-1]
    nuclides,order=build_water(E_keep)
    for n in nuclides:
        print("name=%s,mts=%s,dens=%.4e"%(n["name"],
                                          n["mts"],
                                          n["dens"]))
        # [Claude 判据] build_water 的 heat 接对了没有（讲义 §8 第 ⑤ 条）
    #   路 A：从 build_water 拼好的字典里取 heat，用 sigma_at_ext 查（第二轮第 4 问 甲：当成只有一个道的 XS）
    #   路 B：重新从文件读一次，用 10/22 你 lesson11_mt301.py 里的查法（find_cell + xs_at）
    #   抓：heat 放错核素、被当成道筛短了、sigma_at_ext 传参顺序错、低于 E[0] 时外推没走到
    import math
    from lesson07_freegas import sigma_at_ext
    from xs_lookup import xs_at
    REF_1022 = {"H1": (438.03, 2.126e6), "O16": (0.0508, 8.983e5), "B10": (9.003e9, 1.143e6)}   # 10/22 你打印的数（3–4 位有效数字）
    for n in nuclides:
        assert len(n["heat"]) == len(n["E"]), "%s：heat 有 %d 个点，能量网格 %d 个 —— heat 被筛过了？" % (n["name"], len(n["heat"]), len(n["E"]))
        E_f, heat_f = load_heating(n["name"], T)
        for k, e in enumerate([0.0253, 1.0e6]):
            _, hA = sigma_at_ext(n["E"], [n["heat"]], e)
            hB = xs_at(E_f, heat_f, e, find_cell(E_f, e))
            assert hA == hB, "%s 在 %g eV：路 A %.17g ≠ 路 B %.17g" % (n["name"], e, hA, hB)
            assert abs(hA / REF_1022[n["name"]][k] - 1) < 2e-3, \
                "%s 在 %g eV：%.6g 和 10/22 的 %.6g 对不上（放错核素了？）" % (n["name"], e, hA, REF_1022[n["name"]][k])
        _, hA = sigma_at_ext(n["E"], [n["heat"]], 3.1e-6)
        assert hA == heat_f[0] * math.sqrt(E_f[0] / 3.1e-6), \
            "%s 在 3.1e-6 eV（低于表下限）：%.17g，按 1/v 外推应为 %.17g" % (n["name"], hA, heat_f[0] * math.sqrt(E_f[0] / 3.1e-6))
        print("  %s：heat 三点（0.0253 eV、1 MeV、3.1e-6 eV 外推）路 A = 路 B，和 10/22 对上" % n["name"])
    print("build_water 判据全过")
    pa=os.path.join(OUT_DIR,"water_N%d_s%d"%(N,SEED))
    if os.path.exists(pa + ".json"):
        r=load_result(pa)
    else:
        t1=time.time()
        r=run_slab(N, SEED, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL, src=(edges, shares))
        save_result(r,pa)
        t2=time.time()
        print("花费时间为%.6f"%(t2-t1))
        
