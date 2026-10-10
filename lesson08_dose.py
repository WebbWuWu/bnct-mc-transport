import math
import os
import time
import json
import numpy as np
import matplotlib.pyplot as plt
from xs_lookup import load_nuclide,find_cell
from lesson07_freegas import run_slab,load_kT,table_1overE
# ============================================================
#  常数（Claude 10/10 代写，用户要求；出处写在每行后面）
# ============================================================
T = "294K"                                   # 截面温度，和 load_nuclide / load_kT 用的同一个标签

# ---- 组织（ICRU-44 软组织，第 5 课；只取 H、C、N、O 四种主元素 + 硼）----
RHO = 1.06                                   # g/cm³
NUCS = ["H1", "C12", "N14", "O16", "B10"]    # 核素名 = 库文件名 = 字典里的 "name"
MASS_FRAC = {"H1": 0.102, "C12": 0.143, "N14": 0.034, "O16": 0.708,
             "B10": 15e-6}                    # 质量份额；B10 = 15 ppm（正常组织）
M_N = 1.00866491588                          # 中子摩尔质量 g/mol（awr × M_N = 核素摩尔质量）
N_A = 6.02214076e23                          # 阿伏伽德罗常数 /mol（定义值）

# ---- 换算（讲义 §2.1）----
EV2J = 1.602176634e-19                       # J/eV（定义值）—— 🔴 不是 ρ 的 1.06
G_PER_KG = 1.0e3                             # 1 kg = 1000 g：eV/g × 1000 = eV/kg

# ---- 源：京大 KUR 临床束流（Fujiwara 等 2013, J Radiat Res 54:769, 表 1）----
KUR_FLUX = [3.0e7, 7.3e8, 4.7e7]             # 热 / 超热 / 快 三群通量 cm⁻²s⁻¹（原文数，份额在代码里算）
E_TH_LO = 0.01                               # eV，热群下限（10/10 定，近似）
E_FAST_HI = 2.0e6                            # eV，快群上限（10/10 定，近似）
EPI_EDGES = [0.5, 1.0, 3.0, 10.0, 50.0, 200.0, 1000.0, 3000.0, 1.0e4]   # 超热段箱边界（10/4 定）
KEEP_FACTOR = 2.0                            # E_keep = KEEP_FACTOR × 源最高能量（10/10 丙，余量）

# ---- 几何与运行 ----
D = 15.0                                     # 板厚 cm
NL = 30                                      # 层数 → dz = 0.5 cm
N = 100000                                    # 历史数：先 2 万出草图，交付时改 100000
SEED = 1010                                  # 种子
TUMOR_RATIO = 3.5                            # 肿瘤 52.5 ppm / 正常 15 ppm（10/4 定）

# ---- 输出 ----
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "outputs")

def kur_table():
    """KUR 三群 → 10 箱源表。返回 (edges, shares)：edges 11 个，shares 10 个。"""
    # 1. 三群通量 → 三个份额（各除以三者之和）
    sum_fluk=sum(KUR_FLUX)
    FLUK_share=[]
    for i in range(len(KUR_FLUX)):
        FLUK_share.append(KUR_FLUX[i]/sum_fluk)
    # 2. 超热段：table_1overE(EPI_EDGES) 给段内 8 个箱的占比 → 每个都乘超热份额
    Su_Th=[]
    m1=table_1overE(EPI_EDGES)
    for j in range(len(m1)):
        Su_Th.append(m1[j]*FLUK_share[1])
    # 3. edges = 热下限 + EPI_EDGES + 快上限
    edges=[E_TH_LO]
    for i in range(len(EPI_EDGES)):
        edges.append(EPI_EDGES[i])
    edges.append(E_FAST_HI)
    # 4. shares = [热份额] + 超热 8 个 + [快份额]
    shares=[FLUK_share[0]]
    for i in range(len(Su_Th)):
        shares.append(Su_Th[i])
    shares.append(FLUK_share[2])
    # 5. 返回
        # [Claude assert] 抓：箱边界和份额错位（多 / 少一个）—— sample_E_table 会取到错的箱或越界
    assert len(edges) == len(shares) + 1, "edges 有 %d 个、shares 有 %d 个，应该正好多一个" % (len(edges), len(shares))
    # [Claude assert] 抓：某一段份额没乘对（比如超热段乘错群）—— 份额之和不是 1，make_cdf 末端会被硬改成 1，偏差被藏起来
    assert abs(sum(shares) - 1.0) < 1e-12, "份额之和 = %.15f，不是 1" % sum(shares)
    # [Claude assert] 抓：拼接顺序错（热下限、快上限放错位置）—— 箱边界必须严格递增
    for k in range(len(edges) - 1):
        assert edges[k] < edges[k + 1], "edges 第 %d、%d 个不递增：%r ≥ %r" % (k, k + 1, edges[k], edges[k + 1])
    return edges,shares


def build_tissue(E_keep):
    """真库软组织。返回 (nuclides, order)，和 run_slab 的两个参数一一对应。"""
    nuclides=[]
    order=[]
    # 对 NUCS 里每个名字：
    for name in NUCS:
        lef=[]
        # 1. load_nuclide(名字, T) → 5 个返回值
        E,mts,XS,awr,Q=load_nuclide(name,T)
        # 2. 筛道：能量 < E_keep 的格点上截面不全为 0 的道才留；mts、XS、Q 用同一组下标一起筛
        for i in range(len(XS)):
            if np.any(XS[i][E<E_keep]>0):
                lef.append(i)
        mts2=[mts[i] for i in lef]
        XS2=[XS[i] for i in lef]
        Q2=[Q[i] for i in lef]
        # 3. dens = 质量份额 × ρ × N_A / (awr × M_N)
        dens=MASS_FRAC[name]*RHO*N_A/(awr*M_N)
        # 4. kT = load_kT(名字, T)
        kT=load_kT(name,T)
        # 5. 拼 8 个键的字典 → 放进 nuclides
        dirc={
            "E":E,
            "mts":mts2,
            "XS":XS2,
            "awr":awr,
            "Q_value":Q2,
            "name":name,
            "dens":dens,
            "kT":kT
        }
        loa_e=find_cell(E,0.0253)
        order.append(sorted(range(len(mts2)),key=lambda t:XS2[t][loa_e],reverse=True))
        # 6. 这个核素的 order：留下的道按 0.0253 eV 处截面从大到小排的【下标】 → 放进 order
        nuclides.append(dirc)
    # 7. 返回
    return nuclides,order


def save_result(r, path_stem):
    """把 run_slab 返回的字典存两份：path_stem + ".json" 和 path_stem + ".npz"。"""
    # 1. 存 json
    with open(path_stem+".json","w") as f:
        json.dump(r,f)
    # 2. 存 npz
    np.savez(path_stem+".npz",**r)


def load_result(path_stem):
    """只读 json，返回字典。"""
    # 1. 读 json → 返回
    with open(path_stem+".json","r") as f:
        return json.load(f)


def to_gy(s, sq, N, dz):
    """每层 eV/源中子 → Gy·cm²。返回 (均值列表, SE 列表)。"""
    # 1. 系数 = EV2J × G_PER_KG / (RHO × dz)
    coe=EV2J*G_PER_KG / (RHO*dz)
    # 2. 每层均值、SE（第 4 课公式，和 layer_mean_se 同一个）
    ave=[]
    SE=[]
    for k in range(len(s)):
        m=s[k]/N
        var=(sq[k]/N-m**2)*N/(N-1)
        se=math.sqrt(var/N)
        ave.append(m*coe)
        SE.append(se*coe)
    return ave,SE


if __name__ == "__main__":
    # 1. edges, shares = kur_table()；E_keep = KEEP_FACTOR × 源最高能量；（Claude 给 assert）
    edges,shares=kur_table()
    E_keep=KEEP_FACTOR*edges[-1]
    # [Claude assert] 抓：E_keep 没跟着源变（写死成旧的 2e4 之类）—— 源里最快的中子会碰到被筛掉的道，静默丢反应
    assert E_keep > edges[-1], "E_keep = %r eV 不大于源的最高能量 %r eV，会静默丢掉阈值反应" % (E_keep, edges[-1])
    # 2. nuclides, order = build_tissue(E_keep)
    nuclides, order = build_tissue(E_keep)
    for nuc in nuclides:
        print("name=%s,mts=%s"%(nuc["name"],nuc["mts"]))
    # 3. 文件名 stem = OUT_DIR 下的 "dose_N{N}_s{SEED}"；json 已存在 → load_result；不存在 → run_slab → save_result
    #    run_slab 参数：隐式、D、NL、z_src=0.0、src=(edges, shares)
    stem=os.path.join(OUT_DIR,"dose_N%d_s%d"%(N,SEED))
    if os.path.exists(stem + ".json"):
        r=load_result(stem)
    else:
        t1=time.time()
        r = run_slab(N, SEED, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL, src=(edges, shares))
        save_result(r, stem)
        t2=time.time()
        print("跑一次花费时间为t=%.6f"%(t2-t1))
    
    # 4. to_gy 换算：硼逐道 Br、氮 Nc、D_H DH；肿瘤 D_B = 硼的均值和 SE 都 × TUMOR_RATIO
    dz=r["dz"]
    D_B,DB_SE=to_gy(r["Br_sum"],r["Br_sq"],N,dz)
    D_N,DN_SE=to_gy(r["Nc_sum"],r["Nc_sq"],N,dz)
    D_H,DH_SE=to_gy(r["DH_sum"],r["DH_sq"],N,dz)
    tumer1=[]
    tumer2=[]
    for i in range(len(D_B)):
        tumer1.append(D_B[i]*TUMOR_RATIO)
        tumer2.append(DB_SE[i]*TUMOR_RATIO)
    D_Bc,_=to_gy(r["Bc_sum"], r["Bc_sq"], N, dz)
    com=[]
    for k in range(NL):
        com.append(abs(D_Bc[k]/D_B[k]-1)) 
    print("30层里最大的是%.2e"%max(com))
    
    #    Bc 不画，打印它和 Br 的最大相对差
    # 5. (n,γ)：ag 每层均值和 SE ÷ dz（单位 cm⁻¹，即每 cm³ 反应数 / 单位源注量）
    ng=[]
    ng_se=[]
    for k in range(NL):
        m = r["ag_sum"][k] / N
        var = (r["ag_sq"][k]/N - m**2) * N/(N-1)
        se = math.sqrt(var/N)
        ng.append(m/dz)
        ng_se.append(se/dz)
    # 6. 横轴 = 层中心 (k + 0.5) × dz
    depth=[]
    for k in range(NL):
        depth.append((k + 0.5) * dz)
    # 7. 上下两个子图，共用横轴，纵轴对数：上 = D_B、D_N、D_H 带误差棒 + 肿瘤 D_B 虚线；下 = (n,γ)
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.errorbar(depth, D_B, yerr=DB_SE, label="D_B (15 ppm)", capsize=2)
    ax1.errorbar(depth, D_N, yerr=DN_SE, label="D_N", capsize=2)
    ax1.errorbar(depth, D_H, yerr=DH_SE, label="D_H (recoil)", capsize=2)
    ax1.errorbar(depth, tumer1, yerr=tumer2, label="Tumor D_B (52.5 ppm, upper bound)", capsize=2,linestyle="--")
    ax1.set_yscale("log")  
    ax1.set_ylabel("Dose per unit fluence (Gy cm^2)") 
    ax1.set_title("KUR beam, soft tissue, N = %d" % N)
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(True, which="both", alpha=0.3)
    ax2.errorbar(depth, ng, yerr=ng_se, label="(n,gamma) rate", capsize=2)
    ax2.set_yscale("log")
    ax2.set_ylabel("(n,gamma) rate (cm^-1)")
    ax2.set_xlabel("Depth (cm)")
    ax2.legend()
    ax2.grid(True, which="both", alpha=0.3)
    # 8. 存图：stem + ".png"；show
    fig.tight_layout()   
    fig.savefig(stem + ".png", dpi=150)
    plt.show()
