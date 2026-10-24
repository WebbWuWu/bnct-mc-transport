import openmc,os,json,h5py
from lesson09_openmc import geom,kur,mesh,run,tallies,path1
from lesson08_dose import M_N,OUT_DIR,NL
from lesson11_water import rho_w,share_B,nucs_w

fold="/mnt/e/BNCT_GPU/lib80x_hdf5"

# ---- 10/24 Claude 写（用户 18:38 要求）：两张 tally 清单，格式同 lesson09_openmc.py 的 KUR_TALLIES：(名字, 核素列表, score) ----
#   30 层（带 mesh 过滤器）6 个：右边注释 = 对你 run_slab 返回的哪个键
WATER_TALLIES = [
    ("water flux",          [],      "flux"),        # track_sum
    ("water H1 heating",    ["H1"],  "heating"),     # kerma_sum["H1"]
    ("water O16 heating",   ["O16"], "heating"),     # kerma_sum["O16"]
    ("water B10 heating",   ["B10"], "heating"),     # kerma_sum["B10"]
    ("water H1 n,gamma",    ["H1"],  "(n,gamma)"),   # gH_sum
    ("water B10 n,a",       ["B10"], "(n,a)"),       # Bc_sum / Q_val
]
#   不分层（不加过滤器，整块板一个数；判据 3 方案丙的「总量」）4 个
WATER_TOTAL_TALLIES = [
    ("water total flux",       [],      "flux"),     # track_tot_sum
    ("water total H1 heating", ["H1"],  "heating"),  # kerma_tot_sum["H1"]
    ("water total O16 heating",["O16"], "heating"),  # kerma_tot_sum["O16"]
    ("water total B10 heating",["B10"], "heating"),  # kerma_tot_sum["B10"]
]

def water_tally(lst):
    fil=mesh()
    ts=openmc.Tallies()
    for a,b,c in lst:
        t_a=openmc.Tally(name=a)
        t_a.filters=[fil]
        if len(b)!=0:
            t_a.nuclides=b
        t_a.scores=[c]
        ts.append(t_a)
    return ts

def water_total_tally(lst):
    ts=openmc.Tallies()
    for a,b,c in lst:
        t_a=openmc.Tally(name=a)
        if len(b)!=0:
            t_a.nuclides=b
        t_a.scores=[c]
        ts.append(t_a)
    return ts

def mater_water():
    nuc_awr={}
    for name in nucs_w:
        path=os.path.join(fold,name+".h5")
        with h5py.File(path,"r")as f:
            root=f[name]
            awr=float(root.attrs["atomic_weight_ratio"])
            nuc_awr[name]=awr
    M_H=M_N*nuc_awr["H1"]
    M_O=M_N*nuc_awr["O16"]
    M_water=2*M_H+M_O
    share_H=M_H*2/M_water*(1-share_B)
    share_O=M_O/M_water*(1-share_B)
    frac={
        "H1":share_H,
        "O16":share_O,
        "B10":share_B
    }
    wat=openmc.Material(name='water')
    wat.add_nuclide('H1',share_H,'wo')
    wat.add_nuclide('O16',share_O,'wo')
    wat.add_nuclide('B10',share_B,'wo')
    wat.set_density('g/cm3', rho_w)
    wat.temperature=294
    return wat

if __name__ == "__main__":
    mat=mater_water()
        # [Claude assert] 抓：份额配错核素、漏乘 (1 − 硼份额)、密度单位错、awr 读错文件
    #   两条路：OpenMC 用它自己的原子质量表把 'wo' 换成数密度（atom/b·cm）
    #           vs  你 build_water 用 awr × M_N 算的（10/23 打印：N_H 6.6739e22、N_O 3.3369e22、N_B 9.0035e17 /cm³，4–5 位有效数字）
    #   两边原子质量差约 1e-7（文档 2.3 节），所以 N_H/N_O 容差 1e-5、和 10/23 打印值比容差 1e-4
    nd = mat.get_nuclide_atom_densities()
    assert abs(nd["H1"] / nd["O16"] - 2.0) < 1e-5, "N_H / N_O = %.9g，水里应该是 2" % (nd["H1"] / nd["O16"])
    for nuc, ref in (("H1", 6.6739e22), ("O16", 3.3369e22), ("B10", 9.0035e17)):
        assert abs(nd[nuc] * 1e24 / ref - 1) < 1e-4, \
            "%s：OpenMC 换算出 %.5e /cm³，build_water 是 %.5e /cm³" % (nuc, nd[nuc] * 1e24, ref)
    print("材料过：OpenMC 的数密度和 build_water 对上（N_H/N_O = %.9f）" % (nd["H1"] / nd["O16"]))
    # ---- 10/24 Claude 改（用户 19:16 粘贴 lesson09 第 177–192 行后要求代改）：10 个种子 × 1000 万中子 ----
    for seed in range(1, 11):
        disp=kur()
        tal2=water_tally(WATER_TALLIES)
        for t in water_total_tally(WATER_TOTAL_TALLIES):     # 两个集合合成一个：run() 只收一个
            tal2.append(t)
        path4=os.path.join(path1,"water_s%d"%seed)
        os.makedirs(path4,exist_ok=True)
        stat_path=run(disp,tal2,path4,batch=100,n_bat=100000,seed=seed,mat=mat)
        dirc3={}
        for a,b,c in WATER_TALLIES+WATER_TOTAL_TALLIES:
            m,s=tallies(stat_path,a)
            dirc3[a]={"mean":m,"se":s}
        dirc3["batches"]=100
        dirc3["particles"]=100000
        dirc3["seed"]=seed
        path7=os.path.join(OUT_DIR,"openmc_water_s%d.json"%seed)
        with open(path7,"w") as f:
            json.dump(dirc3,f)
        print("种子 %d 存好：%s" % (seed, path7))    
