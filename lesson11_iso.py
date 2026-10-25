import openmc,os,json,h5py
from lesson09_openmc import geom,kur,mesh,run,tallies,path1
from lesson08_dose import M_N,OUT_DIR,NL
from lesson11_water import rho_w,share_B,nucs_w
import numpy as np
from lesson11_openmc import WATER_TALLIES,WATER_TOTAL_TALLIES,water_tally,water_total_tally,mater_water,fold

ISO_DIR=os.path.join(path1,"iso_lib")
os.makedirs(ISO_DIR,exist_ok=True)
XML_ISO=os.path.join(ISO_DIR,"iso.xml")

def make_iso(name):
    oripath=os.path.join(fold,name+".h5")
    na=openmc.data.IncidentNeutron.from_hdf5(oripath)
    el=na.reactions[2]
    d=el.products[0].distribution[0]
    E_ol=d.angle.energy[0]
    E_oh=d.angle.energy[-1]
    g=openmc.stats.Uniform(-1,1)
    two=openmc.data.AngleDistribution([E_ol,E_oh],[g,g])
    d.angle=two
    path=os.path.join(ISO_DIR,name+".h5")
    na.export_to_hdf5(path,mode="w")
    return path

def make_lib(paths):
    lib=openmc.data.DataLibrary()
    for i in paths:
        lib.register_file(i)
    lib.export_to_xml(XML_ISO)
    return XML_ISO

if __name__=="__main__":
    lis=[]
    for i in nucs_w:
        lis.append(make_iso(i))
    make_lib(lis)
    openmc.config["cross_sections"]=XML_ISO
    # ---- Claude 判据（方案乙）：生成新库之后、开跑 OpenMC 之前 ----
    # ① 抓：换了没存上 / 属性名拼错（d.angel）悄悄新建 / 区间写错。新文件里每个核素每个能点 μ̄_cm 都必须是 0
    for nuc in nucs_w:
        new = openmc.data.IncidentNeutron.from_hdf5(os.path.join(ISO_DIR, nuc + ".h5"))
        ang = new.reactions[2].products[0].distribution[0].angle
        worst = max(abs(m.mean()) for m in ang.mu)
        assert worst < 1e-12, "%s：新文件里弹性 μ̄_cm 最大 %.3g，没换成质心系各向同性" % (nuc, worst)
        lo, hi = ang.mu[0].x[0], ang.mu[0].x[-1]
        assert lo == -1.0 and hi == 1.0, "%s：μ 的区间是 [%g, %g]，余弦应该覆盖 [-1, 1]" % (nuc, lo, hi)
    # ② 抓：原文件本来就各向同性（那方案乙改了等于没改）/ 读错文件。真实 O16 在最高能处强烈前冲
    old = openmc.data.IncidentNeutron.from_hdf5(os.path.join(fold, "O16.h5"))
    ang0 = old.reactions[2].products[0].distribution[0].angle
    print("原库 O16 弹性 μ̄_cm：最低能 %.3g eV 处 %.4f；最高能 %.3g MeV 处 %.4f"
          % (ang0.energy[0], ang0.mu[0].mean(), ang0.energy[-1] / 1e6, ang0.mu[-1].mean()))
    assert ang0.mu[-1].mean() > 0.1, "原库 O16 最高能处 μ̄_cm = %.4f，不像真实的前冲分布" % ang0.mu[-1].mean()
    # ③ 抓：导出时截面或加热系数跟着变了（方案乙只许改角度）。O16 的 MT=2、MT=301 在 294 K 必须逐位相同
    new_O = openmc.data.IncidentNeutron.from_hdf5(os.path.join(ISO_DIR, "O16.h5"))
    for mt in (2, 301):
        a = old.reactions[mt].xs["294K"]
        b = new_O.reactions[mt].xs["294K"]
        assert np.array_equal(a.x, b.x) and np.array_equal(a.y, b.y), "O16 MT=%d 的 294 K 表导出后变了" % mt
    # ④ 抓：新清单里 O16 还指向原库 / 切换那一行被 import 盖回去了
    lib = openmc.data.DataLibrary.from_xml(XML_ISO)
    for nuc in nucs_w:
        p = lib.get_by_material(nuc)["path"]
        assert os.path.samefile(p, os.path.join(ISO_DIR, nuc + ".h5")), "新清单里 %s 指向 %s，不是新库" % (nuc, p)
    assert os.path.samefile(str(openmc.config["cross_sections"]), XML_ISO), \
        "openmc.config 现在是 %s，没切到新清单（切换那一行要在所有 import 之后）" % openmc.config["cross_sections"]
    print("方案乙判据 ①–④ 全过：三个核素弹性已改成质心系各向同性，截面和 MT=301 没动，OpenMC 会读新库")
    mat=mater_water()
    for seed in range(1, 11):
            disp=kur()
            tal2=water_tally(WATER_TALLIES)
            for t in water_total_tally(WATER_TOTAL_TALLIES):     # 两个集合合成一个：run() 只收一个
                tal2.append(t)
            path4=os.path.join(path1,"water_iso_s%d"%seed)
            os.makedirs(path4,exist_ok=True)
            stat_path=run(disp,tal2,path4,batch=100,n_bat=100000,seed=seed,mat=mat)
            dirc3={}
            for a,b,c in WATER_TALLIES+WATER_TOTAL_TALLIES:
                m,s=tallies(stat_path,a)
                dirc3[a]={"mean":m,"se":s}
            dirc3["batches"]=100
            dirc3["particles"]=100000
            dirc3["seed"]=seed
            path7=os.path.join(OUT_DIR,"openmc_water_iso_s%d.json"%seed)
            with open(path7,"w") as f:
                json.dump(dirc3,f)
            print("种子 %d 存好：%s" % (seed, path7))
