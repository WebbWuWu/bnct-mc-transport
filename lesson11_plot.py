import json,os,math
import matplotlib.pyplot as plt
from mcstat import mean_se
from lesson08_dose import EV2J,G_PER_KG,D,NL,OUT_DIR
from lesson11_water import rho_w

n_hi=100000
my_path=os.path.join(OUT_DIR,"water_N100000_s1010.json")
omcnames=[os.path.join(OUT_DIR,"openmc_water_s%d.json"%i) for i in range(1,11)]
omcnames_iso=[os.path.join(OUT_DIR,"openmc_water_iso_s%d.json"%i) for i in range(1,11)]
dz=D/NL
coe=EV2J*G_PER_KG/rho_w/dz
nucs=["B10","H1","O16"]
tallies=["water H1 heating","water O16 heating","water flux","water B10 heating"]

# [Claude 判据] M2 定稿图对拍数（Claude 10/26 另写代码算的；层号从 0 数：0 = 0–0.5 cm，4 = 2–2.5 cm，10 = 5–5.5 cm，20 = 10–10.5 cm，29 = 14.5–15 cm）
#   抓：ρ 用成 1.06（差 6%）、核素配错、忘了除 N / 多除了 N、OpenMC 只读了一个种子、比值分子分母反了、各向同性和原版文件拿反了
REF = {
    "mine B10": {0: 1.490063e-12, 4: 2.954544e-12, 20: 5.590284e-13, 29: 6.167134e-14},
    "mine H1":  {0: 8.180827e-13, 4: 5.298876e-13, 20: 5.502440e-14, 29: 1.383295e-14},
    "mine O16": {0: 4.816720e-14, 4: 3.243521e-14, 20: 3.583942e-15, 29: 9.845119e-16},
    "omc B10":  {0: 1.490882e-12, 4: 2.958639e-12, 20: 5.600350e-13, 29: 6.100738e-14},
    "omc H1":   {0: 8.104074e-13, 4: 5.002720e-13, 20: 5.602501e-14, 29: 1.513968e-14},
    "omc O16":  {0: 4.775545e-14, 4: 3.046984e-14, 20: 3.607275e-15, 29: 1.015442e-15},
    "ratio B10": {0: 9.994504e-01, 4: 9.986161e-01, 20: 9.982027e-01, 29: 1.010883e+00},
    "ratio H1":  {0: 1.009471e+00, 4: 1.059199e+00, 20: 9.821398e-01, 29: 9.136883e-01},
    "ratio O16": {0: 1.008622e+00, 4: 1.064502e+00, 20: 9.935318e-01, 29: 9.695403e-01},
    "ratio_err H1": {0: 1.883996e-02, 4: 3.137634e-02, 20: 8.978728e-02, 29: 1.478110e-01},
    "iso H1":   {0: 1.013872e+00, 10: 9.795540e-01, 20: 9.486843e-01, 29: 9.181527e-01},
    "iso O16":  {0: 1.021101e+00, 10: 9.715534e-01, 20: 9.405091e-01, 29: 9.113590e-01},
    "iso flux": {0: 1.000251e+00, 10: 9.997894e-01, 20: 9.987825e-01, 29: 9.967708e-01},
    "iso B10":  {0: 1.000035e+00, 10: 9.999270e-01, 20: 9.993940e-01, 29: 9.986848e-01},
    "iso_err H1": {0: 7.650903e-04, 10: 1.884047e-03, 20: 3.925756e-03, 29: 6.381620e-03},
}

def check_ref(tag, values):
    """values = 你算出来的那条 30 个数的列表；和 REF[tag] 里那几层比，千分之一以内。"""
    for k, want in REF[tag].items():
        got = values[k]
        assert abs(got / want - 1) < 1e-3, \
            "对拍 %s 第 %d 层：你的 %.6e，Claude 另算 %.6e（差 %.2f%%）" % (tag, k, got, want, (got / want - 1) * 100)
    print("对拍过：%s" % tag)
    
def read(dic,nuclide):
    mean=[]
    se=[]
    for k in range(NL):
        su=dic["kerma_sum"][nuclide][k]
        sq=dic["kerma_sq"][nuclide][k]
        me,_,se1=mean_se(su,sq,n_hi)
        mean.append(me*coe)
        se.append(se1*coe)
    return mean,se

def retally(paths,talname,co):
    mean_seed=[]
    se_seed=[]
    for path in paths:
        with open(path,"r")as f:
            d=json.load(f)
            mean_seed.append(d[talname]["mean"])
            se_seed.append(d[talname]["se"])
    me1=[]
    se1=[]
    for i in range(NL):
        mean1=0.0
        se21=0.0
        for r in range(len(paths)):
            mean1+=mean_seed[r][i]
            se21+=se_seed[r][i]**2
        me1.append(mean1/len(paths)*co)
        se1.append(math.sqrt(se21)/len(paths)*co)
    return me1,se1

def errrtio(me1,se1,me2,se2):
    rat=[]
    rat_err=[]
    for i in range(NL):
        rat.append(me1[i]/me2[i])
        rat_err.append(rat[i]*math.sqrt((se1[i]/me1[i])**2+(se2[i]/me2[i])**2))
    return rat,rat_err

def plot(my_me,my_se,omc_me,omc_se,rat,rat_err):
    deep=[(i+0.5)*dz for i in range(NL)]
    fig1,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=(7,7))
    doses=["D_B", "D_H", "D_O"]
    for j in range (len(nucs)):
        ax1.plot(deep,omc_me[nucs[j]],color="C%d"%j,label="%s (%s) OpenMC"%(doses[j],nucs[j]))
        ax1.errorbar(deep,my_me[nucs[j]],yerr=my_se[nucs[j]],fmt="o",markersize=3,capsize=2,color="C%d"%j,label="%s (%s) this code"%(doses[j],nucs[j]))
        ax2.errorbar(deep,rat[nucs[j]],yerr=rat_err[nucs[j]],fmt="o",markersize=3,capsize=2,color="C%d"%j,label="%s (%s)"%(doses[j],nucs[j]))
    ax1.set_yscale("log")
    ax2.axhline(1,color="gray",linewidth=0.8)
    ax1.grid(True,which="both",alpha=0.3)
    ax2.grid(True,which="both",alpha=0.3)
    ax1.legend(fontsize=8,ncol=2,loc="lower left")
    ax2.legend(fontsize=8)
    ax1.set_title("KUR beam, water + 15 ppm B-10: this code (1e5) vs OpenMC (10 seeds x 1e7)", fontsize=9)
    ax1.set_ylabel("Kerma per unit fluence (Gy cm^2)")
    ax2.set_ylabel("this code / OpenMC")
    ax2.set_xlabel("Depth (cm)")
    ax2.set_title("neutron kerma only (MT=301); photon dose not included; deepest O16 point: R = 0.22 > 0.2, SE unreliable", fontsize=8)
    fig1.tight_layout()
    path1=os.path.join (OUT_DIR,"m2_final.png")
    fig1.savefig(path1,dpi=150)
    
def plot2(ratio,err):
    deep=[(i+0.5)*dz for i in range(NL)]
    fig2,bx=plt.subplots()
    for j in range(len(tallies)):
        bx.errorbar(deep,ratio[tallies[j]],yerr=err[tallies[j]],fmt="o",markersize=3,capsize=2,color="C%d"%j,label="%s"%tallies[j])
    bx.axhline(1,color="gray",linewidth=0.8)
    bx.set_xlabel("Depth (cm)")
    bx.set_ylabel("isotropic / original")
    bx.set_title("Plan B: OpenMC elastic CM-isotropic / original (10 seeds each)")
    bx.grid(True,which="both",alpha=0.3)
    bx.legend(fontsize=8)
    fig2.tight_layout()
    path2=os.path.join (OUT_DIR,"m2_planB.png")
    fig2.savefig(path2,dpi=150)
    
if __name__=="__main__":
    with open(my_path,"r")as f:
        d1=json.load(f)
    my_mean = {}
    my_se = {}
    omc_mean = {}
    omc_se = {}
    ratio_mean = {}
    ratio_err = {}
    for nuc in nucs:
        my_mean[nuc],my_se[nuc]=read(d1,nuc)
        omc_mean[nuc],omc_se[nuc]=retally(omcnames,"water %s heating"%nuc,coe)
        ratio_mean[nuc],ratio_err[nuc]=errrtio(my_mean[nuc],my_se[nuc],omc_mean[nuc],omc_se[nuc])
    ratio2={}
    err2={}
    for i in range(len(tallies)):
        me1,se1=retally(omcnames_iso,tallies[i],1.0)
        me2,se2=retally(omcnames,tallies[i],1.0)
        a,b=errrtio(me1,se1,me2,se2)
        ratio2[tallies[i]]=a
        err2[tallies[i]]=b
    plot(my_mean,my_se,omc_mean,omc_se,ratio_mean,ratio_err)  
    plot2(ratio2,err2)
    plt.show()
