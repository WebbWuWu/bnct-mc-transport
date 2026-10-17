from lesson09_compare import Q_B,Q_N,OUT_DIR,my_run,omc_run,merge_two,N_HIST
import matplotlib.pyplot as plt
import json
from lesson08_dose import D,NL,EV2J,G_PER_KG,RHO
import os
import math
from mcstat import mean_se

E_R=1315.8            #check4 里 4-1 打印的「路 1 中位数」
dz=D/NL
c=EV2J * G_PER_KG / (RHO * dz)
nuc_json=["dose_N100000_s1010_nuc.json",
          "dose_N100000_s1200_nuc.json"]
kur_json=[
    "openmc_kur.json",
    "openmc_kur_s2.json",
    "openmc_kur_s5.json",
    "openmc_kur_s6.json",
    "openmc_kur_s7.json",
    "openmc_kur_s8.json",
    "openmc_kur_s9.json",
    "openmc_kur_s10.json",
    "openmc_kur_s11.json",
    "openmc_kur_s12.json"
]
s_ab=[
    "openmc_kur_sab3.json",
    "openmc_kur_sab4.json"
]
sapa="lesson09_final.png"
save_path=os.path.join(OUT_DIR,sapa)

def merge_n(mean,se):
    n=len(mean)
    res1=[]
    res2=[]
    for i in range(NL):
        mean_n=0.0
        se2_n=0.0
        for r in range(n):
            mean_n+=mean[r][i]
            se2_n+=se[r][i]**2
        res1.append(mean_n/n)
        res2.append(math.sqrt(se2_n)/n)
    return res1,res2

out_class=[]
omc_class=[]
sab_class=[]
for name in nuc_json:
    p=os.path.join(OUT_DIR,name)
    with open(p,"r")as f:
        d=json.load(f)
        out_class.append(d)
for name in kur_json:
    p=os.path.join(OUT_DIR,name)
    with open(p,"r")as f:
        d=json.load(f)
        omc_class.append(d)
for name in s_ab:
    p=os.path.join(OUT_DIR,name)
    with open(p,"r")as f:
        d=json.load(f)
        sab_class.append(d)

DBmean=[]
DBse=[]
DNmean=[]
DNse=[]
DHmean=[]
DHse=[]
for name in out_class:
    DBmean1,DBse1=my_run(name,"Bc",N_HIST)
    DBmean.append(DBmean1)
    DBse.append(DBse1)
    DNmean1,DNse1=my_run(name,"Nc",N_HIST)
    DNmean.append(DNmean1)
    DNse.append(DNse1)
    DHx=[]
    DHx2=[]
    for i in range(NL):
        mean1,_,se1=mean_se(name["sigma_x"]["H1"][i],name["sigma_x2"]["H1"][i],N_HIST)
        DHx.append(mean1)
        DHx2.append(se1)
    DHmean.append(DHx)
    DHse.append(DHx2)
DBmean_nl,DBse_nl=merge_n(DBmean,DBse)
DNmean_nl,DNse_nl=merge_n(DNmean,DNse)
DHmean_nl,DHse_nl=merge_n(DHmean,DHse)

omcDBmean=[]
omcDBse=[]
omcDNmean=[]
omcDNse=[]
omcDHmean=[]
omcDHse=[]
for name in omc_class:
    omcDBmean1,omcDBse1=omc_run(name,"kur B10 n,a",Q_B)
    omcDBmean.append(omcDBmean1)
    omcDBse.append(omcDBse1)
    omcDNmean1,omcDNse1=omc_run(name,"kur N14 n,p",Q_N)
    omcDNmean.append(omcDNmean1)
    omcDNse.append(omcDNse1)
    omcDHx=[]
    omcDHx2=[]
    omcDHmeanhe,omcDHsehe=omc_run(name,"kur H1 heating", 1.0)
    omcDHmeanR,omcDHseR=omc_run(name,"kur n,gamma", E_R)
    for i in range(NL):
        omcDHx.append(omcDHmeanhe[i]-omcDHmeanR[i])
        omcDHx2.append(math.sqrt(omcDHsehe[i]**2+omcDHseR[i]**2))
    omcDHmean.append(omcDHx)
    omcDHse.append(omcDHx2)
omcDBmean_nl,omcDBse_nl=merge_n(omcDBmean,omcDBse)
omcDNmean_nl,omcDNse_nl=merge_n(omcDNmean,omcDNse)
omcDHmean_nl,omcDHse_nl=merge_n(omcDHmean,omcDHse)

SabDBmean=[]
SabDBse=[]
for name in sab_class:
    SabDBmean1,SabDBse1=omc_run(name,"kur B10 n,a",Q_B)
    SabDBmean.append(SabDBmean1)
    SabDBse.append(SabDBse1)
SabDBmean_nl,SabDBse_nl=merge_n(SabDBmean,SabDBse)

for i in range(NL):
    DBmean_nl[i]=c*DBmean_nl[i]
    DBse_nl[i]=c*DBse_nl[i]
    DNmean_nl[i]=c*DNmean_nl[i]
    DNse_nl[i]=c*DNse_nl[i]
    DHmean_nl[i]=c*DHmean_nl[i]
    DHse_nl[i]=c*DHse_nl[i]
    omcDBmean_nl[i]=c*omcDBmean_nl[i]
    omcDBse_nl[i]=c*omcDBse_nl[i]
    omcDNmean_nl[i]=c*omcDNmean_nl[i]
    omcDNse_nl[i]=c*omcDNse_nl[i]
    omcDHmean_nl[i]=c*omcDHmean_nl[i]
    omcDHse_nl[i]=c*omcDHse_nl[i]
    SabDBmean_nl[i]=c*SabDBmean_nl[i]
    SabDBse_nl[i]=c*SabDBse_nl[i]
    
rat_DBmean=[0.0]*NL
rat_DBse=[0.0]*NL
rat_DNmean=[0.0]*NL
rat_DNse=[0.0]*NL
rat_DHmean=[0.0]*NL
rat_DHse=[0.0]*NL
for i in range(NL):
    rat_DBmean[i]=DBmean_nl[i]/omcDBmean_nl[i]
    rat_DBse[i]=rat_DBmean[i]*math.sqrt((DBse_nl[i]/DBmean_nl[i])**2+(omcDBse_nl[i]/omcDBmean_nl[i])**2)
    rat_DNmean[i]=DNmean_nl[i]/omcDNmean_nl[i]
    rat_DNse[i]=rat_DNmean[i]*math.sqrt((DNse_nl[i]/DNmean_nl[i])**2+(omcDNse_nl[i]/omcDNmean_nl[i])**2)
    rat_DHmean[i]=DHmean_nl[i]/omcDHmean_nl[i]
    rat_DHse[i]=rat_DHmean[i]*math.sqrt((DHse_nl[i]/DHmean_nl[i])**2+(omcDHse_nl[i]/omcDHmean_nl[i])**2)

fig,(ax1,ax2)=plt.subplots(2,1,sharex=True)
dep=[]
for i in range(NL):
    dep.append((i+0.5)*dz)
ax1.plot(dep,omcDBmean_nl,color="C0",label="D_B OpenMC")
ax1.plot(dep,omcDNmean_nl,color="C1",label="D_N OpenMC")
ax1.plot(dep,omcDHmean_nl,color="C2",label="D_H OpenMC")
ax1.plot(dep,SabDBmean_nl,color="C0",linestyle="--",label="D_B OpenMC S(a,b)")
ax1.errorbar(dep,DBmean_nl,yerr=DBse_nl,fmt="o",markersize=3,capsize=2,color="C0",label="D_B this code")
ax1.errorbar(dep,DNmean_nl,yerr=DNse_nl,fmt="o",markersize=3,capsize=2,color="C1",label="D_N this code")
ax1.errorbar(dep,DHmean_nl,yerr=DHse_nl,fmt="o",markersize=3,capsize=2,color="C2",label="D_H(H1) this code")
ax2.errorbar(dep,rat_DBmean,yerr=rat_DBse,fmt="o",markersize=3,capsize=2,color="C0",label="D_B")
ax2.errorbar(dep,rat_DNmean,yerr=rat_DNse,fmt="o",markersize=3,capsize=2,color="C1",label="D_N")
ax2.errorbar(dep,rat_DHmean,yerr=rat_DHse,fmt="o",markersize=3,capsize=2,color="C2",label="D_H(H1)")
ax2.axhline(1.0, color="gray", linewidth=0.8)
ax1.grid(True,which="both",alpha=0.3)
ax2.grid(True,which="both",alpha=0.3)
ax1.set_yscale("log")
ax2.set_ylim(0.7,1.3)
ax1.legend(fontsize=8)
ax2.legend(fontsize=8)
ax2.set_xlabel("Depth (cm)")
ax1.set_ylabel("Dose per unit fluence (Gy cm^2)")
ax1.set_title("KUR beam, soft tissue: this code vs OpenMC")
fig.tight_layout()
ax2.set_ylabel("this code / OpenMC")
fig.savefig(save_path,dpi=150)
plt.show()
