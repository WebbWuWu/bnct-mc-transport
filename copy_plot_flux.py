import matplotlib.pyplot as plt
import os

HERE=os.path.dirname(os.path.abspath(__file__))

def read_flux_csv(path):
    z_mid=[]
    phi=[]
    err=[]
    with open(path,"r",encoding="utf-8")as f:
        lines=f.readlines()
    for line in lines[1:]:
        parts=line.strip().split(",")
        z_mid.append(float(parts[3]))
        phi.append(float(parts[4]))
        err.append(float(parts[5]))
    assert len(z_mid) > 0, "csv 里一行数据都没读到（只有表头或文件为空）：%s  总行数=%d" % (path, len(lines))
    return z_mid,phi,err

CASES=[
    (os.path.join(HERE,"flux_case1.csv"),"case 1: d=2.0,$\\Sigma_t$=1.0,$\\Sigma_s$=0.0"),
    (os.path.join(HERE,"flux_case2.csv"),"case 2: d=2.0,$\\Sigma_t$=1.0,$\\Sigma_s$=0.8"),
    (os.path.join(HERE,"flux_case3.csv"),"case 3: d=5.0,$\\Sigma_t$=1.0,$\\Sigma_s$=0.9"),
]

ERR_SCALE=1.0

fig,ax=plt.subplots(figsize=(7,5))

for path,label in CASES:
    z,phi,err=read_flux_csv(path)
    ax.errorbar(x=z,y=phi,yerr=[e*ERR_SCALE for e in err],label=label,marker="o",markersize=3,capsize=2,linewidth=1)

ax.set_xlabel("depth(cm)")
ax.set_ylabel("flux(cm$^{-2}$ per source neutron)")
ax.set_title("track-length flux in a 1D slab (N=1e6,seed=2026)")
ax.legend()
ax.grid(True,alpha=0.3)
ax.set_yscale("log")

fig.savefig(os.path.join(HERE,"copy_lesson03_flux.png"),dpi=150,bbox_inches="tight")
plt.show()
