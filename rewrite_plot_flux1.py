# ============================================================
#  rewrite_plot_flux1.py —— 默写 误差棒通量图（第 1 次）
#  2026-10-14 Claude 开文件（用户 21:25 要求），只写了这段注释；下面全部由你写
#  放在 my_code 根目录：程序按 HERE 找三个 flux_case*.csv
# ============================================================
#  要写：
#    前半直接默：import、HERE、read_flux_csv、CASES、ERR_SCALE
#    画图 12 行对着下面的提纲写，想不起来的空着写 # ?
#      ① 开一张图：拿到「整张图」和「一块坐标区」，尺寸 7×5
#      ② 对三种情况循环，每次取出路径和图例文字
#      ③ 用你上面的函数把这个 csv 读成三个列表
#      ④ 在坐标区上画带误差棒的线：横轴深度、纵轴通量、误差棒 = 每个 SE × 缩放系数、图例文字；
#         外观四项：圆点标记、点大小 3、端帽 2、线宽 1
#      ⑤ 横轴标题  ⑥ 纵轴标题  ⑦ 图标题  ⑧ 显示图例  ⑨ 网格（透明度 0.3）  ⑩ 纵轴改对数
#      ⑪ 存图（分辨率 150、裁掉白边）—— 🔴 文件名别用 lesson03_flux.png（在 git 里，会被覆盖）
#      ⑫ 弹窗显示
#    思路：先想主语 —— 开图、弹窗找 plt；存图找 fig；其余找 ax
#  不写：assert、print
#  外观参数的数字（7×5、3、2、1、0.3、150）写错不扣分（10/14 用户同意）
#  不看：copy_plot_flux.py、plot_lesson03_scaffold.py、讲解_误差棒通量图.md
#  写完：存盘，Ctrl+F5 能出图就出图；把文件名发给 Claude 对比
# ============================================================
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
    
    return z_mid,phi,err

CASES = [
    (os.path.join(HERE, "flux_case1.csv"), "case 1: d=2.0, $\\Sigma_t$=1.0, $\\Sigma_s$=0.0"),
    (os.path.join(HERE, "flux_case2.csv"), "case 2: d=2.0, $\\Sigma_t$=1.0, $\\Sigma_s$=0.8"),
    (os.path.join(HERE, "flux_case3.csv"), "case 3: d=5.0, $\\Sigma_t$=1.0, $\\Sigma_s$=0.9"),
]

ERR_SCALE=1.0
fig,ax=plt.subplots(figsize=(7,5))

for path,label in CASES:
    z,phi,err=read_flux_csv(path)
    ax.errorbar(x=z,y=phi,yerr=[e*ERR_SCALE for e in err],label=label,marker="o", markersize=3, capsize=2, linewidth=1)

ax.set_xlabel("depth (cm)")
ax.set_ylabel("flux (cm$^{-2}$ per source neutron)")
ax.set_title("Track-length flux in a 1D slab (N=1e6, seed=2026)")
ax.legend()
ax.grid(True,alpha=0.3)
ax.set_yscale("log")

fig.savefig(os.path.join(HERE, "copy_lesson03_flux.png"), dpi=150, bbox_inches="tight")
plt.show()
