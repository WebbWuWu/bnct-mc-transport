import math,os,h5py
from xs_lookup import LIB
import numpy as np
from lesson08_dose import T

def qth(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    gra=[]
    with h5py.File(path,"r")as f:
        root=f[nuc]
        E=root["energy"][T][:]
        awr=float(root.attrs["atomic_weight_ratio"])
        names=list(root["reactions"].keys())
        n=0
        for name in sorted(names):
            redundant=root["reactions"][name].attrs["redundant"]
            if redundant==1:
                continue
            mt=int(root["reactions"][name].attrs["mt"])
            Q_value=float(root["reactions"][name].attrs["Q_value"])
            if Q_value>=0:
                n+=1
                continue
            xs_ds=root["reactions"][name][T]["xs"]
            xs=xs_ds[:]
            thr=int(xs_ds.attrs["threshold_idx"])
             # [Claude assert] 抓：threshold_idx 和数据长度对不上（读错温度组，同 load_heating 那条）
            assert len(xs) + thr == len(E), \
            "%s MT=%d：截面 %d 个点 + 阈值下标 %d ≠ 网格 %d 个点" % (nuc, mt, len(xs), thr, len(E))
            E_th=-Q_value*(awr+1)/awr
            E_start=E[thr]
            for i in range(len(xs)):
                if xs[i]>0:
                    E_first=E[i+thr]
                    break
            else:
                print("%s核素MT=%d的道截面全为0"%(name,mt))
                continue
            ratio=E_start / E_th
            # [Claude assert] 抓：E_th 公式写反（A/(A+1)）或用错 A —— 物理上阈值以下不可能有反应，
            #   第一个截面 > 0 的点必须在 E_th 以上（容差 1e-9 只给浮点舍入）。
            #   两条路：库里截面真正开始不为 0 的位置  vs  Q 值 + 运动学。
            #   ⚠️ 写反成 A/(A+1) 时 E_th 变小，这条反而更容易过 —— 所以还有下面第 5 节的表要看。
            assert E_first >= E_th * (1.0 - 1e-9), \
                "%s MT=%d：第一个截面 > 0 的点 %.6e eV 在运动学阈值 %.6e eV 以下 —— 公式或读的量有问题" \
                % (nuc, mt, E_first, E_th)
            gra1=(mt,Q_value,E_th,E_start,E_first,ratio,)
            gra.append(gra1)
    print("\n%s（awr = %.6f）" % (nuc, awr))
    print("%6s %14s %14s %14s %14s %12s" % ("MT", "Q (eV)", "E_th (eV)", "E_start (eV)", "E_first (eV)", "E_start/E_th"))
    for mt, Q, E_th, E_start, E_first, ratio in gra:
        print("%6d %14.6e %14.6e %14.6e %14.6e %12.6f" % (mt, Q, E_th, E_start, E_first, ratio))
    print("Q ≥ 0 跳过了 %d 个道" % n)
    return gra
    
if __name__=="__main__":
    qth1=[0.0]*3
    n=0
    for i in ["H1", "O16", "B10"]:
        qth1[n]=qth(i,T)
        n+=1
