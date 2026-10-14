import os
import h5py
import numpy as np

LIB=r"E:\BNCT_GPU\lib80x_hdf5"
T="294K"
E_TEST=0.0253

def total_xs(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    assert os.path.exists(path), "文件不在这里：%s" % path
    
    with h5py.File(path,"r")as f:
        root=f[nuc]
        E=root["energy"][T][:]
        awr=float(root.attrs["atomic_weight_ratio"])
        names=list(root["reactions"].keys())
        
        total=np.zeros(len(E))
        at_i0={}
        i0=int(np.searchsorted(E,E_TEST,side="right"))-1
        n_add=0
        
        for name in sorted(names):
            rx=root["reactions"][name]
            mt=rx.attrs["mt"]
            redundant=rx.attrs["redundant"]
            
            xs_ds=rx[T]["xs"]
            thr=int(xs_ds.attrs["threshold_idx"])
            xs=xs_ds[:]
            
            at_i0[int(mt)]=0.0 if i0<thr else xs[i0-thr]
            
            if redundant:
                continue
            
            total[thr:]+=xs
            n_add+=1
    print("[%s] %d 个道，加了 %d 个，跳过 %d 个冗余道"
          % (nuc, len(names), n_add, len(names) - n_add))

    assert np.all(total > 0), "总截面出现非正值，最小 %.3e" % np.min(total)          
    return E,total,at_i0,i0,awr

if __name__ == "__main__":
    from xslib import total_xs as total_xs_orig
    for nuc in ["B10", "H1"]:
        E1, t1, a1, i1, w1 = total_xs(nuc, T)
        E2, t2, a2, i2, w2 = total_xs_orig(nuc, T)
        # [Claude assert] 抓：抄写稿和原文在网格、总截面、i0、awr 上有任何一位不同
        assert np.array_equal(E1, E2) and np.array_equal(t1, t2), "%s：E 或 total 和原文不逐位相同" % nuc
        assert i1 == i2 and w1 == w2, "%s：i0 %r vs %r，awr %r vs %r" % (nuc, i1, i2, w1, w2)
        # [Claude assert] 抓：at_i0 少记了道（比如冗余道没记）或记错了值
        assert a1 == a2, "%s：at_i0 不同，抄写稿 %d 个键、原文 %d 个键" % (nuc, len(a1), len(a2))
        print("%s 逐位相同 ✓   total[i0] = %.6f b" % (nuc, t1[i1]))
