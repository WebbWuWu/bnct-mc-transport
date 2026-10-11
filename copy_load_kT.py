import os
from xs_lookup import LIB,sigma_at
import h5py
import math

def load_kT(nuc,T):
    path=os.path.join(LIB,nuc+".h5")
    assert os.path.exists(path),"文件不在这里:%s"%path
    with h5py.File(path,"r")as f:
        root=f[nuc]
        kT=float(root["kTs"][T][()])
        T_label=float(T[:-1])
        assert abs(kT / (8.617333262e-5 * T_label) - 1.0) < 5e-3, \
            "kT=%r eV 与温度标签 %s 对不上：按标签应约为 %.6f eV，差超过 0.5%%" % (kT, T, 8.617333262e-5 * T_label)
    return kT
    
def sigma_at_ext(E,XS,e):
    assert e > 0.0 and math.isfinite(e), "要查的能量 e=%r 不是正的有限数（上游算错了），不能外推" % (e,)
    if e>=E[0]:
            sig,total=sigma_at(E,XS,e)
    else:
        sig1,total1=sigma_at(E,XS,E[0])
        sig=[0.0]*len(sig1)
        for i in range(len(sig1)):
            sig[i]=sig1[i]*math.sqrt(E[0]/e)
        total=total1*math.sqrt(E[0]/e)
    return sig,total

if __name__ == "__main__":
    print(load_kT("H1", "294K"))
    print(load_kT("O16", "294K"))
    E_f = [1e-5, 2e7]; XS_f = [[0.8, 0.8], [0.2, 0.2]]      # 第 6 课那个假库
    print(sigma_at_ext(E_f, XS_f, 3.1e-6))    # 低于表下限，外推
    print(sigma_at_ext(E_f, XS_f, 0.37))      # 在表里，不外推
