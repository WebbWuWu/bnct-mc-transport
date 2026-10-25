import numpy as np
import math,os

HERE = os.path.dirname(os.path.abspath(__file__))
path=os.path.join(HERE,'data','photon','water_mu_muen_nist.txt')

data=np.loadtxt(path,encoding="utf-8")

print(data.shape)

energy=data[:,0]
mu_rho=data[:,1]
muen_rho=data[:,2]

def loglog(x,y,e):
    lo=np.searchsorted(x,e,side="right")-1
    frac=(math.log(e)-math.log(x[lo]))/(math.log(x[lo+1])-math.log(x[lo]))
    lny=math.log(y[lo])+frac*(math.log(y[lo+1])-math.log(y[lo]))
    yv=math.exp(lny)
    return yv

if __name__=="__main__":
    E_lines=[2.2233,0.4778]
    for e in E_lines:
        mu=loglog(energy,mu_rho,e)
        mu_en=loglog(energy,muen_rho,e)
        ratio=mu_en/mu
        lam=1/mu/0.998
        print("e=%.4g,mu/rho=%.4g,muen/rho=%.4g,ratio=%.4g,lambda=%.4g"%(e,mu,mu_en,ratio,lam))
