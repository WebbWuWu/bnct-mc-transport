import math
import os
from  rng import lcg, isotropic_direction
Max=10000
NL=20

def run_slab(N,seed,d,sigma_t,sigma_s):
    rng=lcg(seed)
    n_T=0
    n_A=0
    n_i=0
    n_R=0
    n_stuck=0
    track_sum=[0.0]*NL
    track_sq=[0.0]*NL
    v_sum=0
    v_sq=0
    dz=d/NL
    for i in range(N):
        track_h=[0.0]*NL
        n_hist=0
        z=0.0
        w=1.0
        for step in range(Max):
            s=-math.log(1-rng.random())/sigma_t
            z1=z+s*w
            if z1>d:
                z_end=d
            elif z1<0:
                z_end=0
            else:
                z_end=z1
            k1=int(z_end/dz)
            if k1>=NL:
                k1=NL-1
            if w==0.0:
                track_h[k1]+=s
            else:
                cur=z
                k=int(cur/dz)
                while True:
                    if w<0:
                        k_o=k*dz
                        z_step=max(k_o,z_end)
                    if w>0:
                        k_o=(k+1)*dz
                        z_step=min(k_o,z_end)
                    track_h[k]+=abs(z_step-cur)/abs(w)
                    if k==k1:
                        break
                    if k<k1:
                        k+=1
                    if k>k1:
                        k-=1
                    cur=z_step
            z=z1
            if z<0:
                n_R+=1
                break
            elif z>d:
                n_T+=1
                break
            elif rng.random()<sigma_s/sigma_t:
                n_i+=1
                n_hist+=1
                u,v,w=isotropic_direction(rng)
            else:
                n_hist+=1
                n_i+=1
                n_A+=1
                break
        else:
            n_stuck+=1
        for j in range(NL):
            track_sum[j]+=track_h[j]
            track_sq[j]+=track_h[j]**2
        v_sum+=sigma_t*sum(track_h)-n_hist
        v_sq+=(sigma_t*sum(track_h)-n_hist)**2
    return n_T/N,n_R/N,n_A/N,n_i/N,track_sum,track_sq,dz,v_sum,v_sq,sigma_t
       
        