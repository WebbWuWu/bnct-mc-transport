from lesson11_water import build_water
from lesson07_freegas import run_slab
from lesson08_dose import kur_table,KEEP_FACTOR,save_result,D,NL,OUT_DIR
import os,time

N_hi=20000
seeds=[2001+i for i in range (30)]

if __name__=="__main__":
    edges,shares=kur_table()
    E_keep=KEEP_FACTOR*edges[-1]
    nuclides,order=build_water(E_keep)
    for seed in seeds:
        path=os.path.join(OUT_DIR,"water_N%d_s%d"%(N_hi,seed))
        if os.path.exists(path+".json"):
            print("seed=%d已经有了，跳过"%seed)
            continue
        t1=time.time()
        r=run_slab(N_hi, seed, D, nuclides, None, order, z_src=0.0, implicit=True, NL=NL, src=(edges, shares))
        save_result(r,path)
        t2=time.time()
        print("花费时间为%.6f"%(t2-t1))
        