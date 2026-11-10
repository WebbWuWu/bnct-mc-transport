# ---- 路径（2026-11-10 整理文件夹，Claude 加）：本文件在 my_code 的子文件夹里。
#      下面几行把 my_code 根目录和各子文件夹加进搜索路径，后面的 import 才找得到别的文件。----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

from water import build_water
from transport import run_slab
from dose import kur_table,KEEP_FACTOR,save_result,D,NL,OUT_DIR
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
        