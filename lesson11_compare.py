import os,json,math
from lesson08_dose import OUT_DIR,NL
from mcstat import mean_se
from lesson09_compare import omc_run,compare_z
from lesson07_freegas import Q_val

N_HIST=100000
licom=[
    ("flux","track_sum","track_sq",None,"water flux",1.0),
    ("H1 heating","kerma_sum","kerma_sq","H1","water H1 heating",1.0),
    ("O16 heating","kerma_sum","kerma_sq","O16","water O16 heating",1.0),
    ("B10 heating","kerma_sum","kerma_sq","B10","water B10 heating",1.0),
    ("H1 (n,γ)","gH_sum","gH_sq",None,"water H1 n,gamma",1.0),
    ("B10 (n,α)","Bc_sum","Bc_sq",None,"water B10 n,a",Q_val),
]

def my_run2(su,sq,his):
    mean=[]
    SE=[]
    for i in range(NL):
        mean1,_,SE1=mean_se(su[i],sq[i],his)
        mean.append(mean1)
        SE.append(SE1)
    return mean,SE

if __name__ == "__main__":
    path1=os.path.join(OUT_DIR,"water_N100000_s1010.json")
    with open(path1, "r") as f:
        book1=json.load(f)
    path2=os.path.join(OUT_DIR,"openmc_water_s1.json")
    with open(path2, "r") as f:
        book2=json.load(f)    
    dirc={}
    print("%-16s %-8s %-8s %-8s %-8s"%("量","1si内","2si内","max_z","ka2"))
    for a,b,c,d,e,g in licom:
        if d==None:
            su=book1[b]
            sq=book1[c]
        else:
            su=book1[b][d]
            sq=book1[c][d]
        mean,SE=my_run2(su,sq,N_HIST)
        ommean,omse=omc_run(book2,e,g)
        z,n1,n2,z_ma,ka2=compare_z(mean,SE,ommean,omse)
        print("%-16s %-8s %-8s %-8s %-8s"%(a,n1,n2,z_ma,ka2))
        dirc[a]={
            "z":z,
            "n1":n1,
            "n2":n2,
            "zmax":z_ma,
            "ka2":ka2
        }
    path3=os.path.join(OUT_DIR,"compare_water.json")
    with open(path3,"w")as f:
        json.dump(dirc,f)
