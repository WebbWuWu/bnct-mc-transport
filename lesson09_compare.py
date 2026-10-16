from lesson08_dose import OUT_DIR,NL
from mcstat import mean_se
import math
import os
import json

Q_B = 2341569.0
Q_N = 625876.1
N_HIST=100000
CMP_TABLE=[
    ("track","kur flux",1.0),
    ("Bc","kur B10 n,a",Q_B),
    ("Nc","kur N14 n,p",Q_N),
    ("ag","kur n,gamma",1.0)
]

def my_run(dirc,name,his):
    mean=[]
    se=[]
    su=dirc[name+"_sum"]
    sq=dirc[name+"_sq"]
    for i in range(NL):
        m,v,s=mean_se(su[i],sq[i],his)
        mean.append(m)
        se.append(s)
    return mean,se

def omc_run(dirc,name,coe):
    mean=dirc[name]["mean"]
    se=dirc[name]["se"]
    mean1=[]
    se1=[]
    for i in range(NL):
        mean1.append(mean[i]*coe)
        se1.append(se[i]*coe)
    return mean1,se1
    

def merge_two(mean1,se1,mean2,se2):
    mer1=[]
    mer2=[]
    for i in range(NL):
        mer1.append((mean1[i]+mean2[i])/2)
        mer2.append(math.sqrt(se1[i]**2+se2[i]**2)/2)
    return mer1,mer2

def compare_z(mean,se1,mc_mean,se2):
    z=[]
    z2=[]
    n1=0
    n2=0
    ka2=0.0
    for i in range(NL):
        z1=(mean[i]-mc_mean[i])/math.sqrt(se1[i]**2+se2[i]**2)
        z.append(z1)
        z2.append(abs(z1))
        if abs(z1)<1:
            n1+=1
        if abs(z1)<2:
            n2+=1
        ka2+=z1**2
    return z,n1,n2,max(z2),ka2

if __name__=="__main__":
    path1=os.path.join(OUT_DIR,"dose_N100000_s1010.json")
    with open(path1, "r") as f:
        book1=json.load(f)
    path2=os.path.join(OUT_DIR,"dose_N100000_s1200.json")
    with open(path2, "r") as f:
        book2=json.load(f)
    path3=os.path.join(OUT_DIR,"openmc_kur.json")
    with open(path3, "r") as f:
        book3=json.load(f)
    path4=os.path.join(OUT_DIR,"openmc_kur_s2.json")
    with open(path4, "r") as f:
        book4=json.load(f)
    dirc={}
    print("%-16s %-8s %-8s %-8s %-8s"%("量","1si内","2si内","max_z","ka2"))
    for key,name,coe in CMP_TABLE:
        mm1,ms1=my_run(book1,key,N_HIST)
        mm2,ms2=my_run(book2,key,N_HIST)
        om1,os1=omc_run(book3,name,coe)
        om2,os2=omc_run(book4,name,coe)
        mer1,mer2=merge_two(mm1,ms1,mm2,ms2)
        mer3,mer4=merge_two(om1,os1,om2,os2)
        z,n1,n2,zmax,chi2=compare_z(mer1,mer2,mer3,mer4)
        print("%-16s %8d %8d %8.2f %8.1f"%(name,n1,n2,zmax,chi2))
        dirc[name]={
            "z":z,
            "n1":n1,
            "n2":n2,
            "zmax":zmax,
            "chi2":chi2
        }
    path5=os.path.join(OUT_DIR,"compare_kur.json")
    with open(path5,"w")as f:
        json.dump(dirc,f)
