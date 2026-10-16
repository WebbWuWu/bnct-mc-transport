from lesson08_dose import OUT_DIR,MASS_FRAC,RHO,T,D,kur_table,NL
import os
import openmc
import math
import json
import numpy as np

KUR_TALLIES=[
    ("kur flux",[],"flux"),
    ("kur B10 n,a",["B10"],"(n,a)"),
    ("kur N14 n,p",["N14"],"(n,p)"),
    ("kur n,gamma",[],"(n,gamma)"),
    ("kur H1 heating",["H1"],"heating"),
    ("kur C12 heating",["C12"],"heating"),
    ("kur O16 heating",["O16"],"heating")
]

openmc.config['cross_sections']="/mnt/e/BNCT_GPU/lib80x_hdf5/cross_sections.xml"
path1=os.path.join(OUT_DIR,'openmc')
os.makedirs(path1,exist_ok=True)

def mater():
    soft=openmc.Material(name='soft tissue')
    MASS_sum=0.0
    for nuc in MASS_FRAC.keys():
        soft.add_nuclide(nuc,MASS_FRAC[nuc],'wo')
        MASS_sum+=MASS_FRAC[nuc]
    soft.set_density('g/cm3',RHO*MASS_sum)
    soft.temperature=294#和 lesson08 的 T 一致
    return soft

def geom(mat):  
    xp1=openmc.XPlane(x0=7.3,boundary_type='reflective')
    xp2=openmc.XPlane(x0=-7.3,boundary_type='reflective')
    yp1=openmc.YPlane(y0=7.3,boundary_type='reflective')
    yp2=openmc.YPlane(y0=-7.3,boundary_type='reflective')
    zp1=openmc.ZPlane(z0=D,boundary_type='vacuum')
    zp2=openmc.ZPlane(z0=0,boundary_type='vacuum')

    bou=-xp1&+xp2&-yp1&+yp2&-zp1&+zp2

    c=openmc.Cell(fill=mat,region=bou)
    geo=openmc.Geometry([c])
    return geo

def kur():
    edges,shares=kur_table()

    M=50
    edges1=[0.0]*((len(edges)-1)*M+1)
    shares1=[0.0]*(len(shares)*M)
    dens=[0.0]*(len(shares)*M)
    edges1[0]=edges[0]
    for i in range(len(edges)-1):
        r=(edges[i+1]/edges[i])**(1/M)
        for m in range(i*M+1,(i+1)*M+1):
            edges1[m]=edges[i]*r**(m-i*M)
    for j in range(len(shares1)):
        shares1[j]=shares[int(j/M)]/M
        dens[j]=shares1[j]/(edges1[j+1]-edges1[j])

    dens.append(0.0)

    dsp=openmc.stats.Tabular(edges1,dens,interpolation="histogram")
    return dsp

def mesh():
    grid=openmc.RegularMesh()
    grid.dimension=[1,1,NL]
    grid.lower_left=[-7.3,-7.3,0]
    grid.upper_right=[7.3,7.3,D]
    filt=openmc.MeshFilter(grid)
    return filt

def tally():
    f_grid=mesh()
    f_col=openmc.CollisionFilter([0])
    t=openmc.Tally(name='j1 uncollided flux')
    t.filters=[f_grid,f_col]
    t.scores=['flux']
    ts=openmc.Tallies([t])
    return ts

def tru_tally():
    fil=mesh()
    ts=openmc.Tallies()
    for a,b,c in KUR_TALLIES:
        t_a=openmc.Tally(name=a)
        t_a.filters=[fil]
        if len(b)!=0:
            t_a.nuclides=b
        t_a.scores=[c]
        ts.append(t_a)
    return ts
            
def run(E_dis,tal,inp,batch,n_bat,mat,seed=1):
    os.makedirs(inp,exist_ok=True)
    material=openmc.Materials([mat])
    geo=geom(mat)
    pos=openmc.stats.Point((0,0,1e-6))
    dirn=openmc.stats.Monodirectional((0,0,1))
    S=openmc.IndependentSource(space=pos,angle=dirn,energy=E_dis)

    st=openmc.Settings()
    st.run_mode='fixed source'
    st.particles=n_bat
    st.batches=batch
    st.seed=seed
    st.source=S
    
    path2=os.path.join(inp,'materials.xml')
    path3=os.path.join(inp,'geometry.xml')
    path4=os.path.join(inp,'settings.xml')
    path5=os.path.join(inp,'tallies.xml')
    material.export_to_xml(path2)
    geo.export_to_xml(path3)
    st.export_to_xml(path4)
    tal.export_to_xml(path5)
    
    openmc.run(cwd=inp)
    path6=os.path.join(inp,'statepoint.'+str(batch)+'.h5')
    
    return path6

def tallies(path,name1):
    with openmc.StatePoint(path)as sp:
        t2=sp.get_tally(name=name1)
        m=t2.mean
        s=t2.std_dev
        m1=m.flatten()
        m2=m1.tolist()
        s1=s.flatten()
        s2=s1.tolist()
    return m2,s2

if __name__ == "__main__":
    run_as1=False
    run_Tru=False
    run_Tru2=False
    run_Sab=True
    if run_Tru==True:
        disp=kur()
        tal2=tru_tally()
        path4=os.path.join(path1,"kur")
        os.makedirs(path4,exist_ok=True)
        mat=mater()
        stat_path=run(disp,tal2,path4,batch=100,n_bat=10000,mat=mat)
        dirc3={}
        for a,b,c in KUR_TALLIES:
            m,s=tallies(stat_path,a)
            dirc3[a]={"mean":m,"se":s}
        dirc3["batches"]=100
        dirc3["particles"]=10000
        path7=os.path.join(OUT_DIR,"openmc_kur.json")
        with open(path7,"w") as f:
            json.dump(dirc3,f)
    
    if run_Sab==True:
        for seeds in [3,4]:
            disp=kur()
            tal2=tru_tally()
            path4=os.path.join(path1,"kur_sab%d"%seeds)
            os.makedirs(path4,exist_ok=True)
            mat=mater()
            mat.add_s_alpha_beta("c_H_in_H2O")
            stat_path=run(disp,tal2,path4,batch=100,n_bat=10000,seed=seeds,mat=mat)
            dirc3={}
            for a,b,c in KUR_TALLIES:
                m,s=tallies(stat_path,a)
                dirc3[a]={"mean":m,"se":s}
            dirc3["batches"]=100
            dirc3["particles"]=10000
            path7=os.path.join(OUT_DIR,"openmc_kur_sab%d.json"%seeds)
            with open(path7,"w") as f:
                json.dump(dirc3,f)
            
    if run_Tru2==True:
        disp=kur()
        tal2=tru_tally()
        path4=os.path.join(path1,"kur_s2")
        os.makedirs(path4,exist_ok=True)
        mat=mater()
        stat_path=run(disp,tal2,path4,batch=100,n_bat=10000,seed=2,mat=mat)
        dirc3={}
        for a,b,c in KUR_TALLIES:
            m,s=tallies(stat_path,a)
            dirc3[a]={"mean":m,"se":s}
        dirc3["batches"]=100
        dirc3["particles"]=10000
        path7=os.path.join(OUT_DIR,"openmc_kur_s2.json")
        with open(path7,"w") as f:
            json.dump(dirc3,f)
                
    if run_as1==True:
        tal1=tally()
        E_dis=openmc.stats.Discrete([37.2],[1.0])
        path2=os.path.join(path1,"j1")
        os.makedirs(path2,exist_ok=True)
        mat=mater()
        path3=run(E_dis,tal1,path2,batch=100,n_bat=1000,mat=mat)
        m,s=tallies(path3,'j1 uncollided flux')
        dirc={
            "flux0_mean":m,
            "flux0_se":s,
            "batches":100,
            "particles":1000
        }
        path4=os.path.join(OUT_DIR,"openmc_j1.json")
        with open(path4,"w") as f:
                json.dump(dirc,f)
        



    