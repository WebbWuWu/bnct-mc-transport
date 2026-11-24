import os, sys,openmc,json
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)
from openmc_tissue import geom,mesh,tallies,path1
from openmc_water import mater_water
from dose import OUT_DIR
from photon import E_G_H,E_CUT_G

E0=E_G_H
BATCH=100
N_BAT=100000
SEED=4471

PHOTON_TALLIES=[
    ("photon heating","heating"),
    ("photon total heating","heating")
]

def photon_tally():
    ts=openmc.Tallies()
    for a,c in PHOTON_TALLIES:
        t=openmc.Tally(name=a)
        t.scores=[c]
        if a=="photon heating":
            t.filters=[mesh()]
        ts.append(t)
    return ts

def run_photon(tal,inp,batch,n_bat,mat,seed,electron):
    os.makedirs(inp,exist_ok=True)
    material=openmc.Materials([mat])
    geo=geom(mat)
    pos=openmc.stats.Point((0,0,1e-6))
    dirn=openmc.stats.Monodirectional((0,0,1))
    S=openmc.IndependentSource(space=pos,angle=dirn,energy=openmc.stats.Discrete([E0], [1.0]),particle='photon')

    st=openmc.Settings()
    st.run_mode='fixed source'
    st.particles=n_bat
    st.batches=batch
    st.seed=seed
    st.source=S
    st.photon_transport=True
    st.electron_treatment=electron
    st.cutoff={'energy_photon':E_CUT_G}
    
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

if __name__ == "__main__":
    mat=mater_water()
    for e in ("led","ttb"):
        tal=photon_tally()
        inp=os.path.join(path1,"photon_"+e)
        sp=run_photon(tal,inp,BATCH,N_BAT,mat,SEED,e)
        res={}
        for a,c in PHOTON_TALLIES:
            mean,se=tallies(sp,a)
            res[a]={"mean":mean,"se":se}
        res["batches"]=BATCH
        res["particles"]=N_BAT
        res["seed"]=SEED
        res["electron"]=e
        res["E0"]=E0
        path=os.path.join(OUT_DIR,"openmc_photon_"+e+".json")
        with open(path,"w")as f:
            json.dump(res,f)
            print(path)
            