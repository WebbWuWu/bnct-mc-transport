from lesson08_dose import OUT_DIR,MASS_FRAC,RHO,T,D,kur_table
import os
import openmc
import math

openmc.config['cross_sections']="/mnt/e/BNCT_GPU/lib80x_hdf5/cross_sections.xml"
path1=os.path.join(OUT_DIR,'openmc')
os.makedirs(path1,exist_ok=True)

soft=openmc.Material(name='soft tissue')
MASS_sum=0.0
for nuc in MASS_FRAC.keys():
    soft.add_nuclide(nuc,MASS_FRAC[nuc],'wo')
    MASS_sum+=MASS_FRAC[nuc]
soft.set_density('g/cm3',RHO*MASS_sum)
soft.temperature=294#和 lesson08 的 T 一致
mats=openmc.Materials([soft])
path2=os.path.join(path1,'materials.xml')
mats.export_to_xml(path2)
    
xp1=openmc.XPlane(x0=7.3,boundary_type='reflective')
xp2=openmc.XPlane(x0=-7.3,boundary_type='reflective')
yp1=openmc.YPlane(y0=7.3,boundary_type='reflective')
yp2=openmc.YPlane(y0=-7.3,boundary_type='reflective')
zp1=openmc.ZPlane(z0=D,boundary_type='vacuum')
zp2=openmc.ZPlane(z0=0,boundary_type='vacuum')

bou=-xp1&+xp2&-yp1&+yp2&-zp1&+zp2

c=openmc.Cell(fill=soft,region=bou)
geo=openmc.Geometry([c])
path3=os.path.join(path1,'geometry.xml')
geo.export_to_xml(path3)

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

pos=openmc.stats.Point((0,0,1e-6))
dirn=openmc.stats.Monodirectional((0,0,1))
S=openmc.IndependentSource(space=pos,angle=dirn,energy=dsp)

st=openmc.Settings()
st.run_mode='fixed source'
st.particles=1000
st.batches=100
st.source=S
path4=os.path.join(path1,'settings.xml')
st.export_to_xml(path4)

