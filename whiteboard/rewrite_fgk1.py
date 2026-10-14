import math
from rng import lcg, isotropic_direction
def free_gas_kinematics(E,u,v,w,A,kT,x,mu_T,rng):
    y=math.sqrt(E*A/kT)
    vnx=y*u
    vny=y*v
    vnz=y*w
    for k in range(100):
        uu,vv,ww=isotropic_direction(rng)
        al=uu*u+vv*v+ww*w
        rx=uu-u*al
        ry=vv-v*al
        rz=ww-w*al
        if math.sqrt(rx**2+ry**2+rz**2)>=1e-6:
            break
    else:
        assert False
    e1x=rx/math.sqrt(rx**2+ry**2+rz**2)
    e1y=ry/math.sqrt(rx**2+ry**2+rz**2)
    e1z=rz/math.sqrt(rx**2+ry**2+rz**2)
    
    rbx=mu_T*u+math.sqrt(1-mu_T**2)*e1x
    rby=mu_T*v+math.sqrt(1-mu_T**2)*e1y
    rbz=mu_T*w+math.sqrt(1-mu_T**2)*e1z
    
    vbx=x*rbx
    vby=x*rby
    vbz=x*rbz
    
    vcmx=(vnx+A*vbx)/(1+A)
    vcmy=(vny+A*vby)/(1+A)
    vcmz=(vnz+A*vbz)/(1+A)
    
    vncmx=vnx-vcmx
    vncmy=vny-vcmy
    vncmz=vnz-vcmz
    vncm=math.sqrt(vncmx**2+vncmy**2+vncmz**2)
    
    u1,v1,w1=isotropic_direction(rng)
    vlx=u1*vncm+vcmx
    vly=v1*vncm+vcmy
    vlz=w1*vncm+vcmz
    vl=math.sqrt(vlx**2+vly**2+vlz**2)
    
    El=(kT/A)*(vl**2)
    u2=vlx/vl
    v2=vly/vl
    w2=vlz/vl
    
    return El,u2,v2,w2
            
            