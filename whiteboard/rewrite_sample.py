import math
from rng import lcg, isotropic_direction
import os
import h5py
import numpy as np
from xs_lookup import LIB,sigma_at,find_cell,load_nuclide

def path_A_prob(y):
    return y*math.sqrt(math.pi)/(y*math.sqrt(math.pi)+2)

def sample_target(y,rng):
    for k in range(1000):
        xi2=rng.random()
        xi3=rng.random()
        xi4=rng.random()
        if rng.random()<path_A_prob(y):
            x2=-math.log(1-xi2)-math.log(1-xi3)*math.cos((1-xi4)/2*math.pi)**2
        else:
            x2=-math.log((1-xi2)*(1-xi3))
        x=math.sqrt(x2)
        mu_T=2*rng.random()-1
        arg=x2+y**2-2*mu_T*x*y
        v_rel=math.sqrt(max(arg,0.0))
        pac=v_rel/(x+y)
        if rng.random()<pac:
            return x,mu_T,k+1
    else:
        assert False
