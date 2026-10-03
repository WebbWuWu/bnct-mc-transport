import math
from rng import lcg, isotropic_direction
import os
import h5py
import numpy as np
from xs_lookup import LIB,sigma_at

pi=math.pi
def path_A_prob(y):
    return y*math.sqrt(pi)/(y*math.sqrt(pi)+2)

def sample_target(y,rng):
    for k in range(1000):
        xi2=rng.random()
        xi3=rng.random()
        xi4=rng.random()
        if rng.random()<path_A_prob(y):
            x2=-math.log(1-xi2)-math.log(1-xi3)*math.cos((1-xi4)/2*pi)**2
        else:
            x2=-math.log((1-xi2)*(1-xi3))
        x=math.sqrt(x2)
        mu_T=rng.random()*2-1
        arg=x**2+y**2-2*x*y*mu_T
        V_rel=math.sqrt(max(arg,0.0))
        pac=V_rel/(x+y)
        if rng.random()<pac:
            return x,mu_T,k+1
    else:
        assert False
