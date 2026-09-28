import math
from rng import lcg
def isotropic_direction(rng):
    mu=2*rng.random()-1
    phi=2*rng.random()*math.pi
    sin_theta=math.sqrt(1-mu**2)
    u=sin_theta*math.cos(phi)
    v=sin_theta*math.sin(phi)
    w=mu
    return u,v,w

def mean_se(s, sq, n):
    mean=s/n
    var=n/(n-1)*(sq/n-mean**2)
    se=math.sqrt(var/n)
    return mean,var,se
