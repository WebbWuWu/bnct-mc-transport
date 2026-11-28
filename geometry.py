import math

EPS=1e-9#单位 cm，「在面上」的门槛，离面距离小于它就按方向判 sense


class ZPlane:
    def __init__(self,z0):
        self.z0=z0
    def f(self,x,y,z):
        return z-self.z0
    def normal(self,x,y,z):
        return 0.0,0.0,1.0
    def distance(self,x,y,z,u,v,w,on_surf):
        if on_surf:
            return math.inf
        if w==0.0:
            return math.inf
        t=(self.z0-z)/w
        if t<=0:
            return math.inf
        else:
            return t

class Sphere:
    def __init__(self,x0,y0,z0,R):
        self.x0=x0
        self.y0=y0
        self.z0=z0
        self.R=R
    def f(self,x,y,z):
        dx=x-self.x0
        dy=y-self.y0
        dz=z-self.z0
        return (dx**2+dy**2+dz**2)-self.R**2
    def normal(self,x,y,z):
        dx=x-self.x0
        dy=y-self.y0
        dz=z-self.z0
        return 2*dx,2*dy,2*dz
    def distance(self,x,y,z,u,v,w,on_surf):
        dx=x-self.x0
        dy=y-self.y0
        dz=z-self.z0
        k=dx*u+dy*v+dz*w
        if on_surf:
            if k<0:
                return -2*k
            else:
                return math.inf
        fv=dx**2+dy**2+dz**2-self.R**2
        dirc=k**2-fv
        if dirc<0:
            return math.inf
        sq=math.sqrt(dirc)
        if fv<0:
            return -k+sq
        else:
            if k>=0:
                return math.inf
            else:
                return -k-sq
                
class ZCylinder:
    def __init__(self,x0,y0,R):
        self.x0=x0
        self.y0=y0
        self.R=R
    def f(self,x,y,z):
        dx=x-self.x0
        dy=y-self.y0
        return dx**2+dy**2-self.R**2
    def normal(self,x,y,z):
        dx=x-self.x0
        dy=y-self.y0
        return 2*dx,2*dy,0.0
    def distance(self, x, y, z, u, v, w, on_surf):
        dx=x-self.x0
        dy=y-self.y0
        a=u*u+v*v
        if a==0.0:
            return math.inf
        k=dx*u+dy*v
        if on_surf:
            if k<0:
                return -2*k/a
            else:
                return math.inf
        fv=dx**2+dy**2-self.R**2
        disc=k**2-a*fv
        if disc<0:
            return math.inf
        sq=math.sqrt(disc)
        if fv<0:
            return (-k+sq)/a
        else:
            if k>=0:
                return math.inf
            else:
                return (-k-sq)/a

def sense_of(surf,x,y,z,u,v,w):
    fv=surf.f(x,y,z)
    nx,ny,nz=surf.normal(x,y,z)
    nlen=math.sqrt(nx**2+ny**2+nz**2)
    if abs(fv)<EPS*nlen:
        dot=u*nx+v*ny+w*nz
        if dot>=0:
            return 1
        else:
            return -1
    else:
        if fv>=0:
            return 1
        else:
            return -1
