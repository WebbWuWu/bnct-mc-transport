import math
from whiteboard_05_rewrite import find

def s(e,S):
    k=find(E,e)
    if E[k]==0:
        E[k]=1e-12
    return (S[k]+(e-E[k])*(S[k+1]-S[k])/(E[k+1]-E[k]))
if __name__=="__main__":
    S=[1,2,3,4,5]
    E = [1.0, 2.0, 3.0, 4.0, 5.0]
    e=1
    k=find(E,e)
    e1=E[k]
    assert s(e1)==S[k],"端点自检不合格，S（e）=%.6f"%s(e1)
    e2=E[k+1]
    assert s(e2)==S[k+1],"端点自检不合格，S（e）=%.6f"%s(e2)
    e3=(E[k]+E[k+1])/2
    assert s(e3)==(S[k]+S[k+1])/2,"中点自检不合格，S（e）=%.6f"%s(e3)
    
def s1(e,S):
    k=find(E,e)
    return S[k] * (e/E[k])**(math.log(S[k+1]/S[k])/math.log(E[k+1]/E[k]))
    
if __name__=="__main__":
    k=find(E,e)
    e1=E[k]
    assert s1(e1)==S[k],"端点自检不合格，S（e）=%.6f"%s1(e1)
    e2=E[k+1]
    assert s1(e2)==S[k+1],"端点自检不合格，S（e）=%.6f"%s1(e2)
    e3=math.exp((math.log(E[k]*E[k+1])/2))
    assert abs(s1(e3)-math.exp((math.log(S[k]*S[k+1])/2)))<1e-12,"中点自检不合格，S（e）=%.6f"%s1(e3)