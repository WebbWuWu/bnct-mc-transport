from rng import lcg
Max=1000
def sample_rejection(f,g,g_sample,m,rng):
    for k in range(Max):
        x=g_sample(rng)
        if f(x)>g(x)*m*rng.random():
            return x,k+1
    else:
        return None
