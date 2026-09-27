import math
from rng import lcg

Max = 1000                                    # 第 2 问：p≥0.1 时撞满概率 ≤ 1.8e-46

def sample_rejection(f, g, g_sample, M, rng):
    # 入口防护（Claude 给）
    assert M > 0, "M 必须为正：M=%g" % M
    for k in range(Max):
        x = g_sample(rng)                    # 第 1 步：按 g 抽 x
        # 检查 M 够不够大（Claude 给，第 5 问）
        assert f(x) <= M * g(x), \
            "M 给小了：x=%g 处 f(x)=%g > M*g(x)=%g" % (x, f(x), M * g(x))
        y = rng.random() * M * g(x)          # 第 2 步：y 均匀在 [0, M·g(x))
        if y < f(x):                         # 第 3 步：严格 <（第 3 问）
            return x, k + 1                  # 样本 + 用了几轮（第 4 问）
    else:
        return None                          # 撞满上限（第 2 问：乙）

if __name__ == "__main__":
    # 演示：f(x)=2x, g(x)=1, g_sample=均匀, M=2, 种子 2026
    def f(x):
        return 2*x
    def g(x):
        return 1
    rng=lcg(2026)
    def g_sample(rng):
        return rng.random()
    x,n=sample_rejection(f,g,g_sample,2,rng)
    print("x=%g,n=%g"%(x,n))
    N=100000
    xx=[0.0]*N
    nn=[0.0]*N
    for i in range(N):
        xx[i],nn[i]=sample_rejection(f,g,g_sample,2,rng)
    xx_p=sum(xx)/N
    nn_p=sum(nn)/N
    print("x的平均值=%g,理论值=0.6667, n的平均值=%g,理论值=2 "%(xx_p,nn_p))
    
    # 判据脚本：你写完函数后 Claude 给
