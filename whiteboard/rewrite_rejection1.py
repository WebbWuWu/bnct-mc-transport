# ============================================================
#  rewrite_rejection1.py —— 默写 sample_rejection（第 1 次）
#  2026-10-11 Claude 开文件，只写了这段注释；下面全部由你写
# ============================================================
#  要写：
#    1. import（要哪些自己想）
#    2. sample_rejection（白板题⑧ 舍选法）
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：whiteboard_08_rejection.py、whiteboard_08_checks.py、白板⑧规格书、讲义
#  写完：存盘，把文件名发给 Claude 对拍（同种子跑原函数和你的函数，逐位比）
# ============================================================
from rng import lcg
Max=1000
def sample_rejection(f,g,g_sample,m,rng):
    for k in range(Max):
        x=g_sample(rng)
        y=m*g(x)*rng.random()
        if y<f(x):
            return x,k+1
    else:
        return None
