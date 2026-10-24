# ============================================================
#  rewrite_sample_reaction2.py —— 默写 sample_reaction（第 2 次，10/22）
#  2026-10-21 Claude 开文件，只写了这段注释；下面全部由你写
# ============================================================
#  要写：
#    1. import（要哪些自己想；函数本身用不到的就不写）
#    2. sample_reaction（按截面抽反应道，带 order 排序 + 扫完没返回时的兜底）
#  不写：assert、if __name__ == "__main__": 那一行和下面整块、print
#  不看：sample_reaction.py、whiteboard\copy_sample_reaction.py、讲解_sample_reaction.md、rewrite_sample_reaction1.py
#  水课在纸上手写：晚上照纸敲进来就行（照纸敲不算看原文）
#  写完：存盘，把文件名发给 Claude 对拍（同一组输入跑原函数和你的函数，逐位比）
# ============================================================
def sample_reaction(sigma,total,xi,order):
    n=len(sigma)
    acc=0.0
    for i in range(n):
        rank=order[i]
        acc+=sigma[rank]
        if xi<acc/total:
            return rank
    for i in range(n-1,-1,-1):
        k=order[i]
        if sigma[k]>0:
            return k
