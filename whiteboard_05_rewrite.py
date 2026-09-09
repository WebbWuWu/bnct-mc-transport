import math
def find(E,e):
    n=len(E)
    if e>=E[n-1]:
            return n-2
    if e<=E[0]:
            return 0
    lo=0
    hi=n-1
    step=0
    while True:
        step+=1
        assert step < 100, "二分没有收敛：lo=%d hi=%d（区间没在变小）" % (lo, hi)
        mid=(hi-lo)//2+lo
        if e>E[mid]:
            lo=mid
        if e<E[mid]:
            hi=mid
        if e==E[mid]:
            return mid
        if hi-lo==1:
            return lo
if __name__ == "__main__":
    import bisect, random

    E = [1.0, 2.0, 3.0, 4.0, 5.0]

    # 验证 1：三个边界 —— 抓「循环之前没把不变量建立起来」
    #   这三个数你要先在纸上写出期望值，再跑。不许跑完了再倒推期望。
    expect = {0.5:0, 1.0:0, 3.0:2, 4.999:3, 5.0:3, 7.0:3}
    for e in [0.5, 1.0, 3.0, 4.999, 5.0, 7.0]:
        k = find(E, e)
        print("e=%-7s -> k=%s" % (e, k))
        assert k == expect[e], "e=%s 期望 %s，实得 %s" % (e, expect[e], k)
    # 验证 2：与 bisect 对拍 —— 抓「比较写反了 / 差一格」
    #   两个独立算出来的量：你的 find 是一个，标准库是另一个。
    random.seed(2026)
    n_bad = 0
    for _ in range(2000):
        e = random.uniform(E[0], E[-1] - 1e-12)
        k_ref = bisect.bisect_right(E, e) - 1
        if find(E, e) != k_ref:
            n_bad += 1
    print("验证 2: 2000 次随机对拍，不一致 %d 次" % n_bad)
    assert n_bad == 0, "与 bisect 不一致 %d 次" % n_bad

    # 验证 3：比较次数 —— 抓「循环没有真的在砍一半」（退化成线性也能通过验证 1、2！）
    #   这是三条里唯一能抓住「你写的其实是遍历」的一条。
    N = 100000
    Ebig = [float(i) for i in range(N)]
    calls = [0]

    class Counted(list):                 # 每读一次 E[i] 计一次数
        def __getitem__(self, i):
            calls[0] += 1
            return list.__getitem__(self, i)

    Ec = Counted(Ebig)
    random.seed(7)
    worst = 0
    for _ in range(200):
        calls[0] = 0
        find(Ec, random.uniform(0.0, float(N - 1)))
        worst = max(worst, calls[0])
    # 上界放宽到 4*log2(n)：这条判据抓的是**数量级**，不是常数。
    # 二分无论怎么写都在 ~20 次量级；遍历是 50000 次。中间隔着三个数量级，不会误判。
    bound = 4 * int(math.log2(N))
    print("验证 3: n=%d，最坏读表次数 = %d（上界 %d，线性扫描要 ~%d）"
          % (N, worst, bound, N // 2))
    assert worst <= bound, "读表 %d 次，超过 log2(n)：这不是二分，是遍历" % worst
    
    assert find(E, E[-1]+1) < (len(E)-1),\
        "E(n)不存在"

    print("三条全过。")
