import mcstat
import math
import os
from xs_lookup import sigma_at, BARN
from sample_reaction import sample_reaction
import time

HERE = os.path.dirname(os.path.abspath(__file__))

from rng import lcg, isotropic_direction

Max=10000

NL=20

def run_slab(N,seed,d,nuclides,E_src,order,z_src=None,E_cut=0.0253):
    if z_src is None:
        z_src=d/2
    assert E_src > E_cut, "E_src=%.6g 不高于 E_cut=%.6g，每条历史会零碰撞" % (E_src, E_cut)

    assert len(order) == len(nuclides)

    for i in range(len(nuclides)):
        assert len(order[i]) == len(nuclides[i]["XS"])
    len_nuc=len(nuclides)
    rng=lcg(seed)
    n_T=0
    n_R=0
    n_A=0
    n_stuck=0
    n_cut=0
    n_i=0
    dz=d/NL
    track_sum=[0.0]*NL
    track_sq=[0.0]*NL
    v_sum=v_sq=0.0
    tol_sigt_all=0.0
    for i in range(N):
        z=z_src
        E=E_src
        tol_sigt_s=0.0
        w=1.0
        track_h=[0.0]*NL
        n_hist=0
        for step in range(Max):
            sig_ds=[]
            sig_al=[]
            tot_ls=[]
            for j in range(len(nuclides)):
                Sig_j,tot_j=sigma_at(nuclides[j]["E"],nuclides[j]["XS"],E)
                Sigma_j = nuclides[j]["dens"] * tot_j * BARN
                sig_ds.append(Sig_j)
                sig_al.append(Sigma_j)
                tot_ls.append(tot_j) 
            Sigma_t=sum(sig_al)
            assert Sigma_t > 0.0, "Sigma_t 为零，飞行距离会除零：E=%.6g eV" % E
            # 为什么是 1-rng.random() 而不是 rng.random()：
            #   random() 返回 [0,1) —— 0 取得到，1 取不到。若写 -log(xi)，xi=0 时 log 收到 0 → -inf，程序炸。
            #   取 1-xi 后 log 的参数落在 (0,1]，永远取不到 0。两种写法同分布，但安全性完全不同。
            #   命中概率约 2^-32：跑一亿条历史大概撞一次 —— 调试时永远不出现，演示时出现。
            #   这个选择配套的是本文件里 LCG 的区间，换 RNG 要重新核对。见 whiteboard_02_spec.md。
            s=-math.log(1-rng.random())/Sigma_t
            z0=z
            z1=z0+s*w
            if z1>d:
                z_end=d
            elif z1<0:
                z_end=0.0
            else:
                z_end=z1
            assert 0.0 <= z0 <= d, "z0 跑出板外了: %.15f" % z0
            seg = 0.0                                  # 本段实际记了多少
            if w == 0.0:
                k0 = int(z0 / dz)
                assert 0 <= k0 < NL, "层号越界: k0=%d z0=%.9f" % (k0, z0)
                track_h[k0] += s
                seg    += s
                expect  = s

            else:
                k0 = int(z0 / dz)                      # 起点层
                k1 = int(z_end / dz)                   # 终点层
                if k1 >= NL:                           # z_end 恰好 == d 时夹一下
                    k1 = NL - 1

                cur = z0                               # 笔尖当前位置
                k   = k0                               # 笔尖当前层
                guard=0
                while True:
                    guard+=1
                    assert guard<=NL+2,\
                        "拆层循环未退出"
                    # ① 这一步走到哪：当前层的出口边界 vs z_end，取近的那个
                    #    w>0 出口是 (k+1)*dz ；w<0 出口是 k*dz
                    if w>0:
                        k_o=(k+1)*dz
                        z_step=min(k_o,z_end)
                    else:
                        k_o=k*dz
                        z_step=max(k_o,z_end)
                    # ② 记账：dz_step = abs(这一步终点 - cur)
                    #         track_h[k] += dz_step / abs(w)
                    #         seg        += 同一个值
                    #         记账前先 assert 0 <= k < NL
                    assert 0<=k<NL,\
                        "k越界"
                    dz_step=abs(z_step-cur)
                    track_h[k] += dz_step / abs(w)
                    seg+=dz_step/abs(w)
                    # ③ cur = 这一步的终点
                    cur=z_step
                    # ④ 到 k1 了就 break，否则 k += 1（w>0）或 k -= 1（w<0）
                    if k==k1:
                        break
                    if w > 0:
                        k += 1
                    else:
                        k -= 1

                expect = abs(z_end - z0) / abs(w)

            assert math.isclose(seg, expect, rel_tol=1e-9, abs_tol=1e-12), \
                "段记账不守恒: seg=%.15f expect=%.15f z0=%.6f z_end=%.6f w=%.9f" \
                % (seg, expect, z0, z_end, w)
            tol_sigt_s += Sigma_t * seg
            z = z1                                     
            if z<0:
                n_R+=1
                break
            if z>d:
                n_T+=1
                break
            if len_nuc==1:
                jn=0
            else:
                raise NotImplementedError
            xi=rng.random()
            k=sample_reaction(sig_ds[jn],tot_ls[jn],xi,order[jn])
            if nuclides[jn]["mts"][k]==2:
                u,v,w=isotropic_direction(rng)
                n_i+=1
                n_hist+=1
                if nuclides[jn]["awr"]is not None:
                    A=nuclides[jn]["awr"]
                    al=((A-1)/(A+1))**2
                    xii=rng.random()
                    E=E*(al+(1-al)*xii)
                if E<=E_cut:
                    n_cut+=1
                    break
            else:
                n_A+=1
                n_i+=1
                n_hist+=1
                break
        else:
            n_stuck+=1
        for j in range(NL):
            track_sum[j]+=track_h[j]
            track_sq[j]+=track_h[j]**2
        x = tol_sigt_s - n_hist
        tol_sigt_all += tol_sigt_s
        v_sum+=x
        v_sq+=x**2
    if n_stuck>0:
        print("warning!有%d条历史触顶max_event"% n_stuck)
    assert n_T+n_R+n_A+n_stuck+n_cut == N
    scale = abs(tol_sigt_all) + abs(n_i)
    assert math.isclose(tol_sigt_all - n_i, v_sum,rel_tol=0.0, abs_tol=1e-9*scale), \
                        "验证B总账不一致：累加器=%.6f  n_i=%d  v_sum=%.6f" % (tol_sigt_all, n_i, v_sum)
    return n_T/N,n_R/N,n_A/N,n_i/N,track_sum,track_sq,dz,v_sum,v_sq,n_cut,n_stuck

if __name__ == "__main__":
    from xs_lookup import load_nuclide
    E,mts,XS,awr=load_nuclide("H1","294K")

    nuc_H1={
        "E":E,"mts":mts,"XS":XS,"awr":awr,
        "dens":6.4606e22          # 9/15 实测：软组织里 H-1 的原子数密度 cm^-3
    }

    # order 在【代表能量 0.0253 eV】上排一次，之后整场不变。
    # 已知近似清单 第8条：排序顺序理论上随能量变，这里固定成热能点。
    sig,tot=sigma_at(E,XS,0.0253)
    order_H1=sorted(range(len(sig)),key=lambda t: sig[t],reverse=True)
    order=[order_H1]
    assert mts[order_H1[0]]==2, "排第一的应该是弹性(MT=2)，实际 MT=%d"%mts[order_H1[0]]

    # ---- 判据 B：与 9/15 的 xslib.total_xs 在 0.0246 eV 对拍（两条独立的路）----
    #   9/15 走 xslib + log-log 插值：sigma=30.6562  Sigma=1.980574
    #   今天走 xs_lookup + lin-lin  ：见下
    #   两者差约 2.3e-4，落在 9/17 实测的 log-log vs lin-lin 差异区间内（中位 6.3e-5，最大 9.6e-4）
    #   -> 这是【已解释的差】，不是 bug。
    _sig246,_tot246 = sigma_at(E,XS,0.0246)
    print("判据B  0.0246 eV : sigma=%.4f barn  Sigma=%.6f cm^-1   (9/15: 30.6562 / 1.980574)"
          % (_tot246, nuc_H1["dens"]*_tot246*BARN))
    print("       0.0253 eV : sigma=%.4f barn  Sigma=%.6f cm^-1   (9/18: 30.4138 / 1.964914)"
          % (tot, nuc_H1["dens"]*tot*BARN))
    print()

    nuc_fake={
        # 单能极限用的假库。规格见 2026-09-19 日志。
        #   两个能量点、两点截面相同 -> lin-lin 给常数 -> Sigma_t 恒定（前提2）
        #   道0 必须是散射，才和 lesson03 的 "if xi < sigma_s/sigma_t" 同一个条件
        #   awr=None 是甲方案哨兵 -> 不更新能量、不抽那个随机数（前提3）
        #   dens 取 1e24，使 dens*tot*BARN 正好等于 lesson03 的 sigma_t=1.0
        "E":[1e-5, 2e7], "mts":[2, 102], "XS":[[0.8,0.8],[0.2,0.2]],
        "awr":None, "dens":1e24
    }

    N=100000
    print("%-26s %9s %11s %11s %13s %8s"
          % ("步骤5 计时","耗时(s)","hist/s","碰撞/历史","coll/s","相对A"))

    import lesson03_tally      # 有 __main__ 保护，import 不会跑它的三个算例
    rows=[]
    def timeit(tag, n, fn):
        # 🔴 n 必须显式传：C 用的历史数和 A/B 不一样，
        #    用全局 N 算速度会把 C 的数字放大 10 倍（第一版就是这么错的）。
        t0=time.perf_counter(); r=fn(); dt=time.perf_counter()-t0
        coll_his=r[3]                         # 第 4 个返回值就是 n_i/N，三个算例通用
        rows.append((tag, dt, n/dt, coll_his, n*coll_his/dt))

    # A 基准：Sigma_t 是常数，取它的代价 = 0
    timeit("A lesson03 基准", N,
           lambda: lesson03_tally.run_slab(N, 2026, 2.0, 1.0, 0.8))

    # B 单能极限：物理和 A 逐位相同（判据1），唯一差别是走了查表那套机器
    #   -> A→B 的比值 = 查表机器本身的固定开销
    timeit("B lesson06 单能极限", N,
           lambda: run_slab(N, 2026, 2.0, [nuc_fake], 1.0, [[0,1]],
                            z_src=0.0, E_cut=1e-30))

    # C 真库慢化算例：d=1000 且 z_src=d/2，漏出 0.000%（9/18 实测）
    #   🔴 不能用 d=2 —— 1 MeV 下 H-1 自由程 3.63 cm > 板厚，中子一步飞出，
    #      碰撞/历史 只有 1.48，和 A/B 的 2.25 不是一回事，比值没有意义。
    timeit("C lesson06 H-1 慢化", N//10,
           lambda: run_slab(N//10, 2026, 1000.0, [nuc_H1], 1.0e6, order,
                            z_src=500.0, E_cut=0.0253))

    base=rows[0][4]
    for tag,dt,hs,ch,cs in rows:
        print("%-26s %9.3f %11.1f %11.4f %13.1f %7.3fx"%(tag,dt,hs,ch,cs,cs/base))

    print()
    print("判据A  C 那一行的 碰撞/历史 应落在 18.5 附近（1 MeV -> 0.0253 eV，含过冲）")
    print("       跑出两万 -> alpha=((A-1)/(A+1))**2 那个括号写错了")
    print("🔴 别用调试器(F5)跑这一段。对照：9/18 实测 lesson03 = 103802.6 hist/s，")
    print("   A 那一行差一个量级就说明环境被污染了，整张表不能看。")
