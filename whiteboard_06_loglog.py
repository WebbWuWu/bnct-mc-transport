# -*- coding: utf-8 -*-
# ============================================================
#  whiteboard_06_loglog.py  --  白板题⑥ (3)：log-log 插值的两种形式与互证
#  2026-09-09
# ============================================================
#  【出处声明】下面 exp 形式那一行的两个括号，是用户在草稿纸上推出来、
#  口述给 Claude 的：分子 ln(S[k+1]/S[k])，分母 ln(E[k+1]/E[k])。
#  Claude 只做了打字和验证，没有参与推导。
#
#  为什么单独一个文件：whiteboard_06_interp.py 里 s()/s1() 的签名是
#  (e, S) 两个参数、调用处只给一个，跑到那儿就 TypeError，
#  追加在后面的代码永远执行不到。等那两处修好再合并。
# ============================================================
import math


def loglog_exp(e, E, S, k):
    """exp 形式：把 (1) 的线性插值公式整个搬进 log 空间，算完再 exp 回来。"""
    return math.exp(math.log(S[k])
                    + (math.log(e) - math.log(E[k]))
                    * math.log(S[k+1] / S[k])
                    / math.log(E[k+1] / E[k]))


def loglog_power(e, E, S, k):
    """幂形式：上面那个用指数律化简的结果。p 是双对数斜率。"""
    p = math.log(S[k+1] / S[k]) / math.log(E[k+1] / E[k])
    return S[k] * (e / E[k]) ** p


# ============================================================
#  自测
# ============================================================
if __name__ == "__main__":
    # 真实数据：B-10 的 (n,alpha)，ENDF/B-VIII.0（lib80x_hdf5），294K
    E = [2.53e-2, 1.0]          # eV
    S = [3844.0, 611.3]         # barn
    k = 0

    # ---- 判据 1：两种形式必须逐点一致 ----
    #   这条抓「化简推错了」。两条独立的算路走到同一个数，才说明化简是对的。
    #   容差用相对的：截面跨十个数量级，绝对容差不可能同时适合两头。
    print("   e            exp形式          幂形式         相对差")
    worst = 0.0
    for e in [E[0], 0.05, 0.1, 0.3, E[1]]:
        a = loglog_exp(e, E, S, k)
        b = loglog_power(e, E, S, k)
        rel = abs(a / b - 1.0)
        worst = max(worst, rel)
        print("  %-8.4g   %14.8f  %14.8f   %.2e" % (e, a, b, rel))
    assert worst < 1e-14, "两种形式不一致，最大相对差 %.3e" % worst
    print("判据1：两种形式逐点一致，最大相对差 %.2e" % worst)

    # ---- 判据 2：三条自检（左端点 / 右端点 / 几何中点）----
    #   log 空间里的「中点」是几何中点，「平均」是几何平均 —— 中点判据要跟着换空间。
    #   全部用相对容差。
    d_lo  = abs(loglog_exp(E[k],   E, S, k) / S[k]   - 1.0)
    d_hi  = abs(loglog_exp(E[k+1], E, S, k) / S[k+1] - 1.0)
    g_e   = math.exp((math.log(E[k]) + math.log(E[k+1])) / 2)   # 几何中点
    g_s   = math.exp((math.log(S[k]) + math.log(S[k+1])) / 2)   # 几何平均
    d_mid = abs(loglog_exp(g_e, E, S, k) / g_s - 1.0)
    for name, d in (("左端点", d_lo), ("右端点", d_hi), ("几何中点", d_mid)):
        assert d < 1e-12, "%s自检不合格，相对差 %.3e" % (name, d)
    print("判据2：左端点 %.1e，右端点 %.1e，几何中点 %.1e" % (d_lo, d_hi, d_mid))

    # ---- 判据 3：1/v 律反查 ----
    #   B-10 的 (n,alpha) 在热能区服从 1/v，即 sigma ∝ E^(-1/2)。
    #   所以从真实数据算出来的 p 必须是 -1/2 —— 这是公式自带的判据，
    #   不是我挑的容差：p 偏离 -0.5 太远就说明公式推错了。
    p = math.log(S[k+1] / S[k]) / math.log(E[k+1] / E[k])
    print("判据3：实测 p = %.6f，1/v 律要求 -0.5，偏差 %.2e" % (p, abs(p + 0.5)))
    assert abs(p + 0.5) < 1e-3, "p = %.6f，不是 -1/2，公式推错了" % p

    print("三条判据全过。")
