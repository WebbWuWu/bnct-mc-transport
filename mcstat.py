import math

def mean_se(s,sq,n):
    """由 Σx、Σx²、n 算样本均值、贝塞尔方差、均值标准误。"""
    assert n>1,\
        "n取值不对，n=%d"%n
    mean=s/n
    var=n/(n-1)*(sq/n-mean**2)
    assert var>=0,\
        "var范围出错，var=%.6f"%var
    se=math.sqrt(var/n)
    return mean,var,se
        
def dev_sigma(meas,ref,se):
    """标准化"""
    return (meas-ref)/se

def row(name, meas, ref, se):
    """打印表中每一行"""
    print("  %-12s %14.6f %16.6f %11.2e %+10.2f" % (name, meas, ref, se, dev_sigma(meas,ref,se)))

def header():
    """打印表头"""
    print("  %-12s %14s %16s %11s %10s"%("check","measured","analytic","SE","dev(sigma)"))

def chi2_uniform(counts,expected):
    """计算卡方"""
    sum_chi2=0
    assert expected >= 5, "期望频数 %.2f < 5，卡方近似不成立" % expected
    for k in range(len(counts)):
        sum_chi2+=(counts[k]-expected)**2/expected
    return sum_chi2

def check(name,meas,ref,se,k=3.0):
    """检查每一行值是否正确"""
    row(name, meas, ref, se)
    assert abs(meas-ref)<k*se,\
        "name=%-12s检验不合格,meas=%14.6f,ref=%16.6f,se=%11.2e,偏离=%.6f" % (name, meas, ref, se,dev_sigma(meas,ref,se))


if __name__ == "__main__":
    x = [1, 2, 3, 4]
    s  = 10.       # 你手算的 Σx
    sq = 30       # 你手算的 Σx²
    mean, var, se = mean_se(s, sq, 4)
    assert abs(var - 5/3) < 1e-12, "自测1方差不合格，var=%.6f"%var   # ... 是你手算的方差

    counts = [12, 8, 10, 10]
    expected=10
    chi2_hand=0.8
    chi2_got=chi2_uniform(counts,expected)
    assert abs(chi2_hand-chi2_got)<1e-12,"自测2卡方不合格，真实卡方=%.6f，算出卡方=%.6f"%(chi2_hand,chi2_got)
    print("mcstat 自测通过")
        
