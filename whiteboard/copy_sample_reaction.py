def sample_reaction(sigma,total,xi,order):
    assert len(order)==len(sigma),"order和sigma长度不一致"
    n=len(sigma)
    assert total > 0, "total 必须为正，实际是 %r" % total
    assert abs(sum(sigma) / total - 1) < 1e-12, \
        "sigma 与 total 对不上：sum(sigma)=%.17g  total=%.17g  比值=%.17g" \
        % (sum(sigma), total, sum(sigma) / total)
    assert 0 <= xi < 1, "xi 越界或是 nan: %r" % xi
    
    acc=0.0
    
    for rank in range(n):
        k=order[rank]
        acc+=sigma[k]
        if xi<acc/total:
            return k
    assert 1.0 - acc / total < 1e-12, \
        "扫完未返回，且缺口远大于浮点残差：xi=%.17g  acc=%.17g  total=%.17g  缺口=%.3g" \
        % (xi, acc, total, 1.0 - acc / total)   
        
    for r in range(n-1,-1,-1):
        if sigma[order[r]]>0:
            return order[r]
    assert False, "sigma 全为零而 total>0，不可能：total=%.17g" % total
 
if __name__ == "__main__":
    sigma = [2.6, 0.0, 7.3, 0.45]
    total = sum(sigma)
    order = sorted(range(len(sigma)), key=lambda i: sigma[i], reverse=True)
    print("order =", order, "  total =", total)
    for xi in (0.37, 0.81, 0.97):
        k = sample_reaction(sigma, total, xi, order)
        print("xi = %.2f  ->  k = %d  (sigma[k] = %g)" % (xi, k, sigma[k]))   
