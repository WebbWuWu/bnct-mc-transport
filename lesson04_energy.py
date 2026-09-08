# 写完这一段，扫一遍：
#   [ ] 有没有用 min/max/sum/m/alpha 这类名字覆盖掉还要用的东西
#   [ ] 函数调用：括号有没有？参数给了没有？
#   [ ] 函数返回几个值？接了几个？
#   [ ] assert 有没有漏 abs()？消息里的占位符对不对？
#   [ ] 同一个量有没有写两遍公式
#   [ ] 累加器初始化在哪一层循环外面；判据在哪一层循环外面
#   [ ] 方差公式里的均值，和平方和是同一个量吗
#   [ ] 数组下标有没有夹边界
#   [ ] rng 在哪一层创建？绝不能在循环体内或调用实参里
#   [ ] 这个累加器是 per-sample 还是 per-history？两者不能共用一套变量名
#   [ ] 这条判据是哪一档强度：恒等式(1e-12) / 统计量(3σ) / 工程估计(±%)
#   [ ] 容差是常数还是由量级导出的？统计量里有没有混进随手设的 E0
#
# 2026-09-05 重排：只搬动，未改任何表达式。分块见下方 ===== 标题。
import math
import mcstat

a=1664525
c=1013904223
m=2**32

class lcg:
    def __init__(self,seed):
        self.state=seed
    def uint(self):
        self.state=(self.state*a+c)%m
        return self.state
    def random(self):
        return self.uint()/m

def scatter(E, A, rng):
  mu_cm=rng.random()*2-1
  E_new=E*((A**2+2*A*mu_cm+1)/(A+1)**2)
  assert A**2+2*A*mu_cm+1>0,\
    "mu_lab不存在"
  mu_lab=(A*mu_cm+1)/math.sqrt(A**2+2*A*mu_cm+1)
  return E_new,mu_lab


def alpha(A):
  return (A-1)**2/(A+1)**2

def Ex2(alpha):
  if alpha==0:
    return 2
  else:
    return (2-alpha*((math.log(alpha))**2-2*math.log(alpha)+2))/(1-alpha)

N=1000000
N_hist=50000
N_SCAT=50
E_in=1.0
E0  = 2.0e6
E_TH= 0.025
N_hist_B=20000
print("N=%d  N_hist=%d  N_SCAT=%d  N_hist_B=%d  seed=2026  E_in=%.1f  E0=%.3e eV  E_TH=%.3f eV"
      %(N,N_hist,N_SCAT,N_hist_B,E_in,E0,E_TH))
for A in [1.0,12.0]:

  # ===== 块 0：只依赖 A 的解析量（不碰 rng，不碰任何循环） =====
  rng=lcg(2026)
  al=alpha(A)
  M=20
  E_j=N/M
  if al==1:
    xi=0.0
  elif al==0:
    xi=1
  else:
    xi=1+al*math.log(al)/(1-al)
  L=math.log(E0/E_TH)
  Naive=L/xi
  Renew=(L+Ex2(al)/2/xi)/xi

  # ===== 块 1：循环 1（样本 = 一次散射）验证 1 在循环内 =====
  r_sum=0.0
  r_min=2.0
  r_max=-1.0
  r_sq=0.0
  counts=[0]*M
  x_sum=0.0
  x_sq=0.0
  y_sum=0.0
  y_sq=0.0
  for i in range(N):
    E_new,mu=scatter(E_in,A,rng)
    r_new=E_new/E_in
    if r_new>=r_max:
      r_max=r_new
    if r_new<=r_min:
      r_min=r_new
    assert al<=r_new<=1,\
      "比值不合理"
    t=(r_new-al)/(1-al)
    j=int(t*M)
    if j >= M: j = M-1
    counts[j]+=1
    r_sum+=r_new
    r_sq+=r_new**2
    x = -math.log(r_new)
    x_sum += x
    x_sq  += x*x
    y_sum+=mu
    y_sq+=mu**2

  # ===== 块 2：循环 1 收尾 + 验证 2 / 3 / 5 / 能量比均值 =====
  chi2 = mcstat.chi2_uniform(counts, E_j)
  mean,  var,  se   = mcstat.mean_se(r_sum, r_sq, N)
  x_mean,x_var,x_se = mcstat.mean_se(x_sum, x_sq, N)
  y_mean,y_var,y_se = mcstat.mean_se(y_sum, y_sq, N)
  print()
  print("A=%.0f  alpha=%.6f  xi=%.6f  L=%.6f  Naive=%.6f"%(A,al,xi,L,Naive))
  mcstat.header()
  print("  %-12s [%.6f, %.6f]  in [%.6f, %.6f]   OK"%("1 range",r_min,r_max,al,1.0))
  print("  %-12s %14.6f %16s %11s %10s"%("2 chi2",chi2,"<43.82 df=19","-","E[chi2]=19"))
  mcstat.check("3 xi",       x_mean, xi,        x_se)
  mcstat.check("5 mu_lab",   y_mean, 2/3/A,     y_se)
  mcstat.check("- E'/E mean",mean,   (1+al)/2,  se)
  assert chi2<43.82,\
      "x**2出错,chi2=%.6f"%chi2

  # ===== 块 3：循环 4a（样本 = 一条 N_SCAT 次链） =====
  r50_sum=0.0
  r50_sq=0.0
  for n in range(N_hist):
    E_50=E0
    for step in range(N_SCAT):
      E_50,mu50=scatter(E_50,A,rng)
    r50=math.log(E0/E_50)
    r50_sum+=r50
    r50_sq+=r50**2

  # ===== 块 4：4a 收尾 + 验证 4a（Wald，3σ）/ 4a'（方差恒等式，工程 ±3%） =====
  r50_mean,r50_var,r50_se = mcstat.mean_se(r50_sum, r50_sq, N_hist)
  mcstat.check("4a Wald", r50_mean, N_SCAT*xi, r50_se)
  print("  %-12s %14.6f %16.6f %11s %+9.2f%%"
      %("4a' var", r50_var/N_SCAT/x_var, 1.0, "-", (r50_var/N_SCAT/x_var-1)*100))
  assert abs(r50_var/N_SCAT/x_var-1)<0.03,\
    "r50_var与N_SAT*x_var比值不合格，比值=%.6f"%(r50_var/N_SCAT/x_var)

  # ===== 块 5：循环 4b（样本 = 一条到热能的变长链） =====
  n_sum=0.0
  n_sq=0.0
  for i1 in range(N_hist_B):
    E3=E0
    n=0
    while E3>E_TH:
      E3,mu_cm3=scatter(E3,A,rng)
      n+=1
      assert n<100000,\
        "4b慢化未收敛，n=%.6f,E=%.6f"%(n,E3)
    n_sum+=n
    n_sq+=n**2

  # ===== 块 6：4b 收尾 + 验证 4b（普适界，工程）/ 4b'（更新理论，3σ） =====
  n_mean,n_var,n_se = mcstat.mean_se(n_sum, n_sq, N_hist_B)
  print("  %-12s %14.6f  band=[%.4f, %.4f]  (2/3..1 +-3SE)"
      %("4b over", n_mean-Naive, 2/3-3*n_se, 1+3*n_se))
  mcstat.check("4b' renew", n_mean, Renew, n_se)
  assert 2/3 - 3*n_se <= n_mean - Naive <= 1 + 3*n_se,\
    "n_mean有偏差偏大，Naive=%.6f"%Naive

  # ===== 块 7：打表（下一步要重写成统一表格：实测 / 解析 / σ 偏差） =====
  
