from xslib import total_xs
import matplotlib.pyplot as plt
import numpy as np
import os

# 路径一律从本文件自己的位置算起，不依赖启动时的工作目录。
HERE = os.path.dirname(os.path.abspath(__file__))

T="294K"
E_b, tot_b, at_b, i0_b,awr_b = total_xs("B10", T)
E_h, tot_h, at_h, i0_h ,aer_h= total_xs("H1",  T)

k = i0_b
p = np.log(tot_b[k+1]/tot_b[k]) / np.log(E_b[k+1]/E_b[k])
print("B-10 热能区 log-log 斜率 p = %.6f" % p)
assert abs(p + 0.5) < 1e-2, "p = %.6f，不是 -1/2" % p

plt.plot(E_b, tot_b, label="B-10")        # 第一个参数是横轴，第二个是纵轴
plt.plot(E_h, tot_h, label="H-1")

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Energy (eV)")
plt.ylabel("Total cross section (barn)")
plt.legend()
plt.grid(True, which="both", alpha=0.3)   # which="both" 连次刻度线也画，对数图上有用

plt.savefig(os.path.join(HERE, "total_xs.png"), dpi=150)
plt.show()

k = i0_b
print("E[k]     = %.6e" % E_b[k])
print("E[k+1]   = %.6e" % E_b[k+1])
print("tot[k]   = %.6f" % tot_b[k])
print("tot[k+1] = %.6f" % tot_b[k+1])
print("分子 log(tot比) = %.6e" % np.log(tot_b[k+1]/tot_b[k]))
print("分母 log(E比)   = %.6e" % np.log(E_b[k+1]/E_b[k]))

print("tot_h[0]    = %.4f barn   (E=%.3e eV)" % (tot_h[0], E_h[0]))
print("tot_h[i0_h] = %.4f barn   (E=%.3e eV)" % (tot_h[i0_h], E_h[i0_h]))