# -*- coding: utf-8 -*-
# ============================================================
#  check_geom.py  --  M5 判据 G-1：三种面（ZPlane、Sphere、ZCylinder）的距离、法向、sense 对不对
#  判据归 Claude。2026-11-28 写。运行：Windows Ctrl+F5（openmc (3.10.20) 环境），不读库，几秒（Claude 这边 3 秒）。
#
#  被测：geometry.ZPlane / Sphere / ZCylinder 的 f、normal、distance(x,y,z,u,v,w,on_surf)，
#        geometry.sense_of(surf, x,y,z,u,v,w)
#  两条独立的路（纪律 2）：
#    程序那条：解二次方程 t² + 2kt + f = 0（圆柱 a·t² + 2kt + f = 0），按 sense / surf_on 选根
#    判据这条：**不解方程**。沿射线 r + t·u 每隔 h 直接算面方程 F（判据自己重写的 F，不调程序的 f），
#              找第一次变号，再二分到 10⁻¹⁵ 级。法向用 F 的中心差分。在面上的 sense 用判据自己的 ∇F。
#    比的是同一个量：「从 r 沿 u 到这个面的最短正路程（cm）」，碰不到 = inf。
#
#  四条：
#    [G-1a] 距离：手算的固定例子（讲义那几道 + 刚穿进 / 刚穿出 / 打不中 / 背着走 / 平行）
#           + 每个面 150 条随机射线（一般的点 75 条：球、圆柱一半在里面；刚穿过面 on_surf=True 75 条），
#           |t程序 − t判据| ≤ 1e-12·max(1, t)；一边 inf 另一边不是 = 错
#           抓：公式写错、选错根、忘了除以 a、k 用错量、on_surf 时没丢 t≈0 的根 / 丢错、w=0 或 a=0 没处理
#    [G-1b] 法向：normal vs F 的中心差分（步长 1e-5 cm），每个分量 |差| ≤ 1e-6·max(1, |∇F|)
#           抓：漏了 2 倍、分量串位、圆柱 z 分量不是 0
#    [G-1c] sense：离面远的点（离面 > 1e-6 cm）= F 的正负；离面 3e-10 cm（< EPS）的点，
#           故意放在「方向的反面」那一侧 → 必须按 u·∇F 的正负判，不能按 F 判
#           抓：else 少一支、在面上按 F 判、门槛除以 |∇F| 写错
#    [G-1d] 自检（纪律 4）：判据自己的正确模板必须不响；7 个故意写错的版本每个都必须响
#  门槛（看数据之前定，11/28）：
#    距离 1e-12 是框架 G-1 写的。依据：二分到区间 < 1e-15·t；F 在根附近的舍入误差 ~ 1e-16·R²，
#    换成距离 ~ 1e-16·R，离 1e-12 还有三个多数量级。法向用差分，截断误差对二次函数是 0，
#    舍入 ~ 1e-16·R²/h ~ 1e-10，门槛 1e-6 留四个数量级。
#  抓不到什么（照实写）：
#    · 两个交点挤在同一个 h（1e-4 cm）里的擦边情况：随机样本里概率 ~ 1e-6，不专门测
#    · EPS 取多大才合适（新近似候选，进附录 B 前问你）—— 这里只测「小于 EPS 时按方向判」
#    · 找 cell、穿面换 cell —— 归 G-2、G-3
# ============================================================
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import math
import random
import numpy as np

import geometry as G

INF = math.inf
TOL_T = 1e-12
TOL_N = 1e-6
H_MARCH = 1e-4        # cm，判据走射线的步长
L_MARCH = 120.0       # cm，平面：走这么远还没变号 = 碰不到（随机方向保证 |w| ≥ 0.1、点离面 ≤ 10 cm）
DELTA_ON = 3e-10      # cm，「贴着面」的点离面这么远（< EPS = 1e-9）
SEED = 52817


# ---------------- 判据自己的面方程（独立重写，不调程序） ----------------
def F_of(kind, p, x, y, z):
    if kind == "plane":
        return z - p[0]
    if kind == "sphere":
        return (x - p[0])**2 + (y - p[1])**2 + (z - p[2])**2 - p[3]**2
    return (x - p[0])**2 + (y - p[1])**2 - p[2]**2          # cyl


def gradF(kind, p, x, y, z):
    if kind == "plane":
        return 0.0, 0.0, 1.0
    if kind == "sphere":
        return 2*(x - p[0]), 2*(y - p[1]), 2*(z - p[2])
    return 2*(x - p[0]), 2*(y - p[1]), 0.0


def t_ref(kind, p, x, y, z, u, v, w, on_surf):
    """不解方程：沿射线走，找 F 第一次变号，再二分。碰不到返回 inf。"""
    t0 = 1e-6 if on_surf else 0.0          # 刚穿过的面：跳过 t≈0 那个根（真弦长都远大于 1e-6）
    # 走多远：球、圆柱是有界的，走过「离中心的距离 + R」还没碰到就永远碰不到（只为省时间，不影响判对错）
    if kind == "plane":
        L, h = L_MARCH, 1e-2               # 平面一条线最多交一次，步长大也不会漏
    elif kind == "sphere":
        L, h = math.sqrt((x-p[0])**2 + (y-p[1])**2 + (z-p[2])**2) + p[3] + 1.0, H_MARCH
    else:
        a = u*u + v*v
        if a == 0.0:
            return INF                     # 沿轴飞：离轴距离不变，永远碰不到
        L, h = (math.sqrt((x-p[0])**2 + (y-p[1])**2) + p[2] + 1.0)/math.sqrt(a), H_MARCH
    ts = np.arange(t0, L, h)
    X, Y, Z = x + ts*u, y + ts*v, z + ts*w
    if kind == "plane":
        Fv = Z - p[0]
    elif kind == "sphere":
        Fv = (X - p[0])**2 + (Y - p[1])**2 + (Z - p[2])**2 - p[3]**2
    else:
        Fv = (X - p[0])**2 + (Y - p[1])**2 - p[2]**2
    s0 = np.sign(Fv[0])
    if s0 == 0:
        s0 = np.sign(Fv[1])
    idx = np.nonzero(np.sign(Fv) != s0)[0]
    if len(idx) == 0:
        return INF
    i = idx[0]
    lo, hi = ts[i - 1], ts[i]
    flo = F_of(kind, p, x + lo*u, y + lo*v, z + lo*w)
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if mid <= lo or mid >= hi:
            break
        fm = F_of(kind, p, x + mid*u, y + mid*v, z + mid*w)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5*(lo + hi)


def make(mod, kind, p):
    if kind == "plane":
        return mod.ZPlane(p[0])
    if kind == "sphere":
        return mod.Sphere(p[0], p[1], p[2], p[3])
    return mod.ZCylinder(p[0], p[1], p[2])


def unit(a, b, c):
    n = math.sqrt(a*a + b*b + c*c)
    return a/n, b/n, c/n


def rand_dir(rng, kind):
    while True:
        mu = rng.uniform(-1, 1)
        phi = rng.uniform(0, 2*math.pi)
        s = math.sqrt(1 - mu*mu)
        u, v, w = s*math.cos(phi), s*math.sin(phi), mu
        if kind == "plane" and abs(w) < 0.1:
            continue                       # 太平的话交点在 L_MARCH 外面，判据那条路看不到
        if kind == "cyl" and u*u + v*v < 0.05:
            continue
        return u, v, w


def point_on(kind, p, rng):
    """面上的一个点（判据自己投影上去）。"""
    if kind == "plane":
        return rng.uniform(-6, 6), rng.uniform(-6, 6), p[0]
    if kind == "sphere":
        a, b, c = unit(rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1))
        return p[0] + p[3]*a, p[1] + p[3]*b, p[2] + p[3]*c
    ang = rng.uniform(0, 2*math.pi)
    return p[0] + p[2]*math.cos(ang), p[1] + p[2]*math.sin(ang), rng.uniform(-6, 6)


# ---------------- 用例 ----------------
GEOMS = {
    "plane":  [(-1.5,), (2.37,)],
    "sphere": [(0.5, -1.0, 2.0, 3.0), (1.0, 2.0, -1.0, 1.5)],
    "cyl":    [(0.4, -0.3, 2.5), (-1.3, 0.7, 1.15)],
}
U1 = (0.48, 0.6, 0.64)
FIXED = [   # (kind, p, x,y,z, u,v,w, on_surf, 手算答案；None = 用判据那条路)
    ("plane",  (-1.5,), 0.3, -0.7, 2.3, *unit(0.6, 0.69282, -0.4), False, None),
    ("plane",  (-1.5,), 0.3, -0.7, 2.3, 0.6, 0.8, 0.0, False, INF),                    # 平行
    ("plane",  (-1.5,), 0.3, -0.7, -1.5, *unit(0.6, 0.69282, -0.4), True, INF),        # 刚穿过
    ("sphere", (0.5, -1.0, 2.0, 3.0), 1.5, 0.5, 2.5, *U1, False, 1.196549671592048),  # 讲义：在里面
    ("sphere", (1.0, 2.0, -1.0, 1.5), 3.0, 3.5, 0.0, -0.48, -0.6, -0.64, False, 2.5 - math.sqrt(1.25)),  # 手算第 2 题
    ("sphere", (1.0, 2.0, -1.0, 1.5), 3.0, 3.5, 0.0, *U1, False, INF),                # 背着走
    ("sphere", (1.0, 2.0, -1.0, 1.5), 3.0, 3.5, 0.0, 0.6, -0.8, 0.0, False, INF),     # 打不中
    ("cyl",    (0.4, -0.3, 2.5), 1.0, 0.5, 7.0, *U1, False, 1.9525530802289675),      # 讲义：在里面
    ("cyl",    (0.4, -0.3, 2.5), 1.0, 0.5, 7.0, 0.0, 0.0, 1.0, False, INF),           # 沿轴（a = 0，故意的特殊值）
    ("cyl",    (0.4, -0.3, 2.5), 1.0, 0.5, 7.0, 0.0, 0.0, -1.0, False, INF),
    ("cyl",    (0.4, -0.3, 2.5), 4.0, 2.9, -1.3, 0.6, -0.8, 0.0, False, INF),         # 打不中
]


def close(a, b):
    if a == INF or b == INF:
        return a == b
    return abs(a - b) <= TOL_T*max(1.0, abs(b))


def run_suite(mod, verbose=False):
    """跑 G-1a、b、c，返回出错的条数和前几条说明。mod 要有 ZPlane、Sphere、ZCylinder、sense_of。"""
    bad = []
    n = {"a": 0, "b": 0, "c": 0}

    def fail(tag, msg):
        bad.append("[%s] %s" % (tag, msg))

    def dist(s, *args):
        try:
            return s.distance(*args)
        except Exception as e:                       # 除以 0 之类也算响
            return "报错 %s" % type(e).__name__

    # ---- G-1a 固定例子
    for kind, p, x, y, z, u, v, w, on, ans in FIXED:
        s = make(mod, kind, p)
        ref = ans if ans is not None else t_ref(kind, p, x, y, z, u, v, w, on)
        if ans is not None and ans != INF:           # 手算答案和判据那条路先互相核一遍（核的是判据本身）
            assert close(t_ref(kind, p, x, y, z, u, v, w, on), ans), "判据自己的射线法和手算对不上：%s %s" % (kind, p)
        got = dist(s, x, y, z, u, v, w, on)
        n["a"] += 1
        if isinstance(got, str) or not close(got, ref):
            fail("G-1a", "%s %s 点(%.3g,%.3g,%.3g) 方向(%.3g,%.3g,%.3g) on_surf=%s：程序 %s，应为 %s"
                 % (kind, p, x, y, z, u, v, w, on, got, ref))

    # ---- G-1a 随机射线 + G-1b + G-1c
    rng = random.Random(SEED)
    for kind, plist in GEOMS.items():
        for p in plist:
            s = make(mod, kind, p)
            for i in range(75):
                # 一般的点（里面或外面，离面不贴得太近）
                while True:
                    x, y, z = rng.uniform(-6, 6), rng.uniform(-6, 6), rng.uniform(-6, 6)
                    if kind != "plane" and i % 2 == 0 and F_of(kind, p, x, y, z) > 0:
                        continue                   # 球、圆柱：一半的点放在里面（里面、外面选根规则不一样，两边都要测够）
                    gx, gy, gz = gradF(kind, p, x, y, z)
                    gl = math.sqrt(gx*gx + gy*gy + gz*gz)
                    if abs(F_of(kind, p, x, y, z)) > 1e-6*gl:
                        break
                u, v, w = rand_dir(rng, kind)
                ref = t_ref(kind, p, x, y, z, u, v, w, False)
                got = dist(s, x, y, z, u, v, w, False)
                n["a"] += 1
                if isinstance(got, str) or not close(got, ref):
                    fail("G-1a", "%s %s 随机 #%d：程序 %s，判据 %s" % (kind, p, i, got, ref))
                # G-1b 法向 vs 中心差分
                try:
                    nx, ny, nz = s.normal(x, y, z)
                except Exception as e:
                    fail("G-1b", "%s normal 报错 %s" % (kind, type(e).__name__)); nx = ny = nz = float("nan")
                h = 1e-5
                fd = [(F_of(kind, p, x+h, y, z) - F_of(kind, p, x-h, y, z))/(2*h),
                      (F_of(kind, p, x, y+h, z) - F_of(kind, p, x, y-h, z))/(2*h),
                      (F_of(kind, p, x, y, z+h) - F_of(kind, p, x, y, z-h))/(2*h)]
                n["b"] += 1
                if not all(abs(a - b) <= TOL_N*max(1.0, gl) for a, b in zip((nx, ny, nz), fd)):
                    fail("G-1b", "%s %s 点(%.3g,%.3g,%.3g)：normal %s，差分 %s" % (kind, p, x, y, z, (nx, ny, nz), fd))
                # G-1c 离面远：按 F 的正负
                want = 1 if F_of(kind, p, x, y, z) > 0 else -1
                n["c"] += 1
                try:
                    got_s = mod.sense_of(s, x, y, z, u, v, w)
                except Exception as e:
                    got_s = "报错 %s" % type(e).__name__
                if got_s != want:
                    fail("G-1c", "%s %s 离面远的点：sense %s，应为 %s" % (kind, p, got_s, want))

                # 刚穿过面（on_surf=True）：点就在面上，随机方向
                x, y, z = point_on(kind, p, rng)
                u, v, w = rand_dir(rng, kind)
                ref = t_ref(kind, p, x, y, z, u, v, w, True)
                got = dist(s, x, y, z, u, v, w, True)
                n["a"] += 1
                if isinstance(got, str) or not close(got, ref):
                    fail("G-1a", "%s %s 刚穿过面 #%d：程序 %s，判据 %s" % (kind, p, i, got, ref))

                # G-1c 贴着面：故意放在方向的反面一侧，必须按方向判
                gx, gy, gz = gradF(kind, p, x, y, z)
                gl = math.sqrt(gx*gx + gy*gy + gz*gz)
                dot = u*gx + v*gy + w*gz
                want = 1 if dot >= 0 else -1
                side = -want                                   # 放到「F 的正负和方向相反」的那一侧
                xs, ys, zs = x + side*DELTA_ON*gx/gl, y + side*DELTA_ON*gy/gl, z + side*DELTA_ON*gz/gl
                n["c"] += 1
                try:
                    got_s = mod.sense_of(s, xs, ys, zs, u, v, w)
                except Exception as e:
                    got_s = "报错 %s" % type(e).__name__
                if got_s != want:
                    fail("G-1c", "%s %s 贴着面（离面 %.0e cm）：sense %s，应按方向为 %s" % (kind, p, DELTA_ON, got_s, want))
    if verbose:
        print("  G-1a 距离 %d 条 · G-1b 法向 %d 条 · G-1c sense %d 条" % (n["a"], n["b"], n["c"]))
    return bad


# ---------------- 自检用的正确模板 + 故意写错的版本 ----------------
class _Mod:
    pass


def _template():
    m = _Mod()
    EPS = G.EPS

    class ZPlane:
        def __init__(self, z0): self.z0 = z0
        def f(self, x, y, z): return z - self.z0
        def normal(self, x, y, z): return 0.0, 0.0, 1.0
        def distance(self, x, y, z, u, v, w, on_surf):
            if on_surf or w == 0.0: return INF
            t = (self.z0 - z)/w
            return t if t > 0 else INF

    class Sphere:
        def __init__(self, x0, y0, z0, R): self.x0, self.y0, self.z0, self.R = x0, y0, z0, R
        def f(self, x, y, z): return (x-self.x0)**2 + (y-self.y0)**2 + (z-self.z0)**2 - self.R**2
        def normal(self, x, y, z): return 2*(x-self.x0), 2*(y-self.y0), 2*(z-self.z0)
        def distance(self, x, y, z, u, v, w, on_surf):
            dx, dy, dz = x-self.x0, y-self.y0, z-self.z0
            k = dx*u + dy*v + dz*w
            if on_surf: return -2*k if k < 0 else INF
            f = dx*dx + dy*dy + dz*dz - self.R**2
            d = k*k - f
            if d < 0: return INF
            q = math.sqrt(d)
            if f < 0: return -k + q
            return INF if k >= 0 else -k - q

    class ZCylinder:
        def __init__(self, x0, y0, R): self.x0, self.y0, self.R = x0, y0, R
        def f(self, x, y, z): return (x-self.x0)**2 + (y-self.y0)**2 - self.R**2
        def normal(self, x, y, z): return 2*(x-self.x0), 2*(y-self.y0), 0.0
        def distance(self, x, y, z, u, v, w, on_surf):
            dx, dy = x-self.x0, y-self.y0
            a = u*u + v*v
            if a == 0.0: return INF
            k = dx*u + dy*v
            if on_surf: return -2*k/a if k < 0 else INF
            f = dx*dx + dy*dy - self.R**2
            d = k*k - a*f
            if d < 0: return INF
            q = math.sqrt(d)
            if f < 0: return (-k + q)/a
            return INF if k >= 0 else (-k - q)/a

    def sense_of(s, x, y, z, u, v, w):
        fv = s.f(x, y, z)
        nx, ny, nz = s.normal(x, y, z)
        if abs(fv) < EPS*math.sqrt(nx*nx + ny*ny + nz*nz):
            return 1 if u*nx + v*ny + w*nz >= 0 else -1
        return 1 if fv > 0 else -1

    m.ZPlane, m.Sphere, m.ZCylinder, m.sense_of = ZPlane, Sphere, ZCylinder, sense_of
    return m


def _wrongs():
    """7 个故意写错的版本（多是你写过或差点写的错）。"""
    out = []
    T = _template()

    # W1 圆柱忘了除以 a
    m = _template()
    class C1(T.ZCylinder):
        def distance(self, x, y, z, u, v, w, on_surf):
            return super().distance(x, y, z, u, v, w, on_surf) * (u*u + v*v) if (u*u + v*v) > 0 else INF
    m.ZCylinder = C1; out.append(("W1 圆柱忘了除以 a", m))

    # W2 球在里面取错根（−k − √）
    m = _template()
    class S2(T.Sphere):
        def distance(self, x, y, z, u, v, w, on_surf):
            if not on_surf and self.f(x, y, z) < 0:
                dx, dy, dz = x-self.x0, y-self.y0, z-self.z0
                k = dx*u + dy*v + dz*w
                t = -k - math.sqrt(k*k - self.f(x, y, z))
                return t if t > 0 else INF
            return super().distance(x, y, z, u, v, w, on_surf)
    m.Sphere = S2; out.append(("W2 球在里面取错根", m))

    # W3 on_surf 时取 t = 0 那个根（今天手算第 3 题的错）
    m = _template()
    class S3(T.Sphere):
        def distance(self, x, y, z, u, v, w, on_surf):
            return 0.0 if on_surf else super().distance(x, y, z, u, v, w, on_surf)
    m.Sphere = S3; out.append(("W3 刚穿过时取 t=0", m))

    # W4 不管 on_surf，靠噪声 f 选根（会选到 t≈0 → 卡死）
    m = _template()
    class S4(T.Sphere):
        def distance(self, x, y, z, u, v, w, on_surf):
            return super().distance(x, y, z, u, v, w, False)
    m.Sphere = S4; out.append(("W4 球不管 on_surf", m))

    # W5 平面没挡 w = 0（除以 0）
    m = _template()
    class P5(T.ZPlane):
        def distance(self, x, y, z, u, v, w, on_surf):
            if on_surf: return INF
            t = (self.z0 - z)/w
            return t if t > 0 else INF
    m.ZPlane = P5; out.append(("W5 平面没挡 w=0", m))

    # W6 sense 在面上也按 f 判（不看方向）；同时 else 只返回 −1（今天的错）也会被同一批点抓到
    m = _template()
    m.sense_of = lambda s, x, y, z, u, v, w: 1 if s.f(x, y, z) > 0 else -1
    out.append(("W6 在面上按 f 判 sense", m))

    # W7 圆柱 k 用了位置 y、z（今天的错）
    m = _template()
    class C7(T.ZCylinder):
        def distance(self, x, y, z, u, v, w, on_surf):
            dx, dy = x-self.x0, y-self.y0
            a = u*u + v*v
            if a == 0.0: return INF
            k = dx*y + dy*z
            if on_surf: return -2*k/a if k < 0 else INF
            f = dx*dx + dy*dy - self.R**2
            d = k*k - a*f
            if d < 0: return INF
            q = math.sqrt(d)
            if f < 0: return (-k + q)/a
            return INF if k >= 0 else (-k - q)/a
    m.ZCylinder = C7; out.append(("W7 圆柱 k 用了位置", m))
    return out


def check_G1():
    print("[G-1] 三种面：距离 / 法向 / sense（程序 vs 不解方程的射线法 + 手算）")
    bad = run_suite(G, verbose=True)
    for b in bad[:10]:
        print("   ", b)
    assert not bad, "G-1 没过：%d 条不对（上面列了前 10 条）" % len(bad)
    print("  ✅ G-1 过")


def check_G1_selftest():
    print("[G-1 自检] 正确模板必须不响，7 个写错的版本必须每个都响")
    bad = run_suite(_template())
    assert not bad, "判据自己的模板没过（判据写错了，是 Claude 的锅）：%s" % bad[:3]
    for name, m in _wrongs():
        b = run_suite(m)
        print("   %-22s 响了 %3d 条 %s" % (name, len(b), "✅" if b else "❌ 没响"))
        assert b, "自检：%s 没响 —— G-1 抓不住这类错" % name
    print("  ✅ G-1 自检过")


if __name__ == "__main__":
    check_G1()
    check_G1_selftest()
