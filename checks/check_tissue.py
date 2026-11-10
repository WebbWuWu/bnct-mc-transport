# -*- coding: utf-8 -*-
# ============================================================
#  check_tissue.py  --  追查 ① 第 5 步（软组织 D_H 深处低约 20%，换成 kerma 记账以后还差不差）
#                       + 技术报告 4.5（M1 D_B 12.75–14.25 cm 低 2–3%，新发生器下还在不在）
#  判据归 Claude，2026-11-06 写（跑之前写下门槛）。运行：Windows Ctrl+F5（只读 json，几秒钟）。先跑完 explore\m1_tissue_runs.py。
#
#  比的是同一个量（纪律 2）——都是「每个源中子」：
#    你的键（两个种子 1010、1200 平均）      OpenMC（M1 的 10 个种子平均：openmc_kur + _s2、_s5–_s12，每个 100 批 × 1 万）
#    track_sum / track_sq                  kur flux                       cm
#    Bc_sum / Bc_sq                        kur B10 n,a × Q_val            eV（Q_val = 2341569 eV，同 check_M2_3）
#    kerma_sum / kerma_sq [H1、C12、O16]   kur H1 / C12 / O16 heating     eV（两边都是 径迹 × N × MT=301；H1 两边都含 (n,γ) 氘核反冲，不用再扣）
#  OpenMC 的 SE 用 10 个种子之间的散布算（不用它自报的）：深处 10 层加起来的 SE 也是直接从 10 个「和」的散布来的，层间相关已经在里面。
#  你的深处之和的 SE 没法精确算（run_slab 只存逐层的 Σx²，不存层与层的交叉项），所以给两个界：
#    SE 小 = √Σse²（假设层间不相关）· SE 大 = Σse（假设完全相关）→ 下结论用 SE 大（z 小，保守）
#
#  门槛（追查方案 §1，10/30 定；这里照搬，不因看到数据再改）
#    深处 = 第 21–30 层（下标 20–29，10–15 cm），10 层加起来比：r = 你 / OpenMC，z = 差 / √(SE你² + SE_OMC²)
#    |z| < 3 → 没看到偏差 · |z| ≥ 3 且 |r − 1| > 0.5% → 偏差还在 · |z| ≥ 3 但 |r − 1| ≤ 0.5% → 记下不追
#    追查 ① 第 5 步：H1 kerma 深处仍差 > 10%（r < 0.90）→ 做第 6 步（按核素拆：下面 C12、O16 已经一起打出来了）
#    技术报告 4.5：D_B 第 26–29 层（下标 25–28，12.5–14.5 cm）加起来，同一套门槛
#
#  只看（不是判据）：同一次模拟里 M1 的老记法（DH_sum：碰撞估计器、只记弹性反冲 H）vs kerma 记法扣掉 (n,γ) 氘核反冲
#    （kerma_H − 1315.8 eV × gH_sum）—— 两边都是你的程序、同一批历史，差别 = 「记账方法」本身，正是追查方案 §3 怀疑的那部分
# ============================================================
# ---- 路径（2026-11-10 整理文件夹，Claude 加）：本文件在 my_code 的子文件夹里。
#      下面几行把 my_code 根目录和各子文件夹加进搜索路径，后面的 import 才找得到别的文件。----
import os, sys
_MY_CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录
for _sub in ("whiteboard", "lessons", "openmc_compare", "checks", ""):
    _p = os.path.join(_MY_CODE, _sub) if _sub else _MY_CODE
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json
import math
import os

from mcstat import mean_se

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # my_code 根目录（11/10：本文件挪进了子文件夹，所以多退一层）
OUT_DIR = os.path.join(HERE, "outputs")
NL = 30
DZ = 0.5
N_HIST = 100000
SEEDS = (1010, 1200)
OMC_FILES = ["openmc_kur.json", "openmc_kur_s2.json"] + ["openmc_kur_s%d.json" % s for s in range(5, 13)]
Q_VAL = 2341569.0
E_RECOIL_NG = 1315.8          # eV，H(n,γ) 每次氘核反冲（技术报告 4.3 第 1 条，库里量的）
DEEP = list(range(20, 30))    # 第 21–30 层
DB_TAIL = list(range(25, 29)) # 第 26–29 层（12.5–14.5 cm，M1 那 4 层）
Z_GATE = 3.0
R_GATE = 0.005
R_MAX = 0.2

# (显示名, 你的 sum 键, 你的 sq 键, 子键, OpenMC tally, 系数)
TABLE = [
    ("通量",      "track_sum", "track_sq", None,  "kur flux",        1.0),
    ("D_B",       "Bc_sum",    "Bc_sq",    None,  "kur B10 n,a",     Q_VAL),
    ("H1 kerma",  "kerma_sum", "kerma_sq", "H1",  "kur H1 heating",  1.0),
    ("C12 kerma", "kerma_sum", "kerma_sq", "C12", "kur C12 heating", 1.0),
    ("O16 kerma", "kerma_sum", "kerma_sq", "O16", "kur O16 heating", 1.0),
]


def pick(d, key, sub):
    v = d[key]
    return v[sub] if sub is not None else v


def mine_one(d, sk, qk, sub):
    """一个种子：逐层均值、SE。"""
    su, sq = pick(d, sk, sub), pick(d, qk, sub)
    m, s = [], []
    for k in range(NL):
        mk, _, sk_ = mean_se(su[k], sq[k], N_HIST)
        m.append(mk)
        s.append(sk_)
    return m, s


def omc_ref(omcs, tally, coe):
    """OpenMC 10 个种子：逐层平均和 SE（种子间散布 / √n）；再给一个函数算任意几层之和的平均和 SE。"""
    n = len(omcs)
    per = [[x * coe for x in o[tally]["mean"]] for o in omcs]
    m = [sum(p[k] for p in per) / n for k in range(NL)]
    s = [math.sqrt(sum((p[k] - m[k]) ** 2 for p in per) / (n - 1) / n) for k in range(NL)]

    def block(layers):
        sums = [sum(p[k] for k in layers) for p in per]
        mb = sum(sums) / n
        sb = math.sqrt(sum((x - mb) ** 2 for x in sums) / (n - 1) / n)
        return mb, sb
    return m, s, block


def verdict(r, z):
    if abs(z) < Z_GATE:
        return "没看到偏差"
    return "偏差还在" if abs(r - 1) > R_GATE else "显著但 ≤ 0.5%，记下不追"


if __name__ == "__main__":
    mines = []
    for seed in SEEDS:
        p = os.path.join(OUT_DIR, "tissue_kerma_N%d_s%d.json" % (N_HIST, seed))
        assert os.path.exists(p), "找不到 %s —— 先跑 explore\\m1_tissue_runs.py" % p
        with open(p, "r") as f:
            mines.append(json.load(f))
    omcs = []
    for name in OMC_FILES:
        with open(os.path.join(OUT_DIR, name), "r") as f:
            o = json.load(f)
        # [Claude assert] 抓：混进了开 S(α,β) 的、或粒子数不一样的 OpenMC 结果
        assert o["batches"] == 100 and o["particles"] == 10000, "%s 不是 100 批 × 1 万" % name
        omcs.append(o)
    print("你的：tissue_kerma 种子 %s（各 %d 条）   vs   OpenMC M1 %d 个种子（共 %d 万）\n"
          % ("、".join(str(s) for s in SEEDS), N_HIST, len(omcs), len(omcs) * 100))

    summary = []
    for name, sk, qk, sub, tal, coe in TABLE:
        each = [mine_one(d, sk, qk, sub) for d in mines]
        m1 = [sum(e[0][k] for e in each) / len(each) for k in range(NL)]
        s1 = [math.sqrt(sum(e[1][k] ** 2 for e in each)) / len(each) for k in range(NL)]
        m2, s2, block = omc_ref(omcs, tal, coe)

        print("【%s】逐层：z 和 你/OMC（「--」= 你那边 R > 0.2）" % name)
        for row0 in range(0, NL, 10):
            zs, rs = [], []
            for k in range(row0, row0 + 10):
                if m1[k] <= 0 or s1[k] / m1[k] > R_MAX or m2[k] <= 0:
                    zs.append("%6s" % "--")
                    rs.append("%6s" % "--")
                else:
                    zs.append("%+6.2f" % ((m1[k] - m2[k]) / math.sqrt(s1[k] ** 2 + s2[k] ** 2)))
                    rs.append("%6.3f" % (m1[k] / m2[k]))
            print("  %4.1f–%4.1f cm  z      %s" % (row0 * DZ, (row0 + 10) * DZ, " ".join(zs)))
            print("  %12s  你/OMC %s" % ("", " ".join(rs)))

        blocks = [("深处 10–15 cm", DEEP)]
        if name == "D_B":
            blocks.append(("12.5–14.5 cm（4.5 那 4 层）", DB_TAIL))
        for bname, layers in blocks:
            a = sum(m1[k] for k in layers)
            se_lo = math.sqrt(sum(s1[k] ** 2 for k in layers))
            se_hi = sum(s1[k] for k in layers)
            b, sb = block(layers)
            r = a / b
            z_cons = (a - b) / math.sqrt(se_hi ** 2 + sb ** 2)
            z_lib = (a - b) / math.sqrt(se_lo ** 2 + sb ** 2)
            per_seed = []
            for e in each:
                per_seed.append(sum(e[0][k] for k in layers) / b)
            v = verdict(r, z_cons)
            print("  %s 加起来：r = %.4f（两个种子各 %s）· z = %+.2f（保守）/ %+.2f（层不相关）· R你 = %.4f → %s"
                  % (bname, r, " / ".join("%.4f" % x for x in per_seed), z_cons, z_lib, se_hi / a, v))
            summary.append((name, bname, r, z_cons, v))
        print()

    # ---- 只看：记账方法本身差多少（同一次模拟，两种记法）----
    print("—— 只看，不是判据 ——  同一次模拟里 D_H 两种记法（深处 10 层之和；两个种子平均）")
    if all("DH_sum" in d for d in mines):
        dh = sum(sum(d["DH_sum"][k] for k in DEEP) for d in mines) / len(mines) / N_HIST
        kh = sum(sum(d["kerma_sum"]["H1"][k] for k in DEEP) for d in mines) / len(mines) / N_HIST
        gh = sum(sum(d["gH_sum"][k] for k in DEEP) for d in mines) / len(mines) / N_HIST
        print("  老记法 DH（碰撞估计器、弹性反冲）= %.5g eV · kerma H1 − 1315.8 × (n,γ) = %.5g eV · 比 = %.4f"
              % (dh, kh - E_RECOIL_NG * gh, dh / (kh - E_RECOIL_NG * gh)))
    else:
        print("  json 里没有 DH_sum，跳过")

    print("\n" + "#" * 60)
    print("总结（门槛：|z| ≥ 3 且 |r − 1| > 0.5% 才算「偏差还在」；z 用保守的那个）")
    for name, bname, r, z, v in summary:
        print("  %-10s %-26s r = %.4f  z = %+6.2f  → %s" % (name, bname, r, z, v))
    hk = [x for x in summary if x[0] == "H1 kerma"][0]
    print("\n追查 ① 第 5 步：H1 kerma 深处 r = %.3f → %s" % (hk[2], "仍差 > 10%，做第 6 步（看上面 C12、O16）" if hk[2] < 0.90 else "差 ≤ 10%：20% 里大部分是记账方法 / 旧发生器造成的，看「只看」那一行"))
    db = [x for x in summary if x[0] == "D_B" and x[1].startswith("12.5")][0]
    print("技术报告 4.5：D_B 12.5–14.5 cm r = %.4f、z = %+.2f → %s" % (db[2], db[3], db[4]))
