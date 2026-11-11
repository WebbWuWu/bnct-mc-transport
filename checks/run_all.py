# -*- coding: utf-8 -*-
# ============================================================
#  run_all.py  --  一条命令把「现在这版程序」的判据跑一遍（2026-11-10 Claude 写）
#
#  为什么要它：以前判据靠去掉 # 号手动跑。11/8 M3 交付核对时 T1、T3、判据 A 就没调到。
#
#  怎么跑（Windows，openmc 环境）：VS Code 里打开本文件 Ctrl+F5；或终端里
#        python checks\run_all.py            默认那一组（下面 QUICK）
#        python checks\run_all.py --all      再加上慢的、按课编号的老判据（要 outputs\ 里的结果文件）
#  每个判据是一个单独的进程：一个挂了不影响下一个。最后打印一张表；全过才打印「全部通过」。
#  完整输出存在 outputs\run_all_log.txt。
#
#  ⚠️ 11/10 写的时候 Claude 那边没有 h5py 和真库，只试过不读库的那几项和本文件的框架；
#     第一次在你电脑上跑的结果要发给 Claude 看。
# ============================================================
import os, sys, time, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))          # my_code\checks
ROOT = os.path.dirname(HERE)                               # my_code

# （名字，怎么跑，这条抓什么）。怎么跑：以 .py 结尾 = 直接跑那个文件；否则 = 本文件里 GROUPS 的一组函数
QUICK = [
    ("check_split",   "check_split.py",  "逐位回归：改结构 / 加开关以后，中子那条路的结果一位都没变"),
    ("check_rng",     "check_rng.py",    "随机数：跳跃 = 一步步走；每条历史的段不重叠"),
    ("thermal",       "group:thermal",   "S(α,β)：读表 T1 · 转方向 T2 · 抽 E′ T3E · 抽 E′+μ T3（各带自检）"),
    ("elastic_angle", "group:elastic",   "弹性角分布：读表（判据 A）· 抽 μ_cm（判据 B）（各带自检）"),
    ("check_lookup",  "check_lookup.py", "T5：开 S(α,β) 时 lookup_xs 换对了截面、总截面跟着改"),
    ("check_photon",  "check_photon.py", "P1：光子截面合成水 vs NIST XCOM"),
    ("check_emit",    "check_emit.py",   "G1、G2：生光子的个数、选线、方向、权重"),
]
SLOW = [
    ("check_kerma",     "check_kerma.py",     "K0–K2：按核素 kerma（径迹）vs 碰撞估计 / 逐道 Q 值"),
    ("check_L07_1",     "check_L07_1.py",     "第 7 课 1：kT → 0 退化到第 6 课"),
    ("check_L07_3_5",   "check_L07_3_5.py",   "第 7 课 3、5：sample_target"),
    ("check_L07_4",     "check_L07_4.py",     "第 7 课 4：截面大小没被动过"),
    ("check_L08_src",   "check_L08_src.py",   "第 8 课：源能量表格抽样"),
    ("check_L08_dose",  "check_L08_dose.py",  "第 8 课 3c、5b：Gy 换算、KUR 源表"),
    ("check_L09_3",     "check_L09_3.py",     "M1：软组织 vs OpenMC 逐层（读 outputs\\ 里的 json）"),
    ("check_L09_5",     "check_L09_5.py",     "M1：OpenMC 加 / 不加 S(α,β)"),
    ("check_M2_3",      "check_M2_3.py",      "M2：水模体 vs OpenMC（读 json）"),
    ("check_M3_seeds",  "check_M3_seeds.py",  "M3：按 10 个种子标定的形状检验（读 json）"),
    ("check_g5",        "check_g5.py",        "M4 G5：生出来的光子总权重 vs 径迹长度 γ 计数（配对，约 3–4 分钟；11/11 加）"),
]
# 不放进来的：check_L07_2（5–10 分钟）、check_L08_1_3（3–10 分钟）、check_L09_2 / check_L09_se（只能在 WSL 跑）、
#            check_L08_4（已过时：假库没有 "heat" 键，M2 以后就跑不了了）、check_M3_3 / _recal / _diag、check_tissue*（诊断用，结论已写进小结）


def group_thermal():
    from thermal import load_thermal, sample_thermal, rotate_dir, SAB_PATH, T_TAG
    from check_thermal import (check_T1, check_T1_selftest, check_T3E, check_T3E_selftest,
                               check_T3, check_T3_selftest, check_T2, check_T2_selftest)
    th = load_thermal(SAB_PATH, T_TAG)
    check_T1(th);                   check_T1_selftest(th)
    check_T2(rotate_dir);           check_T2_selftest()
    check_T3E(th, sample_thermal);  check_T3E_selftest(th, sample_thermal)
    check_T3(th, sample_thermal);   check_T3_selftest(th, sample_thermal)


def group_elastic():
    from elastic_angle import load_elastic_angle, sample_elastic_mu, LIB, NUCS
    from check_elastic import check_A, check_A_selftest, check_B, check_B_selftest
    angle = load_elastic_angle(LIB, NUCS)
    check_A(angle);                     check_A_selftest(angle)
    check_B(angle, sample_elastic_mu);  check_B_selftest(angle)


GROUPS = {"thermal": group_thermal, "elastic": group_elastic}


def run_one(how):
    """另开一个 Python 进程跑一项。返回 (退出码, 全部输出)。"""
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"          # 判据里打印 ✅、α 这些字；接管输出时 Windows 默认编码会报错，所以指定 utf-8
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if how.startswith("group:"):
        cmd = [sys.executable, os.path.abspath(__file__), "--group", how.split(":")[1]]
    else:
        cmd = [sys.executable, os.path.join(HERE, how)]
    p = subprocess.run(cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return p.returncode, p.stdout.decode("utf-8", errors="replace")


if __name__ == "__main__":
    if "--group" in sys.argv:                  # 被上面 run_one 调起来的那个进程走这里
        sys.path.insert(0, ROOT)
        sys.path.insert(0, HERE)
        GROUPS[sys.argv[sys.argv.index("--group") + 1]]()
        sys.exit(0)

    todo = QUICK + (SLOW if "--all" in sys.argv else [])
    os.makedirs(os.path.join(ROOT, "outputs"), exist_ok=True)
    log_path = os.path.join(ROOT, "outputs", "run_all_log.txt")
    rows = []
    with open(log_path, "w", encoding="utf-8") as log:
        for name, how, what in todo:
            print("跑 %-14s …" % name, end="", flush=True)
            t0 = time.time()
            code, out = run_one(how)
            dt = time.time() - t0
            # 退出码不是 0 = 有 assert 响了或程序报错；有的判据只打印 ❌ 不 assert，所以输出里有 ❌ 也算没过
            ok = (code == 0) and ("❌" not in out)
            rows.append((name, ok, dt, what))
            print(" %s（%.0f s）" % ("过" if ok else "没过 ❌", dt))
            log.write("=" * 30 + " %s（退出码 %d，%.0f s）\n%s\n" % (name, code, dt, out))
            if not ok:
                print("\n".join("      " + l for l in out.strip().splitlines()[-12:]))
    print("\n%-15s %-5s %6s  %s" % ("判据", "结果", "用时", "抓什么"))
    for name, ok, dt, what in rows:
        print("%-15s %-5s %5.0fs  %s" % (name, "过" if ok else "没过", dt, what))
    bad = [r[0] for r in rows if not r[1]]
    print("\n完整输出：%s" % log_path)
    if bad:
        print("没过的：%s —— 把上面的输出发给 Claude" % "、".join(bad))
        sys.exit(1)
    print("全部通过（%d 项）" % len(rows))
