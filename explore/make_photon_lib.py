# -*- coding: utf-8 -*-
# ============================================================
#  explore/make_photon_lib.py  --  11/7 第 3 项（用户定丙）：给 lib80x 补 H、B、O 的光子数据
#  Claude 2026-11-07 写（工具脚本，不默写）。
#  做什么：① 从 E:\BNCT_GPU\photon_src\ 两个 NNDC zip（ENDF/B-VIII.0 光原子 + 原子弛豫子库）里只解出 H、B、O 六个文件
#          ② openmc.data.IncidentPhoton.from_endf(光原子, 原子弛豫) → lib80x_hdf5\photon\H.h5、B.h5、O.h5（官方库也是这个函数生成的）
#          ③ cross_sections.xml 先备份成 .bak_20261107，再在末尾 </cross_sections> 前加三行 type="photon"（原来 589 行不动；已经有就不加）
#          ④ 判据 P1：读回写好的 h5，1 MeV 康普顿（MT=504）截面 vs Klein–Nishina × Z（自由电子理论，和 ENDF/EPDL97 两条独立的路），|差| < 1%
#              抓什么：元素对错（Z 搞错 → 差几倍）、单位（b vs cm²/g）、文件没写进去 / 读回来是空的
#              参考：NIST XCOM 1 MeV 总 μ/ρ：H 0.1263、O 0.06372 cm²/g → 0.2114、1.693 b/原子（含相干 + 光电，比 KN 多不到 0.2%）
#  运行：先输 wsl → conda activate /home/siyiwu/miniconda3/envs/openmc → cd /mnt/e/BNCT_GPU/my_code → python explore/make_photon_lib.py
#  ⚠️ 只在 WSL 里跑（Windows 那个 openmc 环境里没有 OpenMC）。
# ============================================================
import os
import math
import shutil
import zipfile
import openmc.data

BASE = "/mnt/e/BNCT_GPU"
SRC = os.path.join(BASE, "photon_src")
LIB = os.path.join(BASE, "lib80x_hdf5")
OUT = os.path.join(LIB, "photon")
XML = os.path.join(LIB, "cross_sections.xml")
ELEMS = [(1, "H"), (5, "B"), (8, "O")]
ZIPS = {"photoat": "ENDF-B-VIII.0_photoat.zip", "atom": "ENDF-B-VIII.0_atomic_relax.zip"}
ZIPDIR = {"photoat": "ENDF-B-VIII.0_photoat", "atom": "ENDF-B-VIII.0_atomic_relax"}


def kn_per_electron(E_eV):
    """Klein–Nishina 总截面（每个自由电子，b）。"""
    re = 2.8179403262e-13            # 经典电子半径 cm
    k = E_eV / 510998.95             # 光子能量 / 电子静能
    a = ((1 + k) / k**2 * (2 * (1 + k) / (1 + 2 * k) - math.log(1 + 2 * k) / k)
         + math.log(1 + 2 * k) / (2 * k) - (1 + 3 * k) / (1 + 2 * k)**2)
    return 2 * math.pi * re**2 * a / 1e-24


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    # ① 只解出需要的 6 个文件
    files = {}
    for Z, sym in ELEMS:
        for kind in ["photoat", "atom"]:
            name = "%s/%s-%03d_%s_000.endf" % (ZIPDIR[kind], kind, Z, sym)
            with zipfile.ZipFile(os.path.join(SRC, ZIPS[kind])) as zf:
                zf.extract(name, SRC)
            files[(sym, kind)] = os.path.join(SRC, name)
    print("① 解出 6 个 ENDF 文件 → %s" % SRC)

    # ② 转 h5
    lines = []
    for Z, sym in ELEMS:
        data = openmc.data.IncidentPhoton.from_endf(files[(sym, "photoat")], files[(sym, "atom")])
        path = os.path.join(OUT, sym + ".h5")
        data.export_to_hdf5(path, "w")
        print("② %s：Z = %d，%d 个反应 → %s（%d B）" % (sym, data.atomic_number, len(data.reactions), path, os.path.getsize(path)))
        lines.append('  <library materials="%s" path="photon/%s.h5" type="photon" />\n' % (sym, sym))

    # ③ xml：备份 + 末尾加三行
    with open(XML, "r", encoding="utf-8") as f:
        txt = f.read()
    if 'type="photon"' in txt:
        print("③ xml 里已经有 photon 行，不改")
    else:
        shutil.copy(XML, XML + ".bak_20261107")
        n_old = txt.count("<library")
        txt = txt.replace("</cross_sections>", "".join(lines) + "</cross_sections>")
        with open(XML, "w", encoding="utf-8") as f:
            f.write(txt)
        print("③ xml：备份 → cross_sections.xml.bak_20261107；<library 行 %d → %d" % (n_old, txt.count("<library")))

    # ④ 判据 P1：读回来比 Klein–Nishina × Z
    E = 1.0e6
    kn = kn_per_electron(E)
    ok_all = True
    for Z, sym in ELEMS:
        back = openmc.data.IncidentPhoton.from_hdf5(os.path.join(OUT, sym + ".h5"))
        sig = float(back.reactions[504].xs(E))
        ref = Z * kn
        r = sig / ref
        ok = abs(r - 1.0) < 0.01
        ok_all = ok_all and ok
        print("④ P1 %s：1 MeV 康普顿 %.4f b vs KN × %d = %.4f b → 比 %.4f %s" % (sym, sig, Z, ref, r, "✅" if ok else "❌"))
    assert ok_all, "P1 不过：1 MeV 康普顿截面和 Klein–Nishina × Z 差 1% 以上 —— 把上面几行发 Claude"
    print("全部 ✅ —— 把整段输出发 Claude")
