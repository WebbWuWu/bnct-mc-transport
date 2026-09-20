# matplotlib 讲义（2026-09-08）

> 只讲你今天要用的。**边读边敲**，不要通读。

---

## 1. 心智模型：两层容器

matplotlib 里画图只有两个概念：

```
Figure（画布）  ——  一整张图片，最后存成 png 的就是它
  └── Axes（坐标系） ——  一套 x 轴 y 轴和里面画的东西
        └── 线、点、误差棒、图例……
```

一张 Figure 里可以有多个 Axes（就是「子图」）。今天只用一个。

**开场白永远是这一句：**

```python
fig, ax = plt.subplots(figsize=(7, 5))
```

它一次给你一个画布 `fig` 和一个坐标系 `ax`。`figsize` 的单位是英寸。

> ⚠️ **网上大量教程用 `plt.plot()` / `plt.xlabel()` 这种写法**（叫 pyplot 接口，它偷偷操作「当前那个图」）。
> **不要学那种。** 用 `ax.plot()` / `ax.set_xlabel()` 这种（叫面向对象接口）：**当你要画两个子图对比时，`plt.` 那套会立刻乱套，而 `ax` 这套天然清楚。** 第 9 课 OpenMC 对拍就是双子图（曲线 + 残差）。

---

## 2. 四个动作，一张图就出来了

```python
import matplotlib.pyplot as plt        # 约定俗成缩写 plt，别改

fig, ax = plt.subplots(figsize=(7, 5))

ax.errorbar(x, y, yerr=err, label="case 1")   # ① 画
ax.set_xlabel("z (cm)")                        # ② 标轴
ax.set_ylabel("flux (cm$^{-2}$ per source neutron)")
ax.legend()                                    # ③ 图例（读 label= 里的字）
fig.savefig("lesson03_flux.png", dpi=150, bbox_inches="tight")   # ④ 存
plt.show()
```

**`errorbar` 的常用参数**（`x, y, yerr` 之后的都可省）：

| 参数 | 作用 |
|---|---|
| `yerr=err` | 每个点的误差棒半长（**1σ**，你的 `err` 就是） |
| `label="case 1"` | 图例里显示的名字 |
| `marker="o"` | 数据点画成圆点；`"s"` 方块、`"^"` 三角 |
| `markersize=3` | 点的大小 |
| `capsize=2` | 误差棒两端那两道小横杠的长度，**0 就没有横杠** |
| `linewidth=1` | 连线粗细 |

**`$...$` 是数学模式**：`cm$^{-2}$` 会渲染成 cm⁻²。这是 LaTeX 语法，matplotlib 内置支持，不用装 LaTeX。

**`bbox_inches="tight"`**：把四周多余的白边裁掉。不加的话轴标签经常被切掉半截。

---

## 3. 三个必踩的坑

**① `savefig` 必须在 `show` 之前。**
`plt.show()` 弹窗关掉之后画布被清空，**再 savefig 存出来是一张白图**。顺序记死：**先存后看**。

**② 中文标签会变成一排方块。**
matplotlib 默认字体没有中文。**今天全部用英文标签**，别去折腾字体，那是个坑。

**③ 误差棒看不见 ≠ 没画上。**
你的 case 1 相对误差只有 1e-3 量级，误差棒比线还细。
**验法**：把 `yerr=err` 临时改成 `yerr=[e*100 for e in err]`，如果误差棒明显变大，说明它一直在，只是太小。**验完改回来。**

> 这是「**没见过它失败的 assert，不算验过**」在画图上的版本：**没见过它长出来的误差棒，不算画上了。**

---

## 4. 顺带：怎么把数据写进 csv

你选了 B 路，所以 lesson03 要多吐一个 csv。**用 `with open(...)`：**

```python
with open(path, "w", encoding="utf-8") as f:
    f.write("k,z_lo,z_hi,z_mid,phi,err,R\n")     # 表头
    f.write("%d,%.6f,%.6f,%.6f,%.8e,%.8e,%.8f\n" % (k+1, z_lo, z_hi, z_mid, phi, err, R))
```

三件事解释一下：

| 写法 | 为什么 |
|---|---|
| `with open(...) as f:` | 出了这个缩进块**自动关文件**，哪怕中途抛异常。不用 `with` 就要自己 `f.close()`，而一旦忘了，数据可能还卡在缓冲区里没落盘——**就是你今天上午撞的那个块缓冲** |
| `"w"` | 写模式，**会清空原文件**。追加是 `"a"`，读是 `"r"` |
| `encoding="utf-8"` | 不写的话 Windows 上默认 GBK，将来存中文会乱码。**Windows 上永远显式写 encoding** |
| `%.8e` 存 phi/err | **存文件要比打印多留几位**。屏幕上 `%.6f` 是给人看的，文件是给下一个程序吃的，精度别在这儿丢掉 |

⚠️ **`z_mid`（层中心）是新加的一列**：画图的 x 轴要用层中心，不是层的左边界，也不是层号。想一下为什么。

---

## 5. 读回来

```python
with open(path, "r", encoding="utf-8") as f:
    lines = f.readlines()          # 一次读成一个列表，每个元素是一行（含结尾的 \n）
```

`lines[0]` 是表头，要跳过。剩下每行 `line.strip().split(",")` 切成字符串列表，
再 `float(...)` 转成数。**csv 里存的是文字，`"0.951626"` 不是数字，不转不能画。**

（Python 有 `csv` 模块和 pandas，都更省事。今天先用最原始的写法，**因为你需要知道文件里到底是什么**。第 5 课读 ENDF 时会换成 h5py。）
