# git 讲义（2026-09-08）

> 用法：**边读边在 VS Code 终端里敲**。不要一口气读完再动手。
> 每一节末尾的「判据」就是这一节做完了的标志。

---

## 0. 先确认 Windows 上有 git

VS Code 里打开终端（`Ctrl + ~`），敲：

```
git --version
```

- 打出 `git version 2.xx.x` → 有，跳到第 1 节。
- 报「不是内部或外部命令」→ 去 https://git-scm.com/download/win 装，一路默认，装完**重开 VS Code**（PATH 要重新加载）。

> ⚠️ 如果 `conda` 也报同样的错，那是 9/7 那个 PowerShell 执行策略问题，先跑：
> `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser`，然后重开终端。

**判据**：`git --version` 有输出。

---

## 1. git 是什么（这一节决定你用得对不对）

**git 不是网盘，不是备份工具。**

git 是一条**带说明的快照序列**。每次 `commit`，它把当前所有被跟踪文件的完整内容存成一个快照，附上一句你写的说明和一个时间戳，接在上一个快照后面。

于是你随时可以问三种问题：

| 问题 | 命令 |
|---|---|
| 这个文件从上次快照到现在，我改了什么？ | `git diff` |
| 这一串快照，每个都干了什么？ | `git log` |
| 三天前那一版长什么样？ | `git show <id>:<文件>` |

### 对你这个项目，git 的价值有三条，按重要性排：

**① `git diff` 是一个调试工具，而且是你现在最缺的那一类。**

你踩过的 bug 里有一整类是「**改坏了但不知道改了哪**」——9/7 你怕调试练习把 `lesson04_energy.py` 弄坏，是我重跑一遍核对的。有了 git，这个问题变成一条命令：`git diff lesson04_energy.py`。**改了什么，一行一行摆在你面前。**

而 9/5 的 `ker → xi` 只改了一半那种 bug，`git diff` 会把改过的每一行都高亮出来，你扫一眼就知道有没有漏。

**② 「跑不通了，但十分钟前是好的」有了退路。**

现在你没有退路。今天开始有。

**③ 12 月的叙事。**

老师点开你的 GitHub，看到的不是一坨最终代码，而是三个月里一次一次的提交记录。**这件事别人抄不走。** 代价是：前 4 课（8/30–9/7）的开发历史已经永久丢失，今天只能一次性导入成一个提交。**从今天起不再欠。**

---

## 2. 三个区（不懂这个，`add` 和 `commit` 的区别就永远是玄学）

```
   工作区                暂存区                  仓库
（你的文件夹里            （待提交清单）        （快照序列，.git/）
  实际的文件）
      |                      |                      |
      |--- git add 文件 ---> |                      |
      |                      |--- git commit -----> |
      |                                             |
      |<----------- git checkout（丢弃修改）---------|
```

**为什么要有暂存区？** 因为你一次干活常常改了五个文件，但它们属于**两件不同的事**。暂存区让你只把属于同一件事的那几个文件放进同一个提交。

> **一个提交 = 一件事。** 这是唯一一条需要背下来的纪律。
> 提交说明写不成一句话，就说明你该拆成两个提交。

`git status` 会同时告诉你这三个区的状态。**卡住的时候先敲 `git status`**，它几乎总会顺带告诉你下一步该敲什么。

---

## 3. 最小命令集（今天只需要这些）

```
git status                  # 我现在在什么状态（最常用，一天敲二十次）
git add <文件>              # 把某个文件放进暂存区
git add .                   # 把当前目录下所有改动放进暂存区
git commit -m "说明"        # 把暂存区打成一个快照
git log --oneline           # 看快照序列
git diff                    # 工作区 vs 暂存区：我改了但还没 add 的
git diff --staged           # 暂存区 vs 上个快照：我 add 了但还没 commit 的
git show <id>               # 看某个快照具体改了什么
```

`<id>` 是 `git log` 里那串十六进制，**打前 7 位就够**。

---

## 4. 第一次配置（一台电脑只做一次）

```
git config --global user.name "你的名字或GitHub用户名"
git config --global user.email "你GitHub账号的邮箱"
```

这两个值会写进**每一个**提交，公开可见。用 GitHub 上那个邮箱，否则 GitHub 认不出提交是你的（贡献格子不会绿）。

顺带两条对 Windows 有用的：

```
git config --global init.defaultBranch main
git config --global core.quotepath false
```

第二条的作用：**让 `git status` 正常显示中文文件名**。不设的话 `lesson04_原理讲义.pdf` 会被打成一串 `\346\234\254...` 转义码。（你目录里正好有中文名文件。）

**判据**：`git config --global --list` 能看到你刚设的四条。

---

## 5. `.gitignore` —— 建仓库前必须先写

**为什么必须在 `git init` 之后、`git add` 之前写**：一旦某个文件被提交过一次，它就永远留在历史里了，后面再加进 `.gitignore` 也删不掉那份历史。

在 VS Code 里：左侧文件树的 `my_code` 上右键 → 新建文件 → 文件名填 `.gitignore`（**点开头，没有后缀**，VS Code 里这样建没问题，Windows 资源管理器里建会被拒绝）。

内容：

```
__pycache__/
*.pyc
.vscode/
*.png
```

逐条解释，别抄了不懂：

| 行 | 挡的是什么 | 为什么 |
|---|---|---|
| `__pycache__/` | Python 的字节码缓存目录 | **它是从 .py 自动生成的**。凡是能自动生成的东西都不进仓库——不然每次跑程序 `git status` 都一堆噪声 |
| `*.pyc` | 同上，散落的单个文件 | 兜底 |
| `.vscode/` | 编辑器的本地配置 | 那是你这台机器的设置，不是项目的一部分 |
| `*.png` | 图片 | ⚠️ **这条今天下午可能要撤销**——画出来的通量图要不要进仓库，看第 8 节 |

> 你目录里那个 `__pycache__/` 我删不掉（挂载不给删权限），但没关系：**`.gitignore` 写在前面，它就从来不会被 git 看见**。

**判据**：写完 `.gitignore` 之后敲 `git status`，输出里**没有** `__pycache__`。

---

## 6. 建仓库 + 第一个提交

```
cd E:\BNCT_GPU\my_code
git init
git status
```

`git status` 现在会把所有文件列成 "Untracked files"（红色）。**先看一遍这个列表，确认里面没有你不想公开的东西。**

然后——**这里不要用 `git add .` 一把梭**。今天第一次，手动分成三个提交，练一下「一个提交 = 一件事」：

```
git add lesson01_lcg.py lesson01_flight.py lesson01_error.py
git commit -m "第1课：LCG、飞行距离抽样、误差棒覆盖率实验"

git add lesson02_direction.py lesson02_slab.py lesson02_slab_v3.py lesson02_rewrite_spec.md
git commit -m "第2课：各向同性方向抽样、1D平板 T/R/A（含无脚手架重写版）"

git add lesson03_tally.py lesson03_N1e6.txt lesson04_energy.py lesson04_spec.md lesson04_原理讲义.pdf
git commit -m "第3-4课：径迹长度估计器与误差棒；散射运动学与慢化，九条验证"
```

剩下的（白板题、讲义、`energy1.py`）：

```
git add whiteboard_02_spec.md whiteboard_05_spec.md whiteboard_05_scaffold.py tool_debugger_guide.md tool_git_guide.md .gitignore
git commit -m "白板题规格书、调试器与git讲义、.gitignore"
```

`energy1.py` **先别提交**，见第 7 节。

然后：

```
git log --oneline
```

**判据**：`git log --oneline` 打出 4 行；`git status` 说 "nothing to commit" 只剩 `energy1.py` 一个 untracked。

---

## 7. `energy1.py` 的去留（你来定）

它是 9/7 调试器练习的副本：N 砍小了，`Ex2` 已经改回 `return 2`。

**我的建议：留，但改名 + 加文件头。**

理由：它是「**故意把一个已知 bug 种回去，然后用调试器抓出来**」这件事的物证。12 月被问「你怎么定位 bug」时，这个文件配上进度日志里那段 `+7.47σ` 的记录，是一个完整的故事。删了就只剩口头描述。

如果留，做两件事：

1. 改名 `debug_practice_energy.py`（`energy1` 这个名字三个月后你自己都不知道是什么）
2. 文件第一行加：
   ```python
   # 2026-09-07 VS Code 调试器练习用的副本，不是主线代码。
   # 主线是 lesson04_energy.py。这里 N 被砍小以便单步。
   ```

改完再 `git add` + `commit -m "调试器练习存档"`。

不想留就直接删，也是干净的选择——**但要在 README 里说一句，别让读的人以为你藏了什么。**

---

## 8. 推到 GitHub

在 GitHub 网页上 New repository：

- 名字建议 `bnct-mc-transport`（不要叫 `my_code`，仓库名是简历的一部分）
- **Public**
- **不要**勾 "Add a README file"、不要选 .gitignore 模板、不要选 license
  （勾了会在远端先生成一个提交，和你本地的历史对不上，第一次 push 就会被拒绝——这是新手第一个坑）

建完 GitHub 会给你两行命令，长这样：

```
git remote add origin https://github.com/<你的用户名>/bnct-mc-transport.git
git branch -M main
git push -u origin main
```

第一次 push 会弹窗要你登录 GitHub 授权，照做。

**判据**：用**手机流量或无痕窗口**打开那个仓库链接（不能是你已登录的浏览器——那样私有仓库也看得见，验不出「公开」），文件列表在。

---

## 9. README.md（这是今天真正的产出）

`git init` 之前那些是操作，README 才是给人看的东西。

**README 不是流水账，是一份「我验证了什么」的清单。** 老师三十秒扫一眼，要能看出这个人写代码是带判据的。

结构建议（你自己写，不要我代笔——这是你的东西）：

```markdown
# BNCT 蒙特卡罗中子输运（Python 参考实现）

一句话说这是什么、为什么写（大创 / 保研）。

## 环境
Python 3.10（conda 环境 openmc），只用标准库 math（第5课起引入 numpy/h5py/matplotlib）。
所有程序 `python lessonXX_xxx.py` 直接跑，无参数，seed 固定 2026。

## 目录

| 文件 | 内容 | **验证了什么** |
|---|---|---|
| lesson01_lcg.py | LCG 伪随机数发生器 | 均匀性；误差棒覆盖率（400 种子收敛到 70.8%，理论 68%） |
| lesson02_slab_v3.py | 1D 平板 T/R/A，无脚手架重写 | 三个 case 的平均碰撞次数对 1/(1-c) |
| lesson03_tally.py | 径迹长度估计器 + 20 层通量 + 误差棒 | A：20 层逐层对解析解，全部 3σ 内；B：Σt·Σ(φ·dz) = 平均碰撞次数（跨课交叉检验） |
| lesson04_energy.py | 弹性散射运动学与慢化 | 九条：范围 / 卡方 / ξ / ⟨μ_lab⟩ / Wald / 方差恒等式 / 超越量普适界 / 更新理论 |

## 验证方法论
（这一节值钱：写你那四档判据强度，恒等式 1e-12 / 整数计数 / 统计量 3σ / 工程估计与普适界。
 再写一条「一条 assert 要有两个独立算出来的量」。）
```

⚠️ **README 里不要写虚的**（"高性能"、"高精度"）。**只写你验过的数字。** 每一句都要是你能被追问的。

写完提交：

```
git add README.md
git commit -m "README：各课验证清单与判据方法论"
git push
```

---

## 10. 从今天起的纪律

任务表里每天 ③ 收尾段那条 **git commit**，从今天起真的要做：

```
git add .
git status          # ← 提交前必看：确认没有多余的东西混进去
git commit -m "今天干了什么"
git push
```

提交说明写中文没问题。**写「这次改了什么」，不要写「更新」「修改」。**

好的：`第4课：4b判据从±5%换成超越量普适界[2/3,1]±3SE`
坏的：`update lesson04`

---

## 11. 三个新手坑（会撞到，提前认脸）

**① `git add` 之后又改了文件** → 那次改动**没在暂存区里**，commit 不会带上它。`git status` 会同时把这个文件列在 "Changes to be committed" 和 "Changes not staged" 两处。**再 `git add` 一次。**

**② `git commit` 没写 `-m`** → 弹进一个全屏编辑器（vim），看着像死机。按 `Esc`，输入 `:q!`，回车，退出来重新敲带 `-m` 的。

**③ 提交说明写错了，还没 push** → `git commit --amend -m "新说明"`。**已经 push 过就别改**，历史会分叉。
