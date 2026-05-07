# 维基
https://deepwiki.com/KaihangShi/Ternary-Phase-Diagram

# Ternary-Phase-Diagram

这段代码是基于**Flory-Huggins高分子溶液理论**，用于计算和绘制**三元相图（聚合物 - 溶剂 - 非溶剂/小分子）**的MATLAB脚本。

为了计算不同溶剂条件下的相线变化，你需要理解这四个脚本的功能，并知道如何修改热力学参数。以下是详细的解析与操作指南：

### 一、 代码功能深度解析

在Flory-Huggins理论中，通常有三个组分：
*   **组分 1 (phi1)**：非溶剂 / 小分子 (Nonsolvent / Small molecule)
*   **组分 2 (phi2)**：溶剂 (Solvent)
*   **组分 3 (phi3)**：聚合物 (Polymer)

这四个脚本分别计算三元相图的四个关键部分，它们都依赖于一个隐藏的核心函数 **`fun.m`**（该函数包含了化学势或吉布斯自由能的计算公式）。

#### 1. `src/main.m` (Binodal Line / 双节点线 / 浊点线)
*   **作用**：计算**双节点线（共存曲线）**和**结线（Tie-lines）**。这是实际相变发生的边界。
*   **原理**：给定稀相中的聚合物浓度（`ix3d`），通过非线性优化 `fmincon` 寻找浓相中的各组分浓度以及稀相中的非溶剂浓度，使得两相中各组分的**化学势相等**。
*   **变量解析**：
    *   `x(1)` 和 `x(2)`：代表浓相中的组分1和组分3的体积分数。
    *   `x(3)`：代表稀相中的组分1的体积分数。
    *   最后将两相的数据对（稀相和浓相）成对写入 `binodal_results.csv`。

#### 2. `Spinodal/main.m` (Spinodal Line / 旋节线)
*   **作用**：计算**旋节线**（系统处于绝对不稳定的边界）。
*   **原理**：给定聚合物浓度 `ix3`，利用 `fsolve` 求解混合吉布斯自由能的二阶导数（Hessian矩阵行列式）等于0的点。
*   **注意**：这里用的是 `fsolve`（解方程），而不是 `fmincon`（求极值）。

#### 3. `Critical_point/main.m` (Critical Point / 临界点)
*   **作用**：计算双节点线和旋节线相切交汇的**临界点**。
*   **原理**：寻找吉布斯自由能三阶导数为0的点。

#### 4. `LLE/main.m` (Liquid-Liquid Equilibrium / 单一液液平衡线)
*   **作用**：不使用循环，针对某一个特定的初始猜测值 `x0`，单次计算一条结线（两相平衡点）。通常用于调试 `fun.m` 和测试初始猜测值是否容易收敛。

---

### 二、 如何使用它计算“不同溶剂条件”下的相线变化？

**“不同溶剂条件”在代码底层的体现是 Flory-Huggins 相互作用参数（$\chi$ 参数，Chi parameter）的改变。** 

系统中存在三个 $\chi$ 参数：
1.  $\chi_{12}$ (非溶剂-溶剂)
2.  $\chi_{23}$ (溶剂-聚合物)
3.  $\chi_{13}$ (非溶剂-聚合物，代码中已提取为全局变量 `iX13`)

#### 实操指导步骤：

#### 第一步：定位并修改 `fun.m` 中的参数
这段代码虽然主程序里只有 `iX13`，但 $\chi_{12}$、$\chi_{23}$ 以及分子体积比（链长 $m_1, m_2, m_3$）必定被硬编码在同级目录的 **`fun.m`** 文件中。
*   **如果要更换溶剂**：你需要查阅文献或使用基团贡献法（如Hansen溶解度参数）计算新的溶剂与聚合物之间的 $\chi_{23}$，以及溶剂与非溶剂之间的 $\chi_{12}$。
*   打开 `fun.m`，找到代表 $\chi_{12}$ 和 $\chi_{23}$ 的变量，并替换为你的新溶剂参数。
    *   *良溶剂 (Good solvent)*：$\chi_{23}$ 较小（通常 < 0.5）。
    *   *劣溶剂 (Poor solvent)*：$\chi_{23}$ 较大。

#### 第二步：在主程序修改 `iX13`
在主程序（如 `src/main.m`）中：
```matlab
% interaction parameter nonsolvent(1)-polymer(3)
iX13 = 1.7; % <-- 如果你的非溶剂和小分子没变，这个值可以保持；如果变了，需修改。
```

#### 第三步：调整初始猜测值 `x0`（最关键、最容易报错的一步）
当溶剂改变（即 $\chi$ 参数改变）后，相图的形状和位置会发生巨大偏移。**以前能够收敛的初始猜测值 `x0` 会导致 `fmincon` 报错或找不到解。**

*   **对于 Binodal (`src/main.m`)**:
    ```matlab
    % x0 =[浓相phi1; 浓相phi3; 稀相phi1]
    x0 = [0.1; 0.7; 0.25]; 
    ```
    如果你换了更强的溶剂，互溶区变大，你需要相应地微调 `x0`。你可以先用 `LLE/main.m` 试错几个 `x0`，看 `f` (残差) 是否能收敛到接近 0（比如 `1e-10`）。

#### 第四步：执行与调试顺序
1.  **先计算旋节线 (Spinodal)**：旋节线的数学性质更好（解方程即可），对初始值相对不那么敏感。先运行 `Spinodal/main.m`，你可以大致看出相分离区域的轮廓。
2.  **再计算双节点线 (Binodal)**：根据旋节线的轮廓，合理推测 Binodal 曲线中 `x0` 的分布，然后运行 `src/main.m`。
3.  **观察 `residual`**：查看生成的 `binodal_results.csv` 中的第一列 `residual`（残差）。如果残差很大（比如大于 `1e-4`），说明该点没有真正收敛，得到的数据是错的，需要缩小 `base` 步长或调整那一点的 `x0`。

#### 第五步：可视化
MATLAB运行完毕后，你会得到各个组分的体积分数序列。你需要借助MATLAB的第三方三元相图工具包（如 `ternplot`）或者将CSV导入 **Origin** 软件中，使用 Origin 的 Ternary Plot 功能绘制相图，直观对比不同溶剂（不同参数下）两相区的扩大或缩小。

**总结：你的工作重心不在于修改上面这四个脚本的逻辑，而在于修改 `fun.m` 中的热力学参数，并在报错不收敛时耐心调整 `x0`。**

# 修改过程


| Parameter | File | Line(s) | Description |
|---|---|---|---|
| `v1`, `v2`, `v3` | `fun.m` | 24–28 | Molar volumes of nonsolvent, solvent, polymer |
| `g23` | `fun.m` | 55 | Solvent–polymer interaction parameter |
| `p1`–`p5` | `fun.m` | 119–123 | Polynomial coefficients for nonsolvent–solvent `g12` |
| `iX13` | `main.m` | 24 | Nonsolvent–polymer interaction parameter |
| `x0` | `main.m` | 10 | Initial guess (adjust if solver doesn't converge) |

> 接下来需要拟合得出小分子-溶液的p1-p5参数，
>
> ## 基团数量收集
groups_Tol = {9: 5, 10: 1}  # Toluene: 5*AC, 1*ACCH3
groups_OXY = {9: 4, 10: 2}  # o-Xylene: 4*AC, 2*ACCH3
groups_L8-Bo = {1: 8, 2: 28, 1: 8, 3: 4, 9: 12, 11: 4, 18: 2, 25: 4, 39: 5, 49: 4} #L8-Bo