"""第二章 矩阵及其运算。"""

THEORY = r"""
## 知识点 2.1 矩阵的概念

**定义 2.1** $m\times n$ 个数排成的 $m$ 行 $n$ 列数表称为 $m\times n$ **矩阵**，
记作 $A=(a_{ij})_{m\times n}$。$m=n$ 时称**方阵**。
两个矩阵**同型**指行数、列数分别相等；同型矩阵 $A=B$ 指对应元素全相等。

常见特殊矩阵：**零矩阵**（元素全为零，记 $O$）；**单位矩阵**
$E=\operatorname{diag}(1,\ldots,1)$；**对角矩阵** $\operatorname{diag}(\lambda_1,\ldots,\lambda_n)$；
**数量矩阵** $kE$；**上（下）三角形矩阵**；只有一行的**行矩阵**（行向量）
与只有一列的**列矩阵**（列向量）。

为什么用矩阵：线性方程组的全部信息（系数与常数项）可压缩为增广矩阵，
线性变换也由矩阵完全刻画（第六章）。矩阵是"把一组数当作一个整体来运算"的语言。

## 知识点 2.2 矩阵的运算

**加法与数乘**（仅同型之间）：$(A+B)_{ij}=a_{ij}+b_{ij}$，$(kA)_{ij}=ka_{ij}$。
满足交换、结合律与两条分配律，证明都是按元素归结为数的运算律。
全体 $m\times n$ 矩阵在加法与数乘下构成第四章意义下的向量空间。

**乘法** 设 $A$ 为 $m\times s$ 矩阵，$B$ 为 $s\times n$ 矩阵
（$A$ 的列数必须等于 $B$ 的行数），则 $C=AB$ 为 $m\times n$ 矩阵，

$$
c_{ij}=\sum_{k=1}^{s}a_{ik}b_{kj},
$$

即 $C$ 的 $(i,j)$ 元是 $A$ 的第 $i$ 行与 $B$ 的第 $j$ 列对应元素乘积之和。

为什么这样定义：把线性方程组 $y_i=\sum_k a_{ik}x_k$ 与
$x_k=\sum_j b_{kj}t_j$ 复合，$y_i=\sum_j\big(\sum_k a_{ik}b_{kj}\big)t_j$，
复合变换的系数恰是上式。乘法是为"复合"而生的，不是逐个元素相乘。

**乘法的性质与三条警告**：结合律 $(AB)C=A(BC)$（由定义求和换序可证）、
分配律 $A(B+C)=AB+AC$ 成立；但
（1）**交换律一般不成立**：$AB$ 有意义时 $BA$ 未必有意义，即使都有意义也未必相等；
（2）$AB=O$ 推不出 $A=O$ 或 $B=O$；
（3）消去律不成立：$AB=AC$ 且 $A\ne O$ 推不出 $B=C$。
这三点都是因为乘法来自复合而非数乘，习题 1、4 会给出实例与补救条件。

**方阵的幂** $A^k=\underbrace{A\cdots A}_{k}$；注意 $(AB)^k\ne A^kB^k$（交换律缺失），
但若 $AB=BA$ 则二项式定理成立：$(A+B)^n=\sum_k \binom{n}{k}A^kB^{n-k}$。

**转置** $(A^T)_{ij}=a_{ji}$。性质：$(A^T)^T=A$，$(A+B)^T=A^T+B^T$，
$(kA)^T=kA^T$，$(AB)^T=B^TA^T$（最后一条由定义直接验算：
$(AB)^T$ 的 $(i,j)$ 元是 $(AB)$ 的 $(j,i)$ 元 $=\sum_k a_{jk}b_{ki}$，
$B^TA^T$ 的 $(i,j)$ 元 $=\sum_k (B^T)_{ik}(A^T)_{kj}=\sum_k b_{ki}a_{jk}$，相等）。
$A^T=A$ 称**对称矩阵**，$A^T=-A$ 称**反对称矩阵**（对角元必为零）。

**方阵的行列式** $n$ 阶方阵 $A$ 的行列式记 $|A|$ 或 $\det A$（第一章）。
性质：$|A^T|=|A|$，$|kA|=k^n|A|$（每行提出一个 $k$，共 $n$ 行），
以及关键的乘法定理：

**定理 2.1** $|AB|=|A|\cdot|B|$。（完整证明在第三章知识点 3.2 末尾，
用初等矩阵给出；二阶情形可由定义直接验算。）

## 知识点 2.3 逆矩阵

**定义 2.2** 设 $A$ 为 $n$ 阶方阵，若存在 $n$ 阶方阵 $B$ 使 $AB=BA=E$，
称 $A$ **可逆**，$B$ 为 $A$ 的**逆矩阵**。

**唯一性**：若 $B,C$ 都是 $A$ 的逆，则
$B=BE=B(AC)=(BA)C=EC=C$。故逆唯一，记 $A^{-1}$。

**定义 2.3（伴随矩阵）** 设 $A_{ij}$ 为 $a_{ij}$ 的代数余子式（第一章），
把矩阵 $(A_{ij})$ **转置**后得到的矩阵称为 $A$ 的**伴随矩阵** $A^*$，
即 $A^*$ 的 $(i,j)$ 元是 $A_{ji}$。**注意下标的转置，这是最易错的一步。**

**引理** $AA^*=A^*A=|A|E$。

证明：$AA^*$ 的 $(i,j)$ 元为 $\sum_k a_{ik}A_{jk}$。$i=j$ 时这是 $|A|$ 的按行展开
（定理 1.3），$i\ne j$ 时是异乘变零（定理 1.4），故 $AA^*=|A|E$；列的方向同理。$\square$

**定理 2.2（可逆判据）** $A$ 可逆当且仅当 $|A|\ne 0$，且此时

$$
A^{-1}=\frac{1}{|A|}A^*.
$$

证明：若 $|A|\ne0$，由引理 $A\cdot\frac{A^*}{|A|}=E$ 且
$\frac{A^*}{|A|}\cdot A=E$，按定义 $A$ 可逆且逆为 $\frac{A^*}{|A|}$。
反之若 $A$ 可逆，$AA^{-1}=E$，由定理 2.1，
$|A||A^{-1}|=|E|=1$，故 $|A|\ne0$。$\square$

$|A|=0$ 的方阵称为**奇异矩阵**，$|A|\ne0$ 称**非奇异矩阵**。

**逆的性质**（均设方阵可逆）：
$(A^{-1})^{-1}=A$；$(kA)^{-1}=\frac1k A^{-1}\ (k\ne0)$；
$(AB)^{-1}=B^{-1}A^{-1}$（验算：$(AB)(B^{-1}A^{-1})=A(BB^{-1})A^{-1}=E$，顺序不可颠倒）；
$(A^T)^{-1}=(A^{-1})^T$。

**推论（解矩阵方程）** 若 $A$ 可逆，则 $AX=B$ 有唯一解 $X=A^{-1}B$；
$XA=B$ 有唯一解 $X=BA^{-1}$。左乘、右乘不可混淆（交换律缺失的直接后果）。

求逆的两条路：低阶（$n\le3$）可用伴随矩阵公式；一般情形用第三章的初等行变换。

## 知识点 2.4 矩阵分块法

把矩阵用纵横线划成若干小矩阵块，以块为"元素"做运算，称为**分块矩阵**。
加法要求分法相同；乘法 $AB$ 要求 $A$ 的列分法与 $B$ 的行分法一致，
乘法公式形式与普通矩阵乘法相同（块的乘积不可交换次序）。

**分块对角矩阵** $\operatorname{diag}(A_1,\ldots,A_s)$（非对角块为零）：
行列式 $=|A_1|\cdots|A_s|$（反复用第一章展开或第三章结论）；
可逆当且仅当每块可逆，逆为 $\operatorname{diag}(A_1^{-1},\ldots,A_s^{-1})$。

按行分块给出方程组的向量形式：$A$ 按行分为 $A_1,\ldots,A_m$（行向量），
$Ax=b$ 即 $A_i x=b_i\ (i=1,\ldots,m)$；按列分块 $A=(\alpha_1,\ldots,\alpha_n)$，
则 $Ax=x_1\alpha_1+\cdots+x_n\alpha_n$——方程组有解当且仅当 $b$ 是
$A$ 的列向量组的线性组合，这是通往第四章的桥。

---

# 解题方法

## 题型一：矩阵乘法与方阵的幂

**方法 1（定义法）** $c_{ij}=$（$A$ 第 $i$ 行）$\cdot$（$B$ 第 $j$ 列）逐元算。
动手前先查型：$m\times s$ 乘 $s\times n$ 得 $m\times n$，型对不上说明题目或抄写有误。

**方法 2（按列/按行看法）** $AB$ 的第 $j$ 列 $=A\times(B$ 的第 $j$ 列$)$；
特别 $Ax$ 是 $A$ 的列向量组以 $x$ 的分量为系数的线性组合。
整体结构题、验算题用这个看法更快，也通向第四章。

**例 2.1** $A=\begin{pmatrix} 1 & 0 & 2 \\ 0 & 1 & 1 \end{pmatrix}$，
$B=\begin{pmatrix} 1 & 2 \\ 0 & 1 \\ 3 & 0 \end{pmatrix}$。
$AB$ 是 $2\times2$ 矩阵：第 1 列 $=A(1,0,3)^T=(1+6,\ 0+3)^T=(7,3)^T$
（即 $A$ 的第 1 列乘 1、第 3 列乘 3 相加），第 2 列 $=A(2,1,0)^T=(2,1)^T$。
故 $AB=\begin{pmatrix} 7 & 2 \\ 3 & 1 \end{pmatrix}$（与方法 1 逐元复核一致）。

**方阵的幂** 按结构选方法：（1）**拆分 $+$ 二项式**：拆成数量矩阵 $+$ 幂零或简单矩阵
（前提：两部分可交换），习题 6 是样板；（2）**找规律 $+$ 归纳证明**：先算 $A^2,A^3$ 猜形式再归纳；
（3）**秩一公式**：$A=uv^T$ 时 $A^2=(v^Tu)A$，故 $A^n=(v^Tu)^{n-1}A$；
（4）**对角化法**：一般情形的主力，在第五章题型三。共同纪律：展开高次幂时
每一步都不许交换因子次序，除非已验证可交换。

## 题型二：求逆矩阵

**方法 1（伴随矩阵法）** 适合 $n\le3$ 或元素含字母、需要公式形式时：
算 $|A|$（为零则不可逆，停）$\to$ 逐个算代数余子式 $\to$ **转置**得 $A^*$
$\to$ $A^{-1}=A^*/|A|$ $\to$ 抽查一行一列验算。二阶有速算公式
$\begin{pmatrix} a & b \\ c & d \end{pmatrix}^{-1}
=\dfrac{1}{ad-bc}\begin{pmatrix} d & -b \\ -c & a \end{pmatrix}$
（主对角换位、副对角变号），它正是伴随公式的特例。习题 2 是完整样板。

**方法 2（初等行变换法）** 一般阶数的主力：$(A\mid E)\to(E\mid A^{-1})$，
原理与步骤在第三章题型二。阶数 $\ge3$ 的数值矩阵优先用它（伴随法要算 $n^2$ 个行列式）。

**方法 3（分块求逆）** 矩阵有明显的分块三角/对角结构时按块套公式，
习题 7（下三角分块）是样板：先求对角块的逆，再由 $YA+B^{-1}C=O$ 解非对角块。

**方法 4（定义验证法）** 题目已给出候选逆或乘积关系（如已知 $AB=E$ 的变形）时，
直接按定义验算乘积为 $E$ 即可，不必从头求逆；由 $A^2+2A-3E=O$ 之类多项式关系
求逆也归此类：移项凑出 $A\cdot(\cdots)=E$，括号内即 $A^{-1}$。

## 题型三：解矩阵方程

**方法 1（乘逆法，$A$ 可逆时）** 先判型、再决定乘逆的方向：
$AX=B\Rightarrow X=A^{-1}B$（**左**乘）；$XA=B\Rightarrow X=BA^{-1}$（**右**乘）；
$AXB=C\Rightarrow X=A^{-1}CB^{-1}$。次序写反是这类题的第一大错误。

**例 2.2** 解 $AXB=C$，其中 $A=\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，
$B=\begin{pmatrix} 2 & 0 \\ 1 & 2 \end{pmatrix}$，
$C=\begin{pmatrix} 12 & 4 \\ 26 & 8 \end{pmatrix}$。
$A^{-1}=\begin{pmatrix} -2 & 1 \\ \frac32 & -\frac12 \end{pmatrix}$，
$B^{-1}=\begin{pmatrix} \frac12 & 0 \\ -\frac14 & \frac12 \end{pmatrix}$（二阶速算）。
先算 $A^{-1}C=\begin{pmatrix} -24+26 & -8+8 \\ 18-13 & 6-4 \end{pmatrix}
=\begin{pmatrix} 2 & 0 \\ 5 & 2 \end{pmatrix}$，
再右乘 $B^{-1}$ 得 $X=\begin{pmatrix} 1 & 0 \\ 2 & 1 \end{pmatrix}$。
验算：$AX=\begin{pmatrix} 5 & 2 \\ 11 & 4 \end{pmatrix}$，再乘 $B$ 得
$\begin{pmatrix} 12 & 4 \\ 26 & 8 \end{pmatrix}=C$ ✓。

**方法 2（待定元素法，$A$ 不可逆或只问部分元素时）** 设 $X$ 的元素为未知量，
按乘法定义展开成普通线性方程组，用第三章题型三的方法求解（可能无解或无穷多解，
此时矩阵方程的解也不唯一，须如实写通解）。习题 3 是方法 1 的样板
（$X=A^{-1}B$），习题 4 说明不可逆时连普通方程的消去律都会失效。

## 题型四：证明题（结构与恒等式）

**方法 1（取转置两边对照）** 涉及对称/反对称的命题，先写出候选分解
$S=\frac{A+A^T}{2}$、$T=\frac{A-A^T}{2}$ 再验算（习题 5）。
**方法 2（按定义验算乘积）** 证明某矩阵是逆、是幂的公式：直接相乘核对（习题 6（1）归纳法同理）。
**方法 3（左乘逆矩阵）** 证明由乘积等式推出因子等式：可逆时左乘 $A^{-1}$ 即消去（习题 4（1））；
不可逆时结论一般不成立，要会举反例（习题 4（2）的思路：找丢弃信息的 $A$）。
"""

EXERCISES = [
    {
        "problem_md": r"""设 $A=\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，
$B=\begin{pmatrix} 2 & 0 \\ 1 & 3 \end{pmatrix}$。计算 $AB$ 与 $BA$，并说明结论。""",
        "solution_md": r"""$$
AB=\begin{pmatrix} 1\cdot2+2\cdot1 & 1\cdot0+2\cdot3 \\ 3\cdot2+4\cdot1 & 3\cdot0+4\cdot3 \end{pmatrix}
=\begin{pmatrix} 4 & 6 \\ 10 & 12 \end{pmatrix},\quad
BA=\begin{pmatrix} 2\cdot1+0\cdot3 & 2\cdot2+0\cdot4 \\ 1\cdot1+3\cdot3 & 1\cdot2+3\cdot4 \end{pmatrix}
=\begin{pmatrix} 2 & 4 \\ 10 & 14 \end{pmatrix}.
$$

$AB\ne BA$：矩阵乘法不满足交换律。此后凡把乘积展开、移项，都必须保持左右次序。""",
    },
    {
        "problem_md": r"""求矩阵 $A=\begin{pmatrix} 1 & 2 & 3 \\ 2 & 2 & 1 \\ 3 & 4 & 3 \end{pmatrix}$ 的逆矩阵。""",
        "solution_md": r"""先算行列式：$|A|=1\cdot(2\cdot3-1\cdot4)-2\cdot(2\cdot3-1\cdot3)+3\cdot(2\cdot4-2\cdot3)=2-6+6=2\ne0$，
故 $A$ 可逆。再逐个算代数余子式：
$A_{11}=2$，$A_{12}=-3$，$A_{13}=2$；$A_{21}=6$，$A_{22}=-6$，$A_{23}=2$；
$A_{31}=-4$，$A_{32}=5$，$A_{33}=-2$。
伴随矩阵是代数余子式矩阵的**转置**：

$$
A^*=\begin{pmatrix} 2 & 6 & -4 \\ -3 & -6 & 5 \\ 2 & 2 & -2 \end{pmatrix},\qquad
A^{-1}=\frac{1}{|A|}A^*=\begin{pmatrix} 1 & 3 & -2 \\ -\frac32 & -3 & \frac52 \\ 1 & 1 & -1 \end{pmatrix}.
$$

验算第一行乘第一列：$1\cdot1+2\cdot(-\frac32)+3\cdot1=1$ ✓；
第一行乘第二列：$1\cdot3+2\cdot(-3)+3\cdot1=0$ ✓。""",
    },
    {
        "problem_md": r"""解矩阵方程 $AX=B$，其中 $A$ 为习题 2 的矩阵，
$B=\begin{pmatrix} 5 \\ 3 \\ 7 \end{pmatrix}$。""",
        "solution_md": r"""$A$ 可逆，$X=A^{-1}B$：

$$
X=\begin{pmatrix} 1 & 3 & -2 \\ -\frac32 & -3 & \frac52 \\ 1 & 1 & -1 \end{pmatrix}
\begin{pmatrix} 5 \\ 3 \\ 7 \end{pmatrix}
=\begin{pmatrix} 5+9-14 \\ -\frac{15}{2}-9+\frac{35}{2} \\ 5+3-7 \end{pmatrix}
=\begin{pmatrix} 0 \\ 1 \\ 1 \end{pmatrix}.
$$

验算：$A(0,1,1)^T=(0+2+3,\ 0+2+1,\ 0+4+3)^T=(5,3,7)^T$ ✓。""",
    },
    {
        "problem_md": r"""（1）设 $A$ 可逆，证明由 $AB=AC$ 可推出 $B=C$；
（2）举例说明 $A$ 不可逆时该结论不成立。""",
        "solution_md": r"""（1）$AB=AC$ 两边左乘 $A^{-1}$：$B=A^{-1}(AB)=A^{-1}(AC)=C$。
消去律成立靠的正是可逆性，不是 $A\ne O$。

（2）取 $A=\begin{pmatrix} 1 & 0 \\ 0 & 0 \end{pmatrix}$（$|A|=0$，不可逆），
$B=\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，
$C=\begin{pmatrix} 1 & 2 \\ 5 & 6 \end{pmatrix}$。
则 $AB=\begin{pmatrix} 1 & 2 \\ 0 & 0 \end{pmatrix}=AC$，但 $B\ne C$。
原因：$A$ 把第二行的信息整个丢弃，$B,C$ 在第二行上的差别无法从乘积中恢复——
这正是"不可逆 $=$ 信息有损失"的直观含义。""",
    },
    {
        "problem_md": r"""证明：任一方阵 $A$ 都可唯一地写成一个对称矩阵与一个反对称矩阵之和。""",
        "solution_md": r"""存在性：记 $S=\dfrac{A+A^T}{2}$，$T=\dfrac{A-A^T}{2}$，则 $A=S+T$，且
$S^T=\dfrac{A^T+A}{2}=S$ 为对称，$T^T=\dfrac{A^T-A}{2}=-T$ 为反对称。

唯一性：若 $A=S_1+T_1=S_2+T_2$（$S_i$ 对称、$T_i$ 反对称），转置得
$A^T=S_1-T_1=S_2-T_2$。两式相加：$2S_1=2S_2$；相减：$2T_1=2T_2$。故分解唯一。$\square$""",
    },
    {
        "problem_md": r"""（1）用归纳法证明 $\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}^n=\begin{pmatrix} 1 & n \\ 0 & 1 \end{pmatrix}$；
（2）求 $\begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}^n$。""",
        "solution_md": r"""（1）$n=1$ 显然。设对 $n$ 成立，则
$\begin{pmatrix} 1 & n \\ 0 & 1 \end{pmatrix}\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}
=\begin{pmatrix} 1 & n+1 \\ 0 & 1 \end{pmatrix}$，归纳完成。$\square$

（2）把矩阵拆为可交换的两部分：$\begin{pmatrix} 3 & 1 \\ 0 & 3 \end{pmatrix}=3E+N$，
$N=\begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$，且 $N^2=O$，$3E$ 与 $N$ 可交换，
可用二项式定理且只有前两项非零：

$$
(3E+N)^n=3^nE+n\cdot3^{n-1}N=\begin{pmatrix} 3^n & n\cdot3^{n-1} \\ 0 & 3^n \end{pmatrix}.
$$

要点：二项式展开对矩阵成立的**前提是两矩阵可交换**，拆成"数量矩阵 $+$ 幂零矩阵"是常用技巧。""",
    },
    {
        "problem_md": r"""用分块求逆：求 $M=\begin{pmatrix} 1 & 2 & 0 & 0 \\ 3 & 4 & 0 & 0 \\ 1 & 0 & 2 & 0 \\ 0 & 1 & 1 & 2 \end{pmatrix}$ 的逆矩阵
（提示：按前两行、前两列分块，$M=\begin{pmatrix} A & O \\ C & B \end{pmatrix}$）。""",
        "solution_md": r"""设 $A=\begin{pmatrix} 1 & 2 \\ 3 & 4 \end{pmatrix}$，
$B=\begin{pmatrix} 2 & 0 \\ 1 & 2 \end{pmatrix}$，$C=E_2$。
先算 $A^{-1}=\begin{pmatrix} -2 & 1 \\ \frac32 & -\frac12 \end{pmatrix}$（$|A|=-2$），
$B^{-1}=\begin{pmatrix} \frac12 & 0 \\ -\frac14 & \frac12 \end{pmatrix}$（$|B|=4$）。
对下三角分块矩阵，逆仍是下三角分块：设
$M^{-1}=\begin{pmatrix} A^{-1} & O \\ Y & B^{-1} \end{pmatrix}$，
由 $M^{-1}M=E$ 的左下块得 $YA+B^{-1}C=O$，即
$Y=-B^{-1}CA^{-1}=-B^{-1}A^{-1}$：

$$
B^{-1}A^{-1}=\begin{pmatrix} \frac12 & 0 \\ -\frac14 & \frac12 \end{pmatrix}
\begin{pmatrix} -2 & 1 \\ \frac32 & -\frac12 \end{pmatrix}
=\begin{pmatrix} -1 & \frac12 \\ \frac54 & -\frac12 \end{pmatrix},\quad
Y=\begin{pmatrix} 1 & -\frac12 \\ -\frac54 & \frac12 \end{pmatrix}.
$$

故 $M^{-1}=\begin{pmatrix} -2 & 1 & 0 & 0 \\ \frac32 & -\frac12 & 0 & 0 \\ 1 & -\frac12 & \frac12 & 0 \\ -\frac54 & \frac12 & -\frac14 & \frac12 \end{pmatrix}$。
验算左下块：$Y A + B^{-1}C$ 应为零：$YA$ 的 $(1,1)$ 元 $=1\cdot1+(-\frac12)\cdot3=-\frac12$，
$B^{-1}C=B^{-1}$ 的 $(1,1)$ 元 $=\frac12$，和为 $0$ ✓。""",
    },
    {
        "problem_md": r"""设 $A$ 为三阶方阵且 $|A|=2$。求：（1）$|3A|$；（2）$|A^*|$；（3）$|(2A)^{-1}|$。""",
        "solution_md": r"""（1）$|3A|=3^3|A|=27\cdot2=54$（$n$ 阶方阵数乘 $k$，行列式乘 $k^n$，不是 $k$）。

（2）由 $AA^*=|A|E$ 取行列式：$|A||A^*|=|A|^3$，$|A|\ne0$ 时约去得
$|A^*|=|A|^{n-1}=2^2=4$（此公式对 $n$ 阶方阵为 $|A^*|=|A|^{n-1}$）。

（3）$|(2A)^{-1}|=\dfrac{1}{|2A|}=\dfrac{1}{2^3\cdot2}=\dfrac{1}{16}$。""",
    },
]
