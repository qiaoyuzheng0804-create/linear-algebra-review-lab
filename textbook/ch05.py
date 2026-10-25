"""第五章 相似矩阵及二次型。"""

THEORY = r"""
## 知识点 5.1 向量的内积

**定义 5.1** 实向量 $\alpha=(a_1,\ldots,a_n)^T$、$\beta=(b_1,\ldots,b_n)^T$ 的
**内积** $[\alpha,\beta]=\sum_i a_ib_i=\alpha^T\beta$；
**长度** $\|\alpha\|=\sqrt{[\alpha,\alpha]}$；$\alpha,\beta$ 满足
$[\alpha,\beta]=0$ 称**正交**。内积满足对称性、线性性、正定性
（$[\alpha,\alpha]\ge0$ 且等号仅在 $\alpha=0$），都由定义直接验算。

**定理 5.1（Cauchy–Schwarz 不等式）** $|[\alpha,\beta]|\le\|\alpha\|\,\|\beta\|$，
等号成立当且仅当 $\alpha,\beta$ 线性相关。

证明：$\beta=0$ 时两边为零。设 $\beta\ne0$。对任意实数 $t$，
$0\le[\alpha+t\beta,\alpha+t\beta]=[\beta,\beta]t^2+2[\alpha,\beta]t+[\alpha,\alpha]$。
右端是 $t$ 的二次三项式且恒非负，故判别式 $\le0$：
$4[\alpha,\beta]^2-4[\alpha,\alpha][\beta,\beta]\le0$，即结论。
等号成立 $\iff$ 判别式为零 $\iff$ 存在 $t$ 使 $\alpha+t\beta=0$（或对称情形），
即线性相关。$\square$

由此可定义夹角 $\cos\varphi=\dfrac{[\alpha,\beta]}{\|\alpha\|\|\beta\|}$。

**Schmidt 正交化** 设 $\alpha_1,\ldots,\alpha_r$ 线性无关，令

$$
\beta_1=\alpha_1,\qquad
\beta_k=\alpha_k-\sum_{i=1}^{k-1}\frac{[\alpha_k,\beta_i]}{[\beta_i,\beta_i]}\beta_i\quad(k=2,\ldots,r),
$$

则 $\beta_1,\ldots,\beta_r$ 两两正交且与 $\alpha$ 组等价
（每步 $\beta_k$ 是 $\alpha_1,\ldots,\alpha_k$ 的组合、$\alpha_k$ 的系数为 1，
反解即得等价）。正交性由归纳验证：设 $\beta_1,\ldots,\beta_{k-1}$ 两两正交，
对 $j<k$：$[\beta_k,\beta_j]=[\alpha_k,\beta_j]-
\dfrac{[\alpha_k,\beta_j]}{[\beta_j,\beta_j]}[\beta_j,\beta_j]=0$。$\square$
再令 $\gamma_i=\beta_i/\|\beta_i\|$ 得**标准正交组**。
方阵 $Q$ 的列（行）构成标准正交组 $\iff Q^TQ=E\iff Q^{-1}=Q^T$，
这样的 $Q$ 称**正交矩阵**；它保持内积与长度：
$[Qx,Qy]=x^TQ^TQy=[x,y]$。

## 知识点 5.2 特征值与特征向量

**定义 5.2** 设 $A$ 为 $n$ 阶方阵，若存在数 $\lambda$ 与**非零**向量 $x$ 使
$Ax=\lambda x$，称 $\lambda$ 为 $A$ 的**特征值**，$x$ 为属于 $\lambda$ 的**特征向量**。

$Ax=\lambda x\iff(\lambda E-A)x=0$，有非零解 $\iff$

$$
f(\lambda)=|\lambda E-A|=0,
$$

$f$ 称为**特征多项式**，方程称**特征方程**。特征值就是特征多项式的根；
属于 $\lambda_i$ 的特征向量是齐次方程组 $(\lambda_iE-A)x=0$ 的非零解，
其解空间称**特征子空间**。展开特征多项式并与
$f(\lambda)=\prod(\lambda-\lambda_i)$ 比较系数得：

$$
\sum_{i=1}^n\lambda_i=\operatorname{tr}A\ (\text{主对角元之和}),\qquad
\prod_{i=1}^n\lambda_i=|A|.
$$

（$\lambda^{n-1}$ 的系数来自展开中取 $n-1$ 个对角因子的一项，核对即得；
常数项取 $\lambda=0$ 即 $|-A|=(-1)^n|A|$。）

**定理 5.2** 不同特征值的特征向量线性无关；更一般地，不同特征值的
特征子空间中各取若干线性无关的特征向量，合起来仍线性无关。

证明（对特征值个数归纳，关键一步）：设 $\lambda_1,\ldots,\lambda_m$ 互不相同，
$x_i$ 是属于 $\lambda_i$ 的非零特征向量。假设 $x_1,\ldots,x_m$ 是最短的
反例：$\sum x_i=0$ 且每 $x_i\ne0$。对等式作用 $A$ 得 $\sum\lambda_i x_i=0$；
再减去原式的 $\lambda_m$ 倍，得
$\sum_{i=1}^{m-1}(\lambda_i-\lambda_m)x_i=0$，
是更短的零组合且系数 $\lambda_i-\lambda_m\ne0$，由最短性每个 $x_i=0$，矛盾。
多向量版本把每个特征子空间的无关组作同样论证（先在子空间内合并成
一个向量再套用）。$\square$

## 知识点 5.3 相似矩阵与对角化

**定义 5.3** 若存在可逆 $P$ 使 $P^{-1}AP=B$，称 $B$ 是 $A$ 的**相似矩阵**，
$A$ 与 $B$ **相似**，记 $A\sim B$。

**定理 5.3** 相似矩阵有相同的特征多项式，从而相同的特征值、行列式与迹。

证明：$|\lambda E-B|=|\lambda E-P^{-1}AP|
=|P^{-1}(\lambda E-A)P|=|P^{-1}||\lambda E-A||P|=|\lambda E-A|$，
用了 $|XY|=|X||Y|$ 与 $|P^{-1}||P|=1$。$\square$
（逆命题不成立：特征多项式相同未必相似，见习题 6。）

**定理 5.4（对角化判据）** $n$ 阶方阵 $A$ 可相似对角化（即存在可逆 $P$ 使
$P^{-1}AP=\operatorname{diag}(\lambda_1,\ldots,\lambda_n)$）
当且仅当 $A$ 有 $n$ 个线性无关的特征向量。

证明：$(\Rightarrow)$ 设 $P^{-1}AP=\Lambda$，即 $AP=P\Lambda$。
把 $P$ 按列写为 $(p_1,\ldots,p_n)$，比较两边第 $i$ 列：$Ap_i=\lambda_ip_i$，
$P$ 可逆保证 $p_i\ne0$ 且 $\{p_i\}$ 线性无关，即 $n$ 个无关特征向量。
$(\Leftarrow)$ 反向构造：以 $n$ 个无关特征向量为列作 $P$，
$AP=P\Lambda$ 逐列成立，$P$ 可逆（列无关），右乘 $P^{-1}$ 即得。$\square$

**推论 1** $A$ 有 $n$ 个互不相同的特征值 $\Rightarrow$ $A$ 可对角化
（由定理 5.2）。这只是充分条件。
**推论 2（重数判据）** $A$ 可对角化 $\iff$ 每个特征值的**几何重数**
（特征子空间维数 $n-R(\lambda_iE-A)$）等于**代数重数**
（特征多项式中根的重数）。（把各特征子空间的基合起来，由定理 5.2 后半
无关，总数 $=\sum$ 几何重数；对角化 $\iff$ 总数 $=n=\sum$ 代数重数，
而几何重数恒 $\le$ 代数重数——后者：取特征子空间的基扩充为全空间基，
在该基下 $A$ 相似于分块上三角 $\begin{pmatrix} \lambda I & * \\ 0 & B_1\end{pmatrix}$，
特征多项式含因子 $(\lambda_i-\lambda)^{\text{几何重数}}\cdot|\lambda E-B_1|$，
故代数重数 $\ge$ 几何重数。）

## 知识点 5.4 实对称矩阵的对角化

**定理 5.5** 实对称矩阵的特征值全为实数；属于不同特征值的特征向量互相正交。

证明：设 $Ax=\lambda x$，$x\ne0$ 允许复向量。取共轭转置：
$\bar x^TA^T=\bar\lambda \bar x^T$。左乘：注意 $A^T=A$，
一方面 $\bar x^TAx=\lambda \bar x^Tx$，另一方面
$\bar x^TAx=\bar x^TA^Tx=(A\bar x)^Tx$，而 $A\bar x=\bar\lambda \bar x$，
故 $=\bar\lambda \bar x^Tx$。$\bar x^Tx=\sum|x_i|^2>0$，约去得
$\lambda=\bar\lambda$，即 $\lambda$ 为实数。
正交性：$\lambda_1\ne\lambda_2$ 时，
$\lambda_1[x_1,x_2]=[Ax_1,x_2]=[x_1,Ax_2]=\lambda_2[x_1,x_2]$，
故 $(\lambda_1-\lambda_2)[x_1,x_2]=0$，得 $[x_1,x_2]=0$。$\square$

**定理 5.6（谱定理）** 设 $A$ 为 $n$ 阶实对称矩阵，则存在正交矩阵 $Q$ 使

$$
Q^{-1}AQ=Q^TAQ=\operatorname{diag}(\lambda_1,\ldots,\lambda_n).
$$

证明（归纳法）：$n=1$ 平凡。设对 $n-1$ 阶成立。
由定理 5.5，$A$ 有实特征值 $\lambda_1$ 与实单位特征向量 $q_1$。
把 $q_1$ 扩充为 $\mathbb R^n$ 的标准正交基 $q_1,\ldots,q_n$
（先扩成基再 Schmidt 化），记 $Q_1=(q_1,\ldots,q_n)$ 为正交矩阵。
$A q_1=\lambda_1q_1$ 给出 $Q_1^{-1}AQ_1$ 的第一列为 $(\lambda_1,0,\ldots,0)^T$；
该矩阵仍对称（$(Q_1^TAQ_1)^T=Q_1^TA^TQ_1$），故第一行也是
$(\lambda_1,0,\ldots,0)$，于是
$Q_1^TAQ_1=\begin{pmatrix} \lambda_1 & 0 \\ 0 & A_1 \end{pmatrix}$，
$A_1$ 为 $n-1$ 阶实对称矩阵。由归纳假设存在正交 $Q_2'$ 使
$Q_2'^TA_1Q_2'$ 对角化，取 $Q=Q_1\begin{pmatrix} 1 & 0 \\ 0 & Q_2'\end{pmatrix}$，
$Q$ 是正交矩阵之积、仍正交，且 $Q^TAQ$ 为对角阵。$\square$

## 知识点 5.5 二次型及其标准形

**定义 5.4** $n$ 元二次齐次多项式
$f(x)=\sum_i a_{ii}x_i^2+2\sum_{i<j}a_{ij}x_ix_j$ 称为**二次型**，
可写为 $f=x^TAx$，其中 $A$ 为实对称矩阵（交叉项系数对半分给 $a_{ij},a_{ji}$）。
$A$ 称为二次型的矩阵，$A$ 的秩称为二次型的秩。

可逆线性变换 $x=Cy$（$C$ 可逆）代入得 $f=y^T(C^TAC)y$，
$C^TAC$ 仍对称。$B=C^TAC$ 称为 $A$ 与 $B$ **合同**。
合同是等价关系；合同不改变秩，但一般改变特征值（区别于相似）。

**定义 5.5** 只含平方项的二次型称**标准形**；系数仅为 $1,-1,0$ 的称**规范形**。

**定理 5.7** 任一二次型都可经可逆线性变换化为标准形。

证明思路（配方法，归纳于 $n$）：若某个 $a_{ii}\ne0$，把含 $x_i$ 的项
凑成完全平方，剩余是 $n-1$ 元二次型，归纳处理；
若 $a_{ii}$ 全为零但某 $a_{ij}\ne0$，先作变换 $x_i=y_i+y_j$，$x_j=y_i-y_j$
制造出平方项再回到前一情形。细节在方法三中以实例演示。$\square$
对实对称矩阵，定理 5.6 给出更强的**正交变换法**：用正交变换化标准形，
系数恰是 $A$ 的全部特征值。

**定理 5.8（惯性定理）** 实二次型的标准形中，正平方项个数 $p$ 与负平方项个数 $q$
由二次型唯一确定（不依赖所用变换），分别称**正、负惯性指数**；
规范形为 $y_1^2+\cdots+y_p^2-y_{p+1}^2-\cdots-y_{p+q}^2$。

证明要点：设两组可逆变换把 $f$ 化为以 $z$ 与 $w$ 为变量的两个标准形，
正项个数分别为 $p_1,p_2$，反设 $p_1>p_2$。
取子空间 $L_1=\{z_1=\cdots\}$：令第一个标准形中正项变量自由、其余为零的
那批方向构成 $p_1$ 维子空间 $U$（在 $U$ 上 $f>0$ 除原点）；
第二个标准形中非负项（正项与零项）方向构成 $n-p_2$ 维子空间 $W$
（在 $W$ 上 $f\ge0$）。$p_1+(n-p_2)>n$ 蕴含 $U\cap W$ 含非零向量
（维数公式），该向量使 $f$ 同时 $>0$ 与 $\le0$，矛盾。故 $p_1\le p_2$；
对称地 $p_2\le p_1$。$\square$

**定义 5.6（正定）** 若对一切 $x\ne0$ 有 $f(x)=x^TAx>0$，
称 $f$ **正定**，$A$ 称**正定矩阵**。

**定理 5.9** 设 $A$ 为实对称矩阵，则下列命题等价：
（1）$A$ 正定；（2）$A$ 的特征值全为正；（3）正惯性指数 $=n$；
（4）$A$ 与单位矩阵合同；（5）存在可逆 $P$ 使 $A=P^TP$；
（6）$A$ 的各阶**顺序主子式**全大于零。

证明要点：（1）$\iff$（3）：由惯性定理，标准形在每个坐标方向取值
即系数符号，$f$ 恒正当且仅当无非正系数。（2）$\iff$（3）：正交变换法的
标准形系数即特征值。（3）$\iff$（4）：规范形即 $E$ 对应的二次型 $\sum y_i^2$。
（4）$\Rightarrow$（5）：$A=C^TEC=(C)^T(C)$ 形如 $P^TP$（取 $P=C$）；
（5）$\Rightarrow$（1）：$x^TP^TPx=\|Px\|^2>0$（$Px\ne0$ 由 $P$ 可逆）。
（1）$\Rightarrow$（6）：前 $k$ 个变量构成的子二次型仍正定
（令其余变量为 $0$ 即得），其矩阵是 $A$ 的 $k$ 阶顺序主子矩阵 $A_k$；
$A_k$ 正定则由已证的（1）$\iff$（2）特征值全正，
故顺序主子式 $|A_k|=\prod\lambda>0$。
（6）$\Rightarrow$（3）：对 $n$ 归纳。$a_{11}=\Delta_1>0$。
用合同变换（第一行/列的倍加，即对 $x_1$ 配方）把 $A$ 化为
$\operatorname{diag}(a_{11},B)$，其中 $B=A_1-\frac{1}{a_{11}}\alpha\alpha^T$
是 Schur 补（$A=\begin{pmatrix} a_{11} & \alpha^T \\ \alpha & A_1 \end{pmatrix}$）。
$B$ 的 $k$ 阶顺序主子式满足 $\det B_k=\dfrac{\Delta_{k+1}}{a_{11}}$：
对 $A$ 的 $(k+1)$ 阶顺序主子矩阵做同样的消元，行列式不变，
得 $a_{11}\det B_k=\Delta_{k+1}$。$\Delta_{k+1}>0$ 且 $a_{11}>0$
给出 $B$ 的各阶顺序主子式全正，由归纳假设 $B$ 正定（惯性指数 $n-1$）。
合同变换不改变惯性指数，而 $f$ 在该变换下 $=a_{11}y_1^2+g(y_2,\ldots,y_n)$
且 $a_{11}>0$，故 $A$ 的正惯性指数 $=1+(n-1)=n$。$\square$

---

# 解题方法

## 方法一：求特征值与特征向量的标准步骤

（1）解特征方程 $|\lambda E-A|=0$ 得全部特征值（注意重数）；
（2）对每个 $\lambda_i$，解齐次方程组 $(\lambda_iE-A)x=0$，
其基础解系即该特征子空间的基，全部特征向量 $=k\cdot$（基础解系组合），$k\ne0$。
**小技巧**：二、三阶先用 $\sum\lambda_i=\operatorname{tr}A$、$\prod\lambda_i=|A|$ 验核。

**例 5.1** 求 $A=\begin{pmatrix} 4 & 2 \\ 1 & 3 \end{pmatrix}$ 的特征值与特征向量。
$|\lambda E-A|=(\lambda-4)(\lambda-3)-2=\lambda^2-7\lambda+10=(\lambda-5)(\lambda-2)$，
$\lambda_1=5$，$\lambda_2=2$（验核：和 $7=\operatorname{tr}A$ ✓，积 $10=|A|$ ✓）。
$\lambda=5$：$(5E-A)=\begin{pmatrix} 1 & -2 \\ -1 & 2 \end{pmatrix}$，
方程 $x_1=2x_2$，基础解 $(2,1)^T$，特征向量 $k(2,1)^T\ (k\ne0)$。
$\lambda=2$：$(2E-A)=\begin{pmatrix} -2 & -2 \\ -1 & -1 \end{pmatrix}$，
$x_1=-x_2$，特征向量 $k(1,-1)^T\ (k\ne0)$。

## 方法二：判别对角化与求 $P$

（1）求全部特征值；（2）对每个特征值算几何重数 $n-R(\lambda_iE-A)$；
（3）几何重数处处等于代数重数 $\iff$ 可对角化，
各特征子空间的基按列拼成 $P$，$P^{-1}AP=\operatorname{diag}(\lambda_i)$（顺序对应）；
（4）求 $A$ 的高次幂等问题时用 $A^k=P\Lambda^kP^{-1}$。

**例 5.2** $A=\begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$。
特征方程 $\lambda^2-2\lambda-3=0$，$\lambda=3,-1$ 互不相同，由推论 1 可对角化。
$\lambda=3$：特征向量 $(1,1)^T$；$\lambda=-1$：特征向量 $(1,-1)^T$。
取 $P=\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，
则 $P^{-1}AP=\operatorname{diag}(3,-1)$（读者可用 $AP=P\Lambda$ 逐列验算代替求逆验算）。

## 方法三：化二次型为标准形

**正交变换法**（要求变换为正交变换、或要读特征值信息时用）：
写出 $A$ $\to$ 求特征值与两两正交的单位特征向量（重根时在特征子空间内
Schmidt 正交化）$\to$ 拼正交阵 $Q$，$x=Qy$ 下 $f=\sum\lambda_iy_i^2$。
**配方法**（只要求标准形时通常更快）：逐个变量配方；
无平方项时先用和差变换造平方项。
**初等变换法**：对 $(A\mid E)$ 同步做行、列成对变换把 $A$ 化对角，
下方即变换矩阵 $C$。

**例 5.3（正交变换法）** 化 $f=2x_1^2+5x_2^2+5x_3^2+4x_1x_2-4x_1x_3-8x_2x_3$ 为标准形。
$A=\begin{pmatrix} 2 & 2 & -2 \\ 2 & 5 & -4 \\ -2 & -4 & 5 \end{pmatrix}$，
特征多项式可算得 $-(\lambda-1)^2(\lambda-10)$（验核：迹 $12=1+1+10$ ✓）。
$\lambda=10$：解 $(10E-A)x=0$ 得特征向量 $(1,2,-2)^T$，单位化 $q_3=\frac13(1,2,-2)^T$。
$\lambda=1$：$(E-A)x=0$ 化为 $x_1+2x_2-2x_3=0$，取基础解
$(-2,1,0)^T$、$(2,0,1)^T$；Schmidt 正交化：
$\beta_1=(-2,1,0)^T$；
$\beta_2=(2,0,1)^T-\frac{-4}{5}(-2,1,0)^T=(\frac25,\frac45,1)^T$；
单位化 $q_1=\frac{1}{\sqrt5}(-2,1,0)^T$，$q_2=\frac{1}{3\sqrt5}(2,4,5)^T$。
取正交阵 $Q=(q_1,q_2,q_3)$，正交变换 $x=Qy$ 下

$$
f=y_1^2+y_2^2+10y_3^2.
$$

**例 5.4（配方法，无平方项）** 化 $f=x_1x_2+x_2x_3+x_3x_1$ 为标准形。
无平方项，先作可逆变换 $x_1=y_1+y_2$，$x_2=y_1-y_2$，$x_3=y_3$
（行列式为 $-2\ne0$，可逆）：
$x_1x_2=y_1^2-y_2^2$，$x_2x_3+x_3x_1=2y_1y_3$，故
$f=y_1^2+2y_1y_3-y_2^2=(y_1+y_3)^2-y_2^2-y_3^2$。
再令 $z_1=y_1+y_3$，$z_2=y_2$，$z_3=y_3$（仍可逆），标准形
$f=z_1^2-z_2^2-z_3^2$，正惯性指数 $1$、负惯性指数 $2$、秩 $3$。

## 方法四：判定正定性

三条等价路线（定理 5.9）：求全部特征值看符号（阶低时快）、
看顺序主子式（数字矩阵最常用）、化标准形看惯性指数。
负定把不等号反向（顺序主子式奇数阶 $<0$、偶数阶 $>0$）；
只要某主子式或某特征值给出非正、非负的反例方向，立即排除正定。
"""

EXERCISES = [
    {
        "problem_md": r"""把 $\alpha_1=(1,1,0)^T$，$\alpha_2=(1,0,1)^T$，$\alpha_3=(0,1,1)^T$ 化为标准正交组。""",
        "solution_md": r"""先 Schmidt 正交化。$\beta_1=\alpha_1=(1,1,0)^T$。
$\beta_2=\alpha_2-\dfrac{[\alpha_2,\beta_1]}{[\beta_1,\beta_1]}\beta_1
=(1,0,1)^T-\frac12(1,1,0)^T=(\frac12,-\frac12,1)^T$。
$\beta_3=\alpha_3-\dfrac{[\alpha_3,\beta_1]}{2}\beta_1-
\dfrac{[\alpha_3,\beta_2]}{[\beta_2,\beta_2]}\beta_2$，
其中 $[\alpha_3,\beta_1]=1$，$[\alpha_3,\beta_2]=\frac12$，
$[\beta_2,\beta_2]=\frac32$，系数 $\frac13$：
$\beta_3=(0,1,1)^T-\frac12(1,1,0)^T-\frac13(\frac12,-\frac12,1)^T
=(-\frac23,\frac23,\frac23)^T$。
长度：$\|\beta_1\|=\sqrt2$，$\|\beta_2\|=\sqrt{\frac32}$，$\|\beta_3\|=\frac{2}{\sqrt3}$。
单位化得标准正交组

$$
\gamma_1=\frac{1}{\sqrt2}(1,1,0)^T,\quad
\gamma_2=\frac{1}{\sqrt6}(1,-1,2)^T,\quad
\gamma_3=\frac{1}{\sqrt3}(-1,1,1)^T.
$$

（$\gamma_2$ 由 $(\frac12,-\frac12,1)$ 同向单位化：乘 $\frac{2}{\sqrt6}$；
验算 $[\gamma_2,\gamma_3]=\frac{1}{\sqrt{18}}(-1-1+2)=0$ ✓。）""",
    },
    {
        "problem_md": r"""求 $A=\begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$ 的特征值和全部特征向量，
并求可逆矩阵 $P$ 使 $P^{-1}AP$ 为对角阵。""",
        "solution_md": r"""$|\lambda E-A|=(\lambda-3)^2-1=(\lambda-4)(\lambda-2)$，$\lambda_1=4$，$\lambda_2=2$。
$\lambda=4$：$(4E-A)=\begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$，
特征向量 $k(1,1)^T$；$\lambda=2$：特征向量 $k(1,-1)^T$（$k\ne0$）。
两特征值不同故可对角化，取 $P=\begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix}$，
$P^{-1}AP=\operatorname{diag}(4,2)$。（$A$ 对称，$P$ 的两列还互相正交，
再单位化即得正交对角化，与定理 5.6 一致。）""",
    },
    {
        "problem_md": r"""设三阶方阵 $A$ 的特征值为 $1,-1,2$，$B=A^3-2A$。求 $B$ 的特征值与 $|B|$。""",
        "solution_md": r"""若 $Ax=\lambda x$，则 $A^3x=\lambda^3x$、$Bx=(\lambda^3-2\lambda)x$，
故 $B$ 的特征值（在同一特征向量上取得）为 $\lambda^3-2\lambda$：
$\lambda=1\to-1$；$\lambda=-1\to1$；$\lambda=2\to4$。
$|B|=(-1)\cdot1\cdot4=-4$（行列式等于特征值之积——三特征值互异，
$A$ 可对角化，$B$ 与 $\operatorname{diag}(-1,1,4)$ 相似，
或者用 $|B|=|A^3-2A|$ 配合 $A$ 的特征多项式论证，此处由相似对角化最直接）。""",
    },
    {
        "problem_md": r"""设 $A=\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$。求 $A^{10}$。""",
        "solution_md": r"""特征值：$|\lambda E-A|=(\lambda-1)(\lambda-3)$，$\lambda=1,3$；
特征向量依次 $(1,-1)^T$、$(1,1)^T$。
取 $P=\begin{pmatrix} 1 & 1 \\ -1 & 1 \end{pmatrix}$，
$P^{-1}=\frac12\begin{pmatrix} 1 & -1 \\ 1 & 1 \end{pmatrix}$，
$A=P\operatorname{diag}(1,3)P^{-1}$，于是

$$
A^{10}=P\begin{pmatrix} 1 & 0 \\ 0 & 3^{10} \end{pmatrix}P^{-1}
=\frac12\begin{pmatrix} 1+3^{10} & -1+3^{10} \\ -1+3^{10} & 1+3^{10} \end{pmatrix}
=\begin{pmatrix} 29525 & 29524 \\ 29524 & 29525 \end{pmatrix}
$$

（$3^{10}=59049$）。验算迹：$2\cdot29525=59050=1+3^{10}=\operatorname{tr}(A^{10})$ ✓。""",
    },
    {
        "problem_md": r"""判别二次型 $f=x_1^2+2x_2^2+3x_3^2+4x_1x_2+2x_2x_3$ 的正定性。""",
        "solution_md": r"""二次型矩阵 $A=\begin{pmatrix} 1 & 2 & 0 \\ 2 & 2 & 1 \\ 0 & 1 & 3 \end{pmatrix}$（注意交叉项系数对半分）。
顺序主子式：$\Delta_1=1>0$；$\Delta_2=\begin{vmatrix} 1 & 2 \\ 2 & 2 \end{vmatrix}=-2<0$。
由定理 5.9（6），顺序主子式未全大于零，$f$ 不是正定的。
实际上配方 $f=(x_1+2x_2)^2-2x_2^2+\cdots$ 显示标准形含负系数，
$f$ 是不定的（取 $x=(−2,1,0)$ 得 $f=-2<0$，取 $x=(1,0,0)$ 得 $f=1>0$）。""",
    },
    {
        "problem_md": r"""判断正误并说明理由：若三阶方阵 $A,B$ 有相同的特征多项式，则 $A$ 与 $B$ 相似。""",
        "solution_md": r"""**错误。** 反例：$A=\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}$ 与 $E=\begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$（二阶即可）。
二者特征多项式同为 $(\lambda-1)^2$。但 $E$ 与任何矩阵相似变换下仍是 $E$
（$P^{-1}EP=E$），而 $A\ne E$，故 $A$ 与 $E$ 不相似。
从判据看：$\lambda=1$ 对 $A$ 的几何重数 $=2-R(E-A)=2-1=1<2=$ 代数重数，
$A$ 不可对角化；特征多项式看不见几何重数，这正是定理 5.3 逆命题失败的原因。""",
    },
    {
        "problem_md": r"""用配方法化二次型 $f=2x_1^2+2x_2^2+2x_1x_2+2x_1x_3$ 为标准形，并写出所作可逆线性变换与正、负惯性指数。""",
        "solution_md": r"""含 $x_1$ 的项：$2x_1^2+2x_1x_2+2x_1x_3$。凑平方：
$2\big(x_1+\frac12x_2+\frac12x_3\big)^2$ 展开为
$2x_1^2+2x_1x_2+2x_1x_3+\frac12x_2^2+\frac12x_3^2+x_2x_3$，
故 $f=2\big(x_1+\frac{x_2+x_3}{2}\big)^2+
\big(2-\frac12\big)x_2^2-\frac12x_3^2-x_2x_3$
$=2(\cdots)^2+\frac32x_2^2-x_2x_3-\frac12x_3^2$。
对余下的 $x_2,x_3$ 部分继续配方：
$\frac32\big(x_2-\frac13x_3\big)^2=\frac32x_2^2-x_2x_3+\frac16x_3^2$，
余 $-\frac12-\frac16=-\frac23x_3^2$。于是

$$
f=2y_1^2+\frac32y_2^2-\frac23y_3^2,
$$

其中变换为 $y_1=x_1+\frac12x_2+\frac12x_3$，$y_2=x_2-\frac13x_3$，$y_3=x_3$
（上三角、主对角元为 $1,1,1$，行列式 $1\ne0$，可逆）。
正惯性指数 $2$，负惯性指数 $1$，秩 $3$。""",
    },
    {
        "problem_md": r"""证明：实对称矩阵 $A$ 正定的当且仅当存在可逆矩阵 $P$ 使 $A=P^TP$。""",
        "solution_md": r"""（$\Leftarrow$）若 $A=P^TP$，对 $x\ne0$：$Px\ne0$（$P$ 可逆），
$x^TAx=x^TP^TPx=\|Px\|^2>0$，故 $A$ 正定。

（$\Rightarrow$）$A$ 正定则由定理 5.6，存在正交阵 $Q$ 使
$A=Q\operatorname{diag}(\lambda_1,\ldots,\lambda_n)Q^T$，且由定理 5.9
特征值全正。令 $D=\operatorname{diag}(\sqrt{\lambda_1},\ldots,\sqrt{\lambda_n})$，
则 $A=(QD)(QD)^T$。取 $P=(QD)^T$（可逆：$Q$ 正交、$D$ 对角元非零），
$P^TP=(QD)(QD)^T=A$。$\square$

这条等价性说明正定矩阵恰是"某个可逆变换的 Gram 矩阵"，也是数值计算中
Cholesky 分解 $A=LL^T$ 存在的前提。""",
    },
]
