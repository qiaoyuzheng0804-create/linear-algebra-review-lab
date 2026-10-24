"""第 8 章：厄米性、酉性与谱定理，连接有限维量子测量。"""
THEORY = r"""
## 动机与直觉

一般特征向量基可能不正交，甚至不够组成基。物理上却希望能量等测量值是实数，
互斥结果对应正交子空间，而且换到测量基后仍可用模平方计算概率。
厄米条件同时解决这三个问题。本章从内积的共轭规则出发证明这些性质，
再用正交补归纳构造谱分解；酉矩阵则描述不改变内积的坐标变换与封闭系统演化。
这里的谱定理是有限维定理，证明中每次正交补都会使维数下降一。

## 定义

沿用第 6 章约定 $\langle x,y\rangle=x^\dagger y$，第一变量共轭线性。
**定义 8.1（伴随与厄米）**：$A^\dagger=\bar A^T$，满足
$\langle x,Ay\rangle=\langle A^\dagger x,y\rangle$。
方阵 $H=H^\dagger$ 称厄米矩阵；实矩阵时就是实对称矩阵。
**定义 8.2（酉与正规）**：方阵 $U^\dagger U=I$ 称酉矩阵，方阵
$A^\dagger A=AA^\dagger$ 称正规矩阵。厄米与酉都是正规矩阵。
正交归一列给 $U^\dagger U=I$；它使 $Ux=0$ 只能有 $x=0$，
方阵的秩因而为 n，可逆且 $U^{-1}=U^\dagger$，故也有 $UU^\dagger=I$。
**定义 8.3（谱投影）**：若 $q_j$ 是本征空间 $E_\lambda$ 的正交归一基，定义
$P_\lambda=\sum_{j:\lambda_j=\lambda}q_jq_j^\dagger$。
$P_\lambda x$ 保留 x 在整个 $E_\lambda$ 中的分量；重根时不能只选其中一列冒充整个投影。

## 核心推导

### 推导 8.1：厄米矩阵的特征值为什么为实数

动机：测量读数必须是实数；要从矩阵条目上的共轭对称推到本征值上的实性。
设 $Hx=\lambda x$ 且 $x\ne0$，暂时允许 $\lambda$ 是复数。令 $s=x^\dagger Hx$，
这是一个 $1\times1$ 标量。按乘积取伴随时逆序的规则，
$$ \bar s=s^\dagger=(x^\dagger Hx)^\dagger=x^\dagger H^\dagger x=x^\dagger Hx=s. $$
最后一个等号恰好使用 $H^\dagger=H$，所以 s 为实数。
另一方面，把本征方程左乘 $x^\dagger$ 得
$s=\lambda x^\dagger x=\lambda\|x\|^2$。
两边取共轭，因为 $\|x\|^2$ 为实数，得 $\bar s=\bar\lambda\|x\|^2$。
两式相减并用 $s=\bar s$，
$$ (\lambda-\bar\lambda)\|x\|^2=0. $$
正定性与 $x\ne0$ 给 $\|x\|^2>0$，可以除去这一因子，故 $\lambda=\bar\lambda\in\mathbb R$。
这里没有先假定特征值实数，也没有把 $Hx$ 的共轭误写成 $H\bar x$。

条件与失效：只要求 H 厄米而不要求正定，因此实特征值可以为负或零。
若 $H^\dagger\ne H$，第一行中 $\bar s=s$ 的等号没有依据；
例如 $H=(i)$ 有特征值 i。若把非零本征向量误换成零向量，最后除以范数的步骤也会失效。

### 推导 8.2：不同本征值对应的方向正交

动机：不同测量读数应可通过互相垂直的子空间区分，先证明两个本征方向的内积为零。
设 $Hx=\lambda x,Hy=\mu y$，由推导 8.1，$\lambda,\mu$ 都实。
先从右侧本征方程计算 $s=x^\dagger Hy=\mu x^\dagger y$。
再对另一个标量 $y^\dagger Hx=\lambda y^\dagger x$ 取共轭：
$$ \overline{y^\dagger Hx}=x^\dagger H^\dagger y=x^\dagger Hy,
\qquad \overline{\lambda y^\dagger x}=\bar\lambda x^\dagger y=\lambda x^\dagger y. $$
因此同一个 s 既等于 $\mu x^\dagger y$ 又等于 $\lambda x^\dagger y$，相减得
$(\mu-\lambda)x^\dagger y=0$。当 $\lambda\ne\mu$，可除以差值得 $x^\dagger y=0$。
当两本征值相等，等式只剩 $0=0$，不能推出正交；需在该本征空间内部使用 Gram–Schmidt。
在这个子空间内任意线性组合仍满足 $Hz=\lambda z$，所以正交化不会破坏本征性质。

条件与失效：两个关键等号分别需要厄米性和本征值实性；互异性用在最后的除法。
重根不会破坏谱定理，只是不能靠本结论直接证明同一本征空间里任选向量正交。
非厄米但可对角化的矩阵也可能有不正交本征基，如 $\begin{pmatrix}1&1\\0&2\end{pmatrix}$
的本征向量 $(1,0)^T,(1,1)^T$ 内积为 1。

### 推导 8.3：厄米谱定理的完整归纳证明

动机：两个不同本征向量正交，还没有证明本征向量足够多；归纳法要制造一套完整的正交本征基。
在复空间 $\mathbb C^n$ 上对维数 n 归纳。n=1 时矩阵是实标量，单位向量即所需基；
n=0 时空基也符合结论。设维数小于 n 的厄米算符都有正交归一本征基。

第一步是保证至少一个本征向量存在。第 4、5 章由行列式定义得到
$p_H(z)=\det(zI-H)$ 是首一 n 次多项式。这里调用**代数学基本定理**：
每个非常数复系数多项式在复数域至少有一个根。它是多项式理论的前置定理，
不是由谱定理反推出来的。取根 $\lambda$，则 $\det(\lambda I-H)=0$，
第 3、4 章的不可逆与非零核等价关系给非零 x 使 $Hx=\lambda x$。
除以正范数取 $q_1=x/\|x\|$，推导 8.1 保证 $\lambda$ 为实数。

第二步构造剩余空间 $W=\{w:q_1^\dagger w=0\}$。
任意 z 可以写成 $z=q_1(q_1^\dagger z)+(z-q_1q_1^\dagger z)$，后项与 $q_1$ 内积为零；
而 $cq_1\in W$ 给 $c=0$。故空间是 $\operatorname{span}(q_1)\oplus W$ 的正交直和，
$\dim W=n-1$。对 $w\in W$，逐项计算
$$ q_1^\dagger Hw=(Hq_1)^\dagger w=(\lambda q_1)^\dagger w
=\bar\lambda q_1^\dagger w=0. $$
第一步使用 H 厄米，末步使用 w 的定义，因而 $Hw\in W$：W 在 H 下不变。

第三步验证限制算符 $H_W:W\to W$ 仍厄米。先从 W 的任意基用第 6 章正交化得到
正交归一基 $w_2,\ldots,w_n$，其矩阵 $B_{ij}=w_i^\dagger Hw_j$。
计算共轭转置条目：
$\overline{B_{ji}}=\overline{w_j^\dagger Hw_i}=w_i^\dagger H^\dagger w_j=w_i^\dagger Hw_j=B_{ij}$。
因此 $B^\dagger=B$。等价地，W 中任意 u、v 满足
$\langle u,H_Wv\rangle=\langle H_Wu,v\rangle$，两端确实都在 W 内。

第四步对维数 n−1 的 $H_W$ 应用归纳假设，得到 W 中正交归一本征基 $q_2,\ldots,q_n$。
它们作为全空间的向量仍满足 $Hq_j=\lambda_jq_j$，因为限制只缩小定义域，没有更改 H 的作用。
将 $q_1$ 拼在前面，$U=(q_1,\ldots,q_n)$ 的列两两正交归一且张满全空间，
所以 $U^\dagger U=I$。逐列 $HU=U\Lambda$，左乘 $U^\dagger$，得
$$ U^\dagger HU=\Lambda,\qquad H=U\Lambda U^\dagger,\quad \Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_n)\in\mathbb R^{n\times n}.\tag{8.1} $$
这完成存在性、完整性和酉性的归纳，未把“有 n 个根”误当成“有 n 个独立本征向量”。
若 H 实对称，已知 $\lambda$ 实，将复本征向量写成 a+ib，实部和虚部分别满足本征方程，
至少一个非零，故首向量可选实；之后在实正交补中重复，得到实正交对角化。

动机：把同一测量值的多个方向合并，避免让概率依赖于简并空间内部的任意选基。
由矩阵乘法逐列展开 (8.1)，$H=\sum_j\lambda_jq_jq_j^\dagger$。
合并相同值，并计算 $q_i^\dagger q_j=\delta_{ij}$，得到
$$ H=\sum_\lambda\lambda P_\lambda,\quad
P_\lambda P_\mu=\delta_{\lambda\mu}P_\lambda,\quad
\sum_\lambda P_\lambda=I.\tag{8.2} $$
因为 $(q_iq_i^\dagger)(q_jq_j^\dagger)=q_i(q_i^\dagger q_j)q_j^\dagger$，
不等指标项为零，相等指标项保留；这也证明每个 P 自伴且幂等。
对任意 x，$P_\lambda x\in E_\lambda$，残差与该空间各 q 内积为零，
由第 6 章正交投影唯一性可知它与本征空间内部的选基无关。

条件与失效：有限维、正定内积和厄米性支撑上述每一步。
对非正规矩阵，虽然特征多项式仍有根，但 $q_1^\dagger Hw=(Hq_1)^\dagger w$ 不能成立，
正交补可能不再不变。例如 $J=\begin{pmatrix}1&1\\0&1\end{pmatrix}$ 的本征向量是 $e_1$，
但 $Je_2=e_1+e_2\notin e_1^\perp$，归纳的第二步就在此断掉。

### 推导 8.4：酉性、保内积与保范数的逐向等价

动机：实验可直接检查长度或概率归一化，代数却用 $U^\dagger U=I$；需要证明这些检验等价。
对复线性方阵 U，若 $U^\dagger U=I$，则
$\langle Ux,Uy\rangle=x^\dagger U^\dagger Uy=x^\dagger y=\langle x,y\rangle$。
若保内积，令 y=x 得 $\|Ux\|^2=\|x\|^2$，两边取非负平方根即保范数。

反向从保范数恢复内积，先推导所需的极化恒等式。令 $z=\langle x,y\rangle$。
分配律和共轭对称给
$$ \|x+y\|^2=\|x\|^2+\|y\|^2+z+\bar z,\quad
\|x-y\|^2=\|x\|^2+\|y\|^2-z-\bar z. $$
相减为 $4\operatorname{Re}z$。又因为第一变量共轭线性，
$$ \|x+iy\|^2=\|x\|^2+\|y\|^2+iz-i\bar z,
\quad\|x-iy\|^2=\|x\|^2+\|y\|^2-iz+i\bar z. $$
后两式相减为 $2i(z-\bar z)=-4\operatorname{Im}z$，故
$$ \langle x,y\rangle=\frac14\left(\|x+y\|^2-\|x-y\|^2-i\|x+iy\|^2+i\|x-iy\|^2\right).\tag{8.3} $$
将 x、y 换成 Ux、Uy；线性性给 $Ux\pm Uy=U(x\pm y)$ 及
$Ux\pm iUy=U(x\pm iy)$，再逐项用保范数，右边恢复原 x、y 的式子，所以保内积。
最后保内积意味着对任意 x、y，$x^\dagger(U^\dagger U-I)y=0$。
分别取 $x=e_i,y=e_j$，得到矩阵每个条目为零，故 $U^\dagger U=I$，闭合等价链：
$$ U^\dagger U=I\iff\langle Ux,Uy\rangle=\langle x,y\rangle\ (\forall x,y)
\iff\|Ux\|=\|x\|\ (\forall x).\tag{8.4} $$

条件与失效：必须是复线性映射。共轭映射 $Kx=\bar x$ 保范数，但
$K(iy)=-iKy$，极化证明中把 i 移进映射的一步失效，它一般把内积变成其共轭。
矩形矩阵也可保内积，但只能推出列等距 $U^\dagger U=I$，不能推出双侧可逆；
例如从二维嵌入三维的等距矩阵并不覆盖第三个方向。
只检查每列长度为一还不够：两列相同的矩阵能保持每个基向量长度，却不能保持它们差的长度。

### 推导 8.5：正规矩阵的 Schur 路线

动机：厄米性保证实谱，但酉矩阵可以有复相位谱；找出它们共享的“存在正交本征基”条件。
先证明复方阵的 Schur 三角化路线。用代数学基本定理取单位本征向量 q，
将它正交扩张成基 $Q=(q,Q_2)$。因为 $Aq=\lambda q$，
$$ Q^\dagger AQ=\begin{pmatrix}\lambda&b\\0&B\end{pmatrix}. $$
左下块为零来自 $Q_2^\dagger Aq=\lambda Q_2^\dagger q=0$，这里不需要正交补不变。
对任意的 (n−1) 阶 B 归纳三角化为 $R^\dagger BR=T_2$，
再以 $Q\operatorname{diag}(1,R)$ 换基得到上三角 T。n=1 起步完成归纳。
这一步允许右上块 b 非零，因此对缺陷矩阵仍成立。

若 A 正规，$T=Q^\dagger AQ$ 满足
$T^\dagger T=Q^\dagger A^\dagger AQ=Q^\dagger AA^\dagger Q=TT^\dagger$。
比较 (1,1) 条目：上三角使第一列只有 $t_{11}$，故
$ (T^\dagger T)_{11}=|t_{11}|^2$；而
$ (TT^\dagger)_{11}=\sum_j|t_{1j}|^2$。
相等迫使 $\sum_{j>1}|t_{1j}|^2=0$。每项非负，所以每个 $t_{1j}=0$。
此时 $T=\operatorname{diag}(t_{11},T_2)$，正规等式的右下块给 $T_2$ 仍正规。
对 T₂ 重复这一条目比较，逐行消去全部超对角项，最终 T 是对角阵。
于是正规 A 可酉对角化。反向若 $A=U\Lambda U^\dagger$，则
$A^\dagger A=U\Lambda^\dagger\Lambda U^\dagger$、
$AA^\dagger=U\Lambda\Lambda^\dagger U^\dagger$，对角标量乘法交换，两者相等。
这证明复方阵正规当且仅当可酉对角化。

条件与失效：非正规矩阵仍能 Schur 三角化，失败的是比较行列平方和时没有相等依据。
上面的 J 满足 $JJ^\dagger=\begin{pmatrix}2&1\\1&1\end{pmatrix}$，
$J^\dagger J=\begin{pmatrix}1&1\\1&2\end{pmatrix}$，第一条目不同，所以不能迫使超对角线的 1 消失。
正规矩阵不必厄米；例如 $(i)$ 正规但谱非实。若 U 酉且 $Uv=\lambda v$，
保范数给 $|\lambda|\|v\|=\|v\|$，因此其谱落在单位圆上。

## 适用条件与失效情形

谱分解中的 $U^{-1}=U^\dagger$ 来自列正交归一，不能对一般特征向量矩阵照抄。
简并本征值允许在整个本征空间中酉旋转，单个本征向量因而不唯一，谱投影却唯一。
数值检查应分别检验 $H-H^\dagger$、$U^\dagger U-I$、$HU-U\Lambda$ 与重构，
不能只看“重构误差很小”就断定候选列为正交本征基。
下方练习给定候选谱数据，手写构造投影、重构与这些残差，避免调用特征求解库替代本章论证。

## 几何直觉

厄米算符在正交轴上作实数伸缩，负特征值会反向，零特征值会压塌方向。
酉算符只改变正交坐标方向及复相位，保留全部夹角和长度。
交互图将实对称二维矩阵的正交本征方向画在圆与椭圆上，负号由箭头方向显示；
实对称例子是复厄米几何在实平面内的截面，复相位不能全塞进普通二维箭头。

## 与量子力学的连接

有限维可观测量由厄米算符 H 表示。谱分解 $H=\sum_\lambda\lambda P_\lambda$ 给出测量公设：
归一化态 $|\psi\rangle$ 测得读数 $\lambda$ 的概率为
$$ p(\lambda)=\|P_\lambda\psi\|^2=\langle\psi|P_\lambda|\psi\rangle.\tag{8.5} $$
等号由 $P_\lambda^\dagger P_\lambda=P_\lambda$ 推出；
求和用 $\sum P_\lambda=I$ 得总概率为 1。测得该值后，若概率非零，
态更新为 $P_\lambda\psi/\sqrt{p(\lambda)}$。这条物理公设把线性代数中的投影指定为测量概率，
不是从矩阵定义单独推出的实验事实。

**完整手算算例：两个自旋分量。** 令 $|\psi\rangle=(\sqrt3/2,1/2)^T$，
其范数平方为 $3/4+1/4=1$。
$\sigma_z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}$，
$P_z^+=\begin{pmatrix}1&0\\0&0\end{pmatrix}$、$P_z^-=\begin{pmatrix}0&0\\0&1\end{pmatrix}$，
所以 $\sigma_z=P_z^+-P_z^-$。
$P_z^+\psi=(\sqrt3/2,0)^T$、$P_z^-\psi=(0,1/2)^T$，概率为 $3/4,1/4$。
$\sigma_x=\begin{pmatrix}0&1\\1&0\end{pmatrix}$ 的本征向量
$|+\rangle=(1,1)^T/\sqrt2$、$|-\rangle=(1,-1)^T/\sqrt2$，直接乘矩阵得到本征值 $+1,-1$。
对应投影为 $P_x^+=\frac12\begin{pmatrix}1&1\\1&1\end{pmatrix}$、
$P_x^-=\frac12\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$，相减重构 $\sigma_x$。
振幅 $\langle+|\psi\rangle=(\sqrt3+1)/(2\sqrt2)$，
$\langle-|\psi\rangle=(\sqrt3-1)/(2\sqrt2)$，模平方逐项展开为
$$ p_x^+=\frac{3+1+2\sqrt3}{8}=\frac{2+\sqrt3}{4}\approx0.933013,
\quad p_x^-=\frac{2-\sqrt3}{4}\approx0.066987. $$
两者相加为 1。期望值分别为 $\langle\sigma_z\rangle=3/4-1/4=1/2$、
$\langle\sigma_x\rangle=p_x^+-p_x^-=\sqrt3/2$。
物理自旋算符是 $S_z=\hbar\sigma_z/2,S_x=\hbar\sigma_x/2$，读数改为 $\pm\hbar/2$，概率不变。

**简并投影的数字核对。** 取 $H=\operatorname{diag}(2,2,5)$、
$\psi=(1,2,i)^T/\sqrt6$，读数 2 的整个本征空间是 $\operatorname{span}(e_1,e_2)$。
$P_2=\operatorname{diag}(1,1,0)$ 给 $p(2)=(1+4)/6=5/6$，
$P_5=\operatorname{diag}(0,0,1)$ 给 $p(5)=1/6$。
测到 2 后的态为 $(1,2,0)^T/\sqrt5$，其范数平方为 1。
若在简并空间改用 $q_\pm=(e_1\pm e_2)/\sqrt2$，
两投影矩阵分别是左上块 $\frac12\begin{pmatrix}1&\pm1\\\pm1&1\end{pmatrix}$，
求和仍为 $P_2$；两方向振幅为 $3/\sqrt{12},-1/\sqrt{12}$，
模平方相加 $9/12+1/12=5/6$，因此内部旋转不改变该能量读数的概率。
期望值也可两路核对：$2(5/6)+5(1/6)=5/2$，
直接算 $\psi^\dagger H\psi=(2+8+5)/6=5/2$。

对时间无关厄米 Hamilton H，由谱分解定义
$U(t)=\sum_\lambda e^{-i\lambda t/\hbar}P_\lambda$。
实 t、实非零 $\hbar$ 时相位模为 1，逐项乘积并用投影正交性得
$U^\dagger U=\sum_\lambda P_\lambda=I$，所以演化保总概率、保任意两态内积。
例如对 $H=\hbar\omega\sigma_z/2$，上例态演化为
$(\sqrt3 e^{-i\omega t/2}/2,e^{i\omega t/2}/2)^T$，范数平方仍为 1；
$\sigma_z$ 两概率保持 $3/4,1/4$，$\sigma_x$ 概率变为
$[1\pm(\sqrt3/2)\cos(\omega t)]/2$，展示保总概率与干涉概率随相对相位变化同时成立。
"""

REFERENCE = '''def spectral_audit(H, values, Q, psi):
    H = np.asarray(H, dtype=complex)
    values = np.asarray(values, dtype=float)
    Q = np.asarray(Q, dtype=complex)
    psi = np.asarray(psi, dtype=complex)
    n = len(values)
    if H.shape != (n, n) or Q.shape != (n, n) or psi.shape != (n,):
        raise ValueError("谱数据维数不匹配")
    reconstructed = np.zeros((n, n), dtype=complex)
    unique = np.unique(values)
    probabilities = np.zeros(len(unique))
    for g in range(len(unique)):
        P = np.zeros((n, n), dtype=complex)
        for j in range(n):
            if values[j] == unique[g]:
                P += Q[:, j:j+1] @ Q[:, j:j+1].conj().T
        reconstructed += unique[g] * P
        projected = P @ psi
        probabilities[g] = np.sum(np.abs(projected)**2)
    residuals = np.array([np.sqrt(np.sum(np.abs(H-H.conj().T)**2)),
                          np.sqrt(np.sum(np.abs(Q.conj().T@Q-np.eye(n))**2)),
                          np.sqrt(np.sum(np.abs(H@Q-Q*values)**2)),
                          np.sqrt(np.sum(np.abs(H-reconstructed)**2))])
    return reconstructed, residuals, probabilities
'''
CASES = '''def make_cases(rng):
    cases = []
    for k in range(32):
        n = 2 + k % 3
        Z = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
        Q, _ = np.linalg.qr(Z)
        values = rng.integers(-3, 4, size=n).astype(float)
        H = (Q*values) @ Q.conj().T
        psi = rng.normal(size=n) + 1j*rng.normal(size=n)
        psi /= np.sqrt(np.sum(np.abs(psi)**2))
        candidate = Q.copy()
        if k % 4 == 1:
            candidate[:, 0] *= 1.1
        elif k % 4 == 2:
            candidate = candidate[:, ::-1]
        cases.append((H, values, candidate, psi))
    cases.append((np.eye(2), np.ones(2), np.eye(2), np.array([1., 0.])))
    return cases
'''
WRITTEN = [
    {"id": "ch8-w1", "level": "基础", "problem_md": r"对 $H=\begin{pmatrix}2&i\\-i&2\end{pmatrix}$ 手算正交谱分解。", "solution_md": r"""
1. 依据二阶行列式，$\det(zI-H)=(z-2)^2-1$，根为 3、1。
2. 依据 $(H-3I)x=0$ 的第一行 $-x_1+ix_2=0$，取 $q_1=(i,1)^T/\sqrt2$；$(H-I)x=0$ 给 $q_2=(-i,1)^T/\sqrt2$。
3. 依据共轭内积，$q_1^\dagger q_2=((-i)(-i)+1)/2=0$，两列范数均为 1，故 U 酉。
4. 依据外积，$P_3=\frac12\begin{pmatrix}1&i\\-i&1\end{pmatrix},P_1=\frac12\begin{pmatrix}1&-i\\i&1\end{pmatrix}$，$3P_3+P_1=\begin{pmatrix}2&i\\-i&2\end{pmatrix}=H$。
""", "pitfalls_md": "复本征向量必须用共轭转置求内积；普通转置会改变正交判断。"},
    {"id": "ch8-w2", "level": "基础", "problem_md": r"证明 $U=\frac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$ 酉，并计算它对 $(1,i)^T/\sqrt2$ 的作用与范数。", "solution_md": r"""
1. 依据共轭转置，U 为实对称，$U^\dagger=U$。
2. 依据矩阵乘法，$U^\dagger U=\frac12\begin{pmatrix}2&0\\0&2\end{pmatrix}=I$，所以 U 酉。
3. 依据行列乘法，像为 $((1+i)/2,(1-i)/2)^T$。
4. 依据复数模平方，像的范数平方为 $(2+2)/4=1$，原态范数平方也为 $(1+1)/2=1$。
""", "pitfalls_md": "检查列长度还需检查交叉内积；长度相同不意味着两列正交。"},
    {"id": "ch8-w3", "level": "进阶", "problem_md": r"设复线性方阵 T 保所有向量范数。用极化恒等式证明 $T^\dagger T=I$。", "solution_md": r"""
1. 依据内积分配律，$D_1=\|x+y\|^2-\|x-y\|^2=4\operatorname{Re}\langle x,y\rangle$。
2. 依据第一变量共轭线性，$D_2=\|x+iy\|^2-\|x-iy\|^2=-4\operatorname{Im}\langle x,y\rangle$，故内积为 $(D_1-iD_2)/4$。
3. 依据 T 的复线性与保范数，将 Tx、Ty 代入上式后四个平方范数逐项恢复为原值，得 $\langle Tx,Ty\rangle=\langle x,y\rangle$。
4. 依据 $x=e_i,y=e_j$，$e_i^\dagger(T^\dagger T-I)e_j=0$ 给每个矩阵元为零，故所求等式成立。
""", "pitfalls_md": "内积第一变量共轭线性时，极化式虚部前的符号为负；不能直接用于反线性共轭映射。"},
    {"id": "ch8-w4", "level": "进阶", "problem_md": r"从正规性证明二阶上三角 $T=\begin{pmatrix}a&b\\0&d\end{pmatrix}$ 必有 b=0，并指出 Schur 与谱定理的差别。", "solution_md": r"""
1. 依据矩阵乘法，$(T^\dagger T)_{11}=|a|^2$。
2. 依据另一顺序的乘法，$(TT^\dagger)_{11}=|a|^2+|b|^2$。
3. 依据正规性两条目相等，$|b|^2=0$，复数模的正定性给 b=0。
4. 依据 Schur 构造，任意复方阵均可酉三角化；只有正规性额外强迫 b=0 才得到对角化，缺陷 Jordan 块的 b=1 正是在此受阻。
""", "pitfalls_md": "Schur 三角化没有承诺每列都是本征向量；非正规时超对角项可以存在。"},
    {"id": "ch8-w5", "level": "考试风格", "problem_md": r"对 $H=\operatorname{diag}(2,2,-1)$，态 $\psi=(1,i,\sqrt2)^T/2$，求每个能量结果的概率、测后态及期望值。", "solution_md": r"""
1. 依据谱投影应覆盖整个简并空间，$P_2=\operatorname{diag}(1,1,0),P_{-1}=\operatorname{diag}(0,0,1)$。
2. 依据投影范数平方，$p(2)=(1+1)/4=1/2$、$p(-1)=2/4=1/2$，总和为 1。
3. 依据投影归一化，测得 2 后态为 $(1,i,0)^T/\sqrt2$，测得 −1 后态为 $(0,0,1)^T$。
4. 依据谱期望值，$\langle H\rangle=2(1/2)-(1/2)=1/2$，直接计算 $\psi^\dagger H\psi=(2+2-2)/4=1/2$ 核对。
""", "pitfalls_md": "能量 2 的概率是前两个方向合计；不能把简并态当作两个不同能量读数。"},
    {"id": "ch8-w6", "level": "考试风格", "problem_md": r"初态 $\psi=(1,i)^T/\sqrt2$，$H=\hbar\omega\sigma_z/2$。求时刻 t 的 $\sigma_x$ 两结果概率。", "solution_md": r"""
1. 依据 H 的对角谱，$U(t)=\operatorname{diag}(e^{-i\omega t/2},e^{i\omega t/2})$。
2. 依据酉演化，$\psi(t)=(e^{-i\omega t/2},i e^{i\omega t/2})^T/\sqrt2$，范数平方为 1。
3. 依据 $\sigma_x$ 本征态 $|\pm\rangle=(1,\pm1)^T/\sqrt2$，振幅为 $(e^{-i\omega t/2}\pm i e^{i\omega t/2})/2$。
4. 依据模平方展开，正号交叉项为 $i e^{i\omega t}-i e^{-i\omega t}=-2\sin(\omega t)$，所以 $p_+=(1-\sin\omega t)/2,p_-=(1+\sin\omega t)/2$。
""", "pitfalls_md": "初态含 i，使干涉项是正弦；把它当作实叠加态会错写为余弦。"},
]
QUIZ = [
    {"q_md": "厄米矩阵的特征值可以是？", "options": ["只能正实数", "任意实数，包括负数与零", "只能整数", "任意复数"], "answer": 1, "explain_md": "x†Hx 等于自身共轭，只推出 λ 实；正定性是另一个条件。"},
    {"q_md": "同一简并本征空间内任选两个本征向量是否必正交？", "options": ["是", "只有三维才是", "否，需在空间内正交化", "简并空间不存在"], "answer": 2, "explain_md": "正交证明要除以不同本征值之差；相等时只能得到 0=0。"},
    {"q_md": "哪个条件使 Schur 上三角阵的超对角项全部为零？", "options": ["正规性", "可逆性", "迹为零", "行列式为一"], "answer": 0, "explain_md": "正规性使第一行与第一列平方范数相等，先消去第一行非对角项，再归纳。"},
    {"q_md": r"测得简并本征值 $\lambda$ 的概率是？", "options": [r"$\lambda^2$", "任取一列的振幅", r"$\operatorname{tr}H$", r"$\|P_\lambda\psi\|^2$"], "answer": 3, "explain_md": "Pλ 投向整个本征空间；概率是投影长度平方，与简并空间内部选基无关。"},
    {"q_md": "保范数推出保内积时，除了极化恒等式还必须使用？", "options": ["复线性", "所有条目为正", "特征值互异", "矩阵有逆的数值库"], "answer": 0, "explain_md": "必须把 Ux±iUy 写成 U(x±iy)；反线性共轭虽保范数，却不保持复内积。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch8-q{_i}"

CHAPTER = {
    "id": "ch8", "num": 8, "title": "厄米矩阵、酉矩阵与谱定理", "depends_on": ["ch5", "ch6"],
    "theory_md": THEORY, "viz": ["hermitian_spectrum"], "written": WRITTEN, "quiz": QUIZ,
    "code_exercises": [{
        "id": "ch8-c1", "title": "从谱投影重构并审计候选本征基",
        "prompt_md": r"实现 `spectral_audit(H, values, Q, psi)`。H 为 n 阶复厄米矩阵（n=2～4），values 为长度 n 的候选实本征值数组，重根用完全相同浮点数表示，Q 为 n×n 候选本征向量列，psi 是归一化态。候选 Q 可能有归一化或本征顺序错误，须用残差揭示。对 values 的每个不同值（升序）手写外积求和构造 P，返回三元组 `(重构H, residuals, weights)`；重构为 ΣλPλ，residuals 长度 4 依次为 H−H†、Q†Q−I、HQ−QΛ、H−重构H 的 Frobenius 范数，weights 为各 ‖Pλ psi‖²。仅在候选数据通过上述验证时，weights 才能解释为测量概率；错误候选仍须返回这些数值供审计。允许 @、共轭、np.unique、逐项绝对值平方求和与开方；不调用 linalg、eig/eigh、SVD 或现成范数。按简并空间聚合，不能把每列模平方直接作为不同读数概率。维数不匹配抛 ValueError。",
        "function_name": "spectral_audit", "starter": "def spectral_audit(H, values, Q, psi):\n    # TODO：分组构造投影，重构，逐项计算四类残差与投影权重\n    return H, np.zeros(4), np.zeros(1)\n",
        "banned": ["np.linalg", "np.linalg.eig", "np.linalg.eigh", "np.linalg.eigvals", "np.linalg.eigvalsh", "np.linalg.svd", "np.poly", "np.roots", "np.polynomial", "np.matrix", "np.asmatrix", "np.mat", "np.bmat", "np.matrixlib", "np.matlib", "scipy"],
        "reference": REFERENCE, "make_cases": CASES,
    }],
}
