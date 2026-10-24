"""第 5 章：特征值、重数与 Cayley–Hamilton 定理。"""
THEORY = r"""
## 动机与直觉

一般矩阵把各坐标耦合在一起。若找到只被乘以一个数的方向，沿这个方向反复作用就化为标量幂。
本章先回答这种方向何时存在，再区分“多项式有几个重根”与“有几个独立方向”，
最后证明矩阵自身满足其特征多项式，为化简高次幂提供工具。

## 定义

**定义 5.1（特征对）**：$A\in\mathbb F^{n\times n}$，若 $Av=\lambda v$ 且 $v\ne0$，
则 $\lambda$ 是特征值，$v$ 是对应特征向量。本征空间 $E_\lambda=\ker(A-\lambda I)$ 包含零向量。
**定义 5.2（特征多项式与重数）**：沿用行列式写
$$ p_A(z)=\det(zI-A),\qquad \det(A-zI)=(-1)^np_A(z).\tag{5.1} $$
$p_A$ 是首一 $n$ 次多项式。若 $p_A(z)=(z-\lambda)^a q(z)$ 且 $q(\lambda)\ne0$，
$a$ 为代数重数；$g=\dim E_\lambda$ 为几何重数。
**定义 5.3（矩阵多项式）**：若 $p(z)=\sum_{k=0}^n c_kz^k$，定义
$p(A)=\sum_{k=0}^n c_kA^k$，其中 $A^0=I$。这是标量系数多项式的代入。

## 核心推导

### 推导 5.1：特征值为何对应行列式的零点

动机：$Av=\lambda v$ 同时含未知数 $\lambda,v$，先把存在非零 $v$ 转成只含 $\lambda$ 的条件。
若 $Av=\lambda v$ 且 $v\ne0$，移项给 $(A-\lambda I)v=0$。
该矩阵有非零核，按第 3 章的核维数公式不满秩；第 4 章可逆等价链给 $\det(A-\lambda I)=0$。
反向，若该行列式为零，矩阵不可逆、秩小于 $n$，所以核维数 $n-r>0$。
取核中的一个非零向量 $v$，则 $(A-\lambda I)v=0$，移项恢复 $Av=\lambda v$。
因此
$$ \lambda\text{是特征值}\iff\det(A-\lambda I)=0\iff p_A(\lambda)=0.\tag{5.2} $$
$p_A$ 的首项来自置换定义中所有对角位置都取 $z$，系数是 $1$；任何其他选择最多含 $n-1$ 个 $z$，
所以它确实为首一 $n$ 次。沿每行提出 $-1$ 则给式 (5.1) 的符号关系。

条件与失效：特征向量必须非零，零向量满足所有 $\lambda$ 的齐次方程却不能证明有特征值。
数域也必须说明：实旋转矩阵可能没有实特征值，在复数域上特征多项式才保证完全分裂。

### 推导 5.2：不同特征值对应的向量线性无关

动机：找到多个特征值后，需要知道它们是否真的提供不同的独立方向。
设 $Av_j=\lambda_jv_j$，$v_j\ne0$，且 $\lambda_1,\ldots,\lambda_m$ 两两不同。
用 $m$ 归纳。$m=1$ 时，$c_1v_1=0$，取 $v_1$ 的某个非零坐标即得 $c_1=0$。
假设任意 $m-1$ 个这样的向量无关。对关系 $\sum_{j=1}^m c_jv_j=0$，
作用 $A-\lambda_mI$，利用线性性与特征方程得到
$$ \sum_{j=1}^{m-1} c_j(\lambda_j-\lambda_m)v_j=0. $$
归纳假设迫使每个 $c_j(\lambda_j-\lambda_m)=0$；因差不为零，$c_j=0$（$j<m$）。
回到原关系剩 $c_mv_m=0$，非零向量条件再给 $c_m=0$，归纳完成。
这个证明还适用于 $w_j\in E_{\lambda_j}$ 的和为零：若有非零 $w_j$，把这些非零项作为特征向量，
同一结论迫使它们全为零。因此不同本征空间的和是直和。

条件与失效：相同特征值时 $\lambda_j-\lambda_m=0$，不能除以它；
同一特征空间里可以选出无关向量，但不能仅凭相同特征方程保证无关。

### 推导 5.3：几何重数不超过代数重数

动机：多项式的重根数是否等于可用的方向数？先证明方向数至多有那么多，再给严格小于的例子。
取 $E_\lambda$ 的基 $v_1,\ldots,v_g$，用第 1 章基扩充方法补成全空间基，列矩阵记 $P$。
前 $g$ 个向量的像为 $\lambda v_j$，所以按第 2 章逐列写表示矩阵，得到
$$ B=P^{-1}AP=\begin{pmatrix}\lambda I_g&C\\0&D\end{pmatrix}. $$
由第 4 章乘法公式，
$\det(zI-B)=\det(P^{-1}(zI-A)P)=\det(zI-A)$，因为 $\det P^{-1}\det P=1$。
对分块上三角矩阵，置换求和中的非零项不能把下块行配给左块列；
下块行恰好占满右块列，于是上块行只能占左块列。置换分成两块，块间没有逆序，符号与乘积都分离。
因此
$$ p_A(z)=(z-\lambda)^g\det(zI-D). $$
代数重数是 $(z-\lambda)$ 的总因子数，上式已经贡献 $g$ 个，故 $1\le g\le a$。

**完整手算算例**：$J=\begin{pmatrix}2&1\\0&2\end{pmatrix}$。
$zI-J=\begin{pmatrix}z-2&-1\\0&z-2\end{pmatrix}$，所以 $p_J(z)=(z-2)^2$，代数重数为 $2$。
$J-2I=\begin{pmatrix}0&1\\0&0\end{pmatrix}$，齐次方程给 $v_2=0$，$v_1=s$。
本征空间为 $\operatorname{span}(1,0)^T$，几何重数为 $1$。
再算 $N=J-2I$，逐项乘法给 $N^2=0$，所以 $p_J(J)=J^2-4J+4I=0$。
实际 $J^2=\begin{pmatrix}4&4\\0&4\end{pmatrix}$，减去
$4J=\begin{pmatrix}8&4\\0&8\end{pmatrix}$ 再加 $4I$，四个条目都为零。

条件与失效：$g=a$ 需要额外结构；上三角矩阵也可能像 $J$ 一样只有一个方向。
不等式证明中，余块 $D$ 还可能贡献同一个根，所以不能从已有因子数直接断言等号。

### 推导 5.4：Cayley–Hamilton，从可对角化补到一般矩阵

动机：特征多项式能否约束矩阵本身？直接把行列式里的变量替成矩阵会混淆标量与矩阵，
应先按定义 5.3 证明 $p_A(A)=0$。

第一步，设 $A=P\Lambda P^{-1}$，$\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_n)$。
$A^0=PI P^{-1}$；若 $A^k=P\Lambda^kP^{-1}$，则
$A^{k+1}=P\Lambda^k(P^{-1}P)\Lambda P^{-1}=P\Lambda^{k+1}P^{-1}$。
归纳给所有非负整数幂。由有限求和分配律，
$$ p_A(A)=P\left(\sum_k c_k\Lambda^k\right)P^{-1}
=P\operatorname{diag}(p_A(\lambda_1),\ldots,p_A(\lambda_n))P^{-1}=0. $$
最后一步用每个 $\lambda_j$ 都是 $p_A$ 的根，即式 (5.2)。重根但可对角化的情况也在其中。

第二步，在复数域构造上三角表示。这里只调用代数基本定理：非恒定复多项式有根。
取 $p_A$ 的一个根及对应非零特征向量，把它扩成基，得到
$P_1^{-1}AP_1=\begin{pmatrix}\lambda&*\\0&B\end{pmatrix}$。
对 $(n-1)$ 阶 $B$ 重复这一步，归纳起点为一阶矩阵。
每次用 $\operatorname{diag}(1,Q)$ 只改变下块基，左下零块保留；有限次后得到 $A=STS^{-1}$，$T$ 上三角。
这给出了三角化的构造，而非假设所有矩阵都能对角化。

第三步，令 $d_1,\ldots,d_n$ 为互不相同的实数，定义
$T_\epsilon=T+\epsilon\operatorname{diag}(d_1,\ldots,d_n)$。
它的对角条目为 $t_{ii}+\epsilon d_i$。对每对 $i\ne j$，二者相等至多排除一个
$\epsilon=-(t_{ii}-t_{jj})/(d_i-d_j)$。
全部只有有限多个排除值，因此可选非零实数列 $\epsilon_m\to0$ 避开它们。
上三角行列式给 $p_{T_{\epsilon_m}}(z)=\prod_i(z-t_{ii}-\epsilon_m d_i)$，根两两不同。
推导 5.2 给 $n$ 个无关特征向量，所以 $A_m=ST_{\epsilon_m}S^{-1}$ 可对角化，且 $A_m\to A$。
按第一步，$p_{A_m}(A_m)=0$。

第四步，写 $p_X(z)=\sum_{k=0}^n c_k(X)z^k$。
置换定义说明每个 $c_k(X)$ 都是 $X$ 条目的有限乘积之和，因此连续；
$X^k$ 的每个条目也是条目的多项式。因此 $F(X)=\sum_k c_k(X)X^k$ 是矩阵条目的多项式映射。
由 $A_m\to A$，逐条目取极限给
$$ p_A(A)=F(A)=\lim_m F(A_m)=0.\tag{5.3} $$
注意极限时特征多项式的系数也随矩阵变化，不能始终使用某个固定 $p_{A_m}$。
对实矩阵，先把它视为复矩阵完成证明；$p_A$ 系数与 $A$ 的幂仍为实矩阵，所以结论回到实数域。
这完成由稠密性与多项式连续性补足的一般证明。

条件与失效：一般矩阵不能直接写出可逆特征向量矩阵，缺陷例 $J$ 在第一步就失效，
但第二至四步仍成立。该论证依赖有限维多项式与复数域三角化，不直接推广到无界微分算符。

### 推导 5.5：幂迭代为什么选择模最大的特征值

动机：高阶矩阵不宜展开特征多项式求根，希望只靠矩阵乘向量找到一个主方向。
设 $A$ 有特征向量基，$|\lambda_1|>|\lambda_j|$（$j>1$），且
$x_0=\sum_jc_jv_j$ 中 $c_1\ne0$。则
$$ A^kx_0=\lambda_1^k\left(c_1v_1+\sum_{j>1}c_j(\lambda_j/\lambda_1)^kv_j\right). $$
括号内第二项逐项趋零，方向趋向 $v_1$；每步归一化只改变整体标量，防止模长过大或过小。
若单位向量 $x$ 趋向该方向，$\lambda(x)=x^\dagger Ax/(x^\dagger x)$ 趋向 $\lambda_1$，
因为分子分母都是连续函数且分母非零。负或复特征值会导致相位变化，应比较射线与残差，而非未经对齐的坐标。

条件与失效：若初始向量没有主方向分量，展开中 $c_1=0$，迭代无法创造它；
若最大模并列，其他比值的模不小于一，衰减步骤失效，可能振荡而不收敛到单方向。
练习使用给定步数，输出当前近似，不把有限次迭代自动宣称为精确本征对。

## 适用条件与失效情形

实矩阵的非实根成共轭对，但在实向量空间中没有相应的实特征向量。
“有 $n$ 个根，按重数计”仅是多项式结论；能否凑齐 $n$ 个独立向量由几何重数控制。
Cayley–Hamilton 在缺陷矩阵上照样成立，不能把它误当成可对角化的判据。

## 几何直觉

实特征方向是一条不变直线：$Av=\lambda v$，$\lambda<0$ 时还翻转朝向。
`eigen_demo` 展示单位圆的像与这些直线；平面转动 $90^\circ$ 没有实不变直线，
但复向量 $(1,-i)^T$ 是 $\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ 的 $i$ 特征向量。
零特征值把一个方向压成零，这与第 4 章行列式为零一致。

## 与量子力学的连接

时间无关 Hamilton 算符的定态 Schrödinger 方程是
$\hat H|n\rangle=E_n|n\rangle$；在正交归一基下变为 $Hc_n=E_nc_n$。
例如自旋 $1/2$ 在磁场中的 $H=(\hbar\omega/2)\sigma_z$，
$\sigma_z=\operatorname{diag}(1,-1)$，能量为 $\pm\hbar\omega/2$，本征态为上下自旋。
若 $[H,S]=0$ 且 $H\psi=E\psi$，则 $H(S\psi)=S(H\psi)=E(S\psi)$，
所以对称操作 $S$ 把能量 $E$ 的本征空间映回自身；若产生与 $\psi$ 无关的态，便出现简并。
对易本身不强迫简并，例如 $S=I$ 也与任何 $H$ 对易。
简并微扰中的小矩阵也来自本征方程。设
$H(\epsilon)=H_0+\epsilon V$，$H_0,V$ 厄米，$E_0$ 是孤立的简并能级，
其本征空间有正交归一基 $|a\rangle$。取某条可微本征分支
$E(\epsilon)=E_0+\epsilon E^{(1)}+o(\epsilon)$，
$|\psi(\epsilon)\rangle=|\psi^{(0)}\rangle+\epsilon|\psi^{(1)}\rangle+o(\epsilon)$，
其中 $|\psi^{(0)}\rangle=\sum_b c_b|b\rangle$。
把两展开代回 $H(\epsilon)\psi(\epsilon)=E(\epsilon)\psi(\epsilon)$，
常数阶是已知方程，逐项比较一次阶得
$$ (H_0-E_0I)|\psi^{(1)}\rangle+V|\psi^{(0)}\rangle
=E^{(1)}|\psi^{(0)}\rangle. $$
左乘 $\langle a|$，由厄米性和 $H_0|a\rangle=E_0|a\rangle$ 得
$\langle a|(H_0-E_0I)=0$，所以未知的一阶态修正消去。
将零阶态展开代入后，剩下
$$ \sum_b\underbrace{\langle a|V|b\rangle}_{W_{ab}}c_b=E^{(1)}c_a. $$
这就是 $Wc=E^{(1)}c$ 的来源：W 的本征值给一阶能量修正，
其本征向量才给应选择的零阶组合，不能在简并空间中先任意固定一个基态就各算对角元。
$\overline{W_{ba}}=\langle a|V^\dagger|b\rangle=W_{ab}$，故 W 也是厄米矩阵。

具体取 $H_0=\operatorname{diag}(2,2,5)$、
$V=\begin{pmatrix}1&2&0\\2&1&0\\0&0&0\end{pmatrix}$，在 $E_0=2$ 空间中
$W=\begin{pmatrix}1&2\\2&1\end{pmatrix}$。
$\det(\mu I-W)=(\mu-1)^2-4=(\mu-3)(\mu+1)$，修正为 3、-1。
对 3，方程 $-2c_1+2c_2=0$ 给单位组合 $(1,1)^T/\sqrt2$；
对 -1，$2c_1+2c_2=0$ 给 $(1,-1)^T/\sqrt2$。
代回完整三阶矩阵直接核对
$H(\epsilon)(1,\pm1,0)^T/\sqrt2=(2+\epsilon\pm2\epsilon)(1,\pm1,0)^T/\sqrt2$，
所以本例两能量恰为 $2+3\epsilon$、$2-\epsilon$。
它们在 $\epsilon=0$ 合并为同一能级，非零小 $\epsilon$ 时分裂；两对角元同为 1 并没有给出分裂后的修正。

条件与失效：上述比较针对 $\epsilon=0$ 附近的可微本征分支与一阶系数；一般模型被略去的余项不能视为零。
本数字例中 V 不耦合第三能级，才使小矩阵答案同时成为精确答案。
若 W 本身仍简并，一阶方程不能选择其简并空间内的唯一组合，需要继续检查更高阶方程。
"""

REFERENCE = '''def power_iteration(A, x0, steps):
    A = np.asarray(A)
    x = np.array(x0, dtype=np.result_type(A, x0, float), copy=True)
    norm = np.sqrt(np.sum(np.abs(x)**2))
    if norm == 0:
        raise ValueError("初始向量非零")
    x /= norm
    for _ in range(steps):
        y = A @ x
        norm = np.sqrt(np.sum(np.abs(y)**2))
        if norm == 0:
            break
        x = y / norm
    index = int(np.argmax(np.abs(x)))
    x = x / (x[index] / abs(x[index]))
    value = np.sum(np.conj(x) * (A @ x)) / np.sum(np.abs(x)**2)
    return value, x
'''
CASES = '''def make_cases(rng):
    cases = [(np.array([[-3.]]), np.array([2.]), 15),
             (np.diag([4., -1.]), np.array([1., 2.]), 0),
             (np.diag([-4., 1.]), np.array([1., 2.]), 80),
             (np.zeros((2, 2)), np.array([1., -2.]), 8)]
    for index in range(28):
        n = int(rng.integers(2, 6))
        u = rng.normal(size=n)
        if index % 2:
            u = u + 1j*rng.normal(size=n)
        u = u / np.sqrt(np.sum(np.abs(u)**2))
        dominant = (-1 if index % 3 == 0 else 1) * 4.0
        A = np.eye(n) + (dominant-1) * np.outer(u, np.conj(u))
        x0 = u + .2*rng.normal(size=n)
        cases.append((A, x0, 80))
    return cases
'''
WRITTEN = [
    {"id": "ch5-w1", "level": "基础", "problem_md": r"求 $A=\begin{pmatrix}3&1\\0&2\end{pmatrix}$ 的特征值与各本征空间。", "solution_md": r"""
1. 依据定义 5.2，$p_A(z)=\det\begin{pmatrix}z-3&-1\\0&z-2\end{pmatrix}=(z-3)(z-2)$。
2. 依据式 (5.2)，特征值为 $3,2$。
3. 依据核定义，$(A-3I)v=0$ 给 $v_2=0$，所以 $E_3=\operatorname{span}(1,0)^T$。
4. 依据同一核定义，$(A-2I)v=0$ 给 $v_1+v_2=0$，所以 $E_2=\operatorname{span}(1,-1)^T$。
""", "pitfalls_md": "求本征空间要解齐次方程；特征值本身不是特征向量。"},
    {"id": "ch5-w2", "level": "基础", "problem_md": r"求 $J=\begin{pmatrix}2&1\\0&2\end{pmatrix}$ 的代数、几何重数。", "solution_md": r"""
1. 依据三角行列式，$p_J(z)=(z-2)^2$。
2. 依据代数重数定义，$2$ 的代数重数是 $2$。
3. 依据齐次方程，$(J-2I)v=0$ 只有条件 $v_2=0$，解为 $(s,0)^T$。
4. 依据本征空间维数定义，几何重数为 $1$；所以这里严格小于代数重数。
""", "pitfalls_md": "重根数统计多项式因子，几何重数统计独立向量，两者不要混用。"},
    {"id": "ch5-w3", "level": "进阶", "problem_md": r"若 $Au=2u,Av=5v$ 且 $u,v\ne0$，逐步证明二者无关。", "solution_md": r"""
1. 依据无关性定义，设 $au+bv=0$，目标是证明 $a=b=0$。
2. 依据线性性，对关系作用 $A-5I$，得 $-3au+0bv=0$。
3. 依据 $u\ne0$ 与 $-3\ne0$，有 $a=0$。
4. 依据原关系，剩下 $bv=0$，再由 $v\ne0$ 得 $b=0$，故无关。
""", "pitfalls_md": "必须作用在线性关系的整个等式上；不能从两个向量不相等就断言无关。"},
    {"id": "ch5-w4", "level": "进阶", "problem_md": r"对 $A=\begin{pmatrix}1&1\\0&2\end{pmatrix}$ 用 Cayley–Hamilton 求 $A^3$。", "solution_md": r"""
1. 依据三角行列式，$p_A(z)=(z-1)(z-2)=z^2-3z+2$。
2. 依据式 (5.3)，$A^2-3A+2I=0$，所以 $A^2=3A-2I$。
3. 依据矩阵乘法分配律，$A^3=3A^2-2A=3(3A-2I)-2A=7A-6I$。
4. 依据逐条目运算，$7A-6I=\begin{pmatrix}1&7\\0&8\end{pmatrix}$，即所求。
""", "pitfalls_md": "常数项代入矩阵多项式后是常数乘 I，不能把标量加到每个条目。"},
    {"id": "ch5-w5", "level": "考试风格", "problem_md": r"$A=\operatorname{diag}(-3,1)$，$x_0=(1,1)^T$。求幂迭代的极限射线与 Rayleigh 商；若 $x_0=(0,1)^T$ 呢？", "solution_md": r"""
1. 依据对角矩阵幂，$A^kx_0=((-3)^k,1)^T$。
2. 依据归一化定义，其单位向量为 $((-3)^k,1)^T/\sqrt{9^k+1}$，首分量交替变号，射线趋向 $\operatorname{span}e_1$。
3. 依据 Rayleigh 商定义，值为 $(-3\cdot9^k+1)/(9^k+1)$，除以 $9^k$ 后极限为 $-3$。
4. 依据 $Ae_2=e_2$，若初始为 $e_2$，各步保持 $e_2$，商始终为 $1$；主方向系数为零使推导 5.5 的条件失效。
""", "pitfalls_md": "主特征值指模最大，不是数值最大；负值会使向量符号交替，但射线可收敛。"},
    {"id": "ch5-w6", "level": "考试风格", "problem_md": r"二重简并子空间中的微扰矩阵 $W=\begin{pmatrix}0&\delta\\\delta&0\end{pmatrix}$，$\delta>0$。求一阶修正与归一化态。", "solution_md": r"""
1. 依据特征多项式定义，$p_W(z)=z^2-\delta^2$，根为 $\delta,-\delta$。
2. 依据 $(W-\delta I)c=0$，有 $c_1=c_2$；由归一化 $|c_1|^2+|c_2|^2=1$，可取 $(1,1)^T/\sqrt2$。
3. 依据 $(W+\delta I)c=0$，有 $c_1=-c_2$，同样可取 $(1,-1)^T/\sqrt2$。
4. 依据简并微扰的子空间本征问题，若 $H=H_0+\epsilon V$，这两个组合的能量为 $E_0\pm\epsilon\delta+o(\epsilon)$，简并在一阶分裂。
""", "pitfalls_md": "简并情况下不能只读取原基中的对角矩阵元 0；必须先在简并子空间对角化 W。"},
]
QUIZ = [
    {"q_md": "特征向量的必要要求是？", "options": ["必须非零", "必须实数", "长度必须为 1", "每个分量非零"], "answer": 0, "explain_md": "归一化可选，但非零是定义要求，否则任意 λ 都满足方程。"},
    {"q_md": "几何重数 g 与代数重数 a 的关系是？", "options": [r"$g>a$", r"$1\le g\le a$", r"$g=a$ 恒成立", r"$g=0$"], "answer": 1, "explain_md": "选本征空间基扩充后，特征多项式至少含 g 个 (z-λ) 因子。"},
    {"q_md": "不同特征值的非零特征向量必然？", "options": ["正交", "长度相同", "线性无关", "每个分量相同"], "answer": 2, "explain_md": "推导 5.2 的归纳消去证明无关；一般矩阵没有正交性保证。"},
    {"q_md": "Cayley–Hamilton 适用于？", "options": ["仅对角矩阵", "仅可逆矩阵", "仅无重根矩阵", "任意有限方阵"], "answer": 3, "explain_md": "先在可对角化情形证明，再由复三角化构造稠密逼近，并用多项式连续性取极限。"},
    {"q_md": "幂迭代可能不能选出主方向的原因是？", "options": ["主特征值为负数", "初始态在主方向上的系数为零", "每步归一化", "使用复数"], "answer": 1, "explain_md": "矩阵幂只乘已有的本征分量，不能凭空产生缺失的主分量。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch5-q{_i}"

CHAPTER = {
    "id": "ch5", "num": 5, "title": "特征值与特征向量", "depends_on": ["ch3", "ch4"],
    "theory_md": THEORY, "viz": ["eigen_demo"], "written": WRITTEN, "quiz": QUIZ,
    "code_exercises": [{
        "id": "ch5-c1", "title": "幂迭代与相位对齐",
        "prompt_md": r"实现 `power_iteration(A, x0, steps)`，返回 `(value, vector)`。先用手写模长归一化非零 x0，再恰好迭代 steps 次：计算 y=A@x，若 y 模长为零则停止并保留当前 x，否则令 x=y/模长。最后选模最大的分量（并列取最先）将整体相位除去，使该分量为正实数；返回 Rayleigh 商及此单位向量。模长用 abs、sum、sqrt 实现。A 为尺度适中的实或复 Hermitian 方阵，steps 为非负整数；一般用例具有唯一最大模特征值和非零主分量，另含一维、零步、零矩阵边界。不得调用 linalg（含 eig/eigvals/eigh/eigvalsh/norm）、np.poly/roots/polynomial 或 scipy，禁用 numpy.matrix 系列封装。固定步数及相位约定用于自动比较，收敛还需检查残差。",
        "function_name": "power_iteration", "starter": "def power_iteration(A, x0, steps):\n    # TODO：矩阵乘向量、手写归一化、相位对齐、Rayleigh 商\n    return 0.0, np.zeros_like(x0)\n",
        "banned": ["np.matrix", "np.asmatrix", "np.mat", "np.bmat", "np.matrixlib", "np.matlib", "np.matrix.I", "np.matrix.getI", "np.poly", "np.roots", "np.polynomial", "np.linalg", "np.linalg.eig", "np.linalg.eigvals", "np.linalg.eigh", "np.linalg.eigvalsh", "np.linalg.norm", "scipy"],
        "reference": REFERENCE, "make_cases": CASES,
    }],
}
