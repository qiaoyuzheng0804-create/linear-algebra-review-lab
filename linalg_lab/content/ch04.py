"""第 4 章：行列式，从置换到可逆性。"""
THEORY = r"""
## 动机与直觉

消元回答能否解方程，行列式把方阵是否压塌一个方向浓缩成一个标量。
本章不先假设面积公式或乘法公式，而从置换求和构造它们；符号记录定向，绝对值记录实空间的体积比。

## 定义

**定义 4.1（置换与奇偶性）**：$S_n$ 是 $1,\ldots,n$ 的双射全体。
$N(\sigma)=\#\{(i,j):i<j,\sigma(i)>\sigma(j)\}$ 称逆序数，
$\operatorname{sgn}\sigma=(-1)^{N(\sigma)}$。置换复合从右向左作用。
**定义 4.2（Leibniz 定义）**：对 $\mathbb F=\mathbb R$ 或 $\mathbb C$ 上的方阵，
$$ \det A=\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)\prod_{i=1}^n a_{i,\sigma(i)}.\tag{4.1} $$
约定零阶行列式为 $1$。删去第 $i$ 行第 $j$ 列得到 $M_{ij}$，
$C_{ij}=(-1)^{i+j}\det M_{ij}$ 为代数余子式；$\operatorname{adj}A=(C_{ji})$。

## 核心推导

### 推导 4.1：置换符号如何产生交错多重线性

动机：消元用行操作，要从定义计算每种操作如何影响行列式，不能把换行变号当作未经证明的规则。
先看交换置换列表中相邻的两个位置：二者的相互次序翻转，贡献的逆序数改变 $1$ 或 $-1$；
任何第三个位置与它们产生的逆序数总和不变，故总奇偶性翻转。
交换任意位置 $r<s$ 可以先把第 $r$ 项相邻移动到 $s$，再把原第 $s$ 项移动回 $r$，
共 $(s-r)+(s-r-1)=2(s-r)-1$ 次相邻交换，次数为奇数。
所以对位置对换 $\tau=(r\ s)$，$\operatorname{sgn}(\sigma\circ\tau)=-\operatorname{sgn}\sigma$。
把任意置换 $\rho$ 用相邻交换排回恒等列表，每次逆序数奇偶翻转；若用了 $m$ 次，
则 $(-1)^m=\operatorname{sgn}\rho$。记这些交换为 $\tau_1,\ldots,\tau_m$，于是
$\rho\tau_1\cdots\tau_m=\mathrm{id}$，且 $\sigma\rho\tau_1\cdots\tau_m=\sigma$。
从 $\sigma\rho$ 开始逐次复合，也翻转 $m$ 次，故
$\operatorname{sgn}\sigma=(-1)^m\operatorname{sgn}(\sigma\rho)$，再乘 $(-1)^m$ 得
$\operatorname{sgn}(\sigma\rho)=\operatorname{sgn}\sigma\operatorname{sgn}\rho$。
取 $\rho=\sigma^{-1}$，得逆置换与原置换同号。

令 $B$ 为 $A$ 交换第 $r,s$ 行的结果，故 $b_{i,j}=a_{\tau(i),j}$。
定义 (4.1) 的一项乘积变为
$$ \prod_i a_{\tau(i),\sigma(i)}=\prod_k a_{k,\sigma(\tau(k))}.
$$
令 $\pi=\sigma\circ\tau$，这是 $S_n$ 到自身的双射，且
$\operatorname{sgn}\sigma=-\operatorname{sgn}\pi$。逐项重编号求和，得到 $\det B=-\det A$。
若两行相同，交换后矩阵不变，所以 $d=-d$；实数或复数中 $2\ne0$，因此 $d=0$。
这就是交错性。

固定其他行，第 $r$ 行写成 $\alpha u+\beta v$。每个置换项恰含一个第 $r$ 行元素，
其因子 $\alpha u_{\sigma(r)}+\beta v_{\sigma(r)}$ 按分配律展开，剩余因子不变。
有限求和可拆成两和，故
$$ D(\ldots,\alpha u+\beta v,\ldots)=\alpha D(\ldots,u,\ldots)+\beta D(\ldots,v,\ldots). $$
逐行成立，称多重线性。因此一行乘 $c$ 使行列式乘 $c$；
$R_r\leftarrow R_r+cR_s$ 时，多出的项有两行相同，为零，行列式保持。
再在 (4.1) 中令 $\pi=\sigma^{-1}$，将 $\prod_i a_{\sigma(i),i}$ 重编号为
$\prod_j a_{j,\pi(j)}$，并用逆置换同号，得 $\det A^T=\det A$。
于是转置后的行性质逐条给出列的多重线性与交错性。

条件与失效：线性只对单独一行或一列成立，不表示 $\det(A+B)=\det A+\det B$；
同时改变两行会在分配律中产生交叉项，不能只保留两个行列式。

### 推导 4.2：按行展开与伴随恒等式

动机：置换有 $n!$ 项，希望把它分成更小的行列式，并把计算工具接到解方程上。
先把第 $i$ 行用 $i-1$ 次相邻交换移到首行；对固定列 $j$，再用 $j-1$ 次相邻列交换把该列移到首列。
其余行列的相对顺序保留。新矩阵中所有满足 $\sigma(1)=1$ 的项，删掉首行首列后
与 $S_{n-1}$ 一一对应，首个值 $1$ 不产生逆序，符号就是小置换的符号。
把两次搬动的符号还原，原矩阵所有满足 $\sigma(i)=j$ 的项总和是
$(-1)^{i-1+j-1}a_{ij}\det M_{ij}=a_{ij}C_{ij}$。
每个置换恰有一个 $j=\sigma(i)$，所以这些组互不重叠且覆盖全部项：
$$ \det A=\sum_{j=1}^n a_{ij}C_{ij}.\tag{4.2} $$

现在 $(A\operatorname{adj}A)_{ik}=\sum_j a_{ij}C_{kj}$。
若 $i=k$，按第 $i$ 行展开等于 $\det A$；若 $i\ne k$，把 $A$ 的第 $k$ 行替换为第 $i$ 行，
其第 $k$ 行余子式未变，上式就是这个新矩阵的行列式，而它有重复行，所以为零。
对转置矩阵作同样的行展开，再转置回去，得到另一侧恒等式。因此
$$ A\operatorname{adj}A=\operatorname{adj}A\,A=(\det A)I.\tag{4.3} $$

条件与失效：余子式符号由行列搬动次数决定，不能漏掉 $(-1)^{i+j}$。
恒等式对奇异矩阵仍成立，但此时不能除以 $\det A$ 构造逆。

### 推导 4.3：经初等矩阵分解证明乘法公式

动机：复合变换的缩放比应该相乘；先证明代数版本，再解释几何，避免循环论证。
对上三角阵 $U$，非零置换项必须满足 $\sigma(i)\ge i$，否则会取到下三角零元。
因为 $\sum_i\sigma(i)=\sum_i i$，这些不等式只有全部取等才能成立，故
$\det U=\prod_i u_{ii}$，特别 $\det I=1$。
任一初等矩阵 $E$ 是对 $I$ 作一次行操作，左乘 $E$ 对任何 $X$ 作相同行操作。
推导 4.1 给出的因子分别为 $-1,c,1$，而作用在 $I$ 上时它们恰为 $\det E$，所以
$$ \det(EX)=\det E\det X. $$
这只用单次操作性质，还未使用一般乘法公式。

若 $A$ 可逆，第 3 章消元得到 $E_m\cdots E_1A=I$，故
$A=E_1^{-1}\cdots E_m^{-1}$。每个逆仍是初等矩阵。
把单次公式连续用于 $AB=E_1^{-1}\cdots E_m^{-1}B$，得到
$$ \det(AB)=\left(\prod_{j=1}^m\det E_j^{-1}\right)\det B. $$
在同一串公式里令 $B=I$，括号就是 $\det A$，乘法公式成立。
若 $A$ 奇异，消元所得阶梯形有零行，行操作因子非零，故 $\det A=0$。
记消元总矩阵为 $E$，则 $EAB=(EA)B$ 也有同一零行，所以 $\det(EAB)=0$。
逐个初等操作的因子仍非零，倒推 $\det(AB)=0=\det A\det B$。
两种情形合起来证明
$$ \det(AB)=\det A\det B.\tag{4.4} $$

条件与失效：初等分解为可逆因子之积只适用于可逆 $A$；奇异情形需要单独的零行论证。
乘法公式要求同阶方阵，不能直接把一般矩形乘积拆成两个未定义的行列式。

### 推导 4.4：可逆等价链与克拉默法则

动机：把数值判据与唯一解公式联系起来，并说明为何分母不能为零。
若 $\det A\ne0$，式 (4.3) 两边除以它，给双侧逆
$A^{-1}=\operatorname{adj}A/\det A$。若 $A$ 可逆，式 (4.4) 给
$1=\det(AA^{-1})=\det A\det A^{-1}$，因此 $\det A\ne0$。
结合第 3 章的主元与核维数结论，得到完整等价链
$$ \det A\ne0\iff A\text{可逆}\iff\operatorname{rank}A=n
\iff\ker A=\{0\}\iff\text{任意 }b\text{ 的 }Ax=b\text{ 有唯一解}. $$
其中满秩给 $n$ 个主元，可把 $[A\mid I]$ 化为 $[I\mid A^{-1}]$；
反向可逆映射核为零，核维数 $n-\operatorname{rank}A$ 迫使满秩。
“任意右端唯一解”正是线性映射双射：存在性给满射，$b=0$ 的唯一性给单射。

设 $A_j(b)$ 将第 $j$ 列改成 $b$，而 $Ax=b=\sum_k x_kA_{:k}$。
对被替换列用多重线性，
$\det A_j(b)=\sum_k x_k\det A_j(A_{:k})$。
当 $k\ne j$，第 $j,k$ 列重复，该项为零；当 $k=j$，该项为 $x_j\det A$。
所以在 $\det A\ne0$ 时，
$$ x_j=\frac{\det A_j(b)}{\det A}.\tag{4.5} $$

**完整手算算例**：$A=\begin{pmatrix}2&1\\1&3\end{pmatrix}$，$b=(5,5)^T$。
二阶置换为恒等与对换，符号分别 $+1,-1$，故 $\det A=2\cdot3-1\cdot1=5$。
替换第一列给 $\det\begin{pmatrix}5&1\\5&3\end{pmatrix}=15-5=10$；
替换第二列给 $\det\begin{pmatrix}2&5\\1&5\end{pmatrix}=10-5=5$。
由式 (4.5) 得 $x_1=2,x_2=1$。代回两行：$2\cdot2+1=5$、$2+3\cdot1=5$，核对完成。

条件与失效：$\det A=0$ 时分母不可除；此时可能无解，也可能有自由变量，不能由一个零行列式区分二者。
克拉默法则适合理论分析与小阶手算，大阶实现应采用消元而非反复计算余子式。

## 适用条件与失效情形

本章代数结论适用于实或复方阵；普通面积与体积意义使用实矩阵。
非零但很小的行列式仍代表精确可逆，浮点消元却可能把小主元当零；差别发生在主元阈值判断。
行列式数值大小随整体单位缩放成 $\det(cA)=c^n\det A$，不能单凭其大小衡量条件数。

## 几何直觉

令两列 $u=(a,c)^T,v=(b,d)^T$。当 $u\ne0$，底长为 $\sqrt{a^2+c^2}$，
$v$ 在单位法向量 $(-c,a)^T/\sqrt{a^2+c^2}$ 上的有向高度为
$(ad-bc)/\sqrt{a^2+c^2}$，底乘有向高等于 $ad-bc=\det A$。
$u=0$ 时两者都是零。以上面算例为例，底长 $\sqrt5$，有向高 $(6-1)/\sqrt5=\sqrt5$，面积为 $5$。
交换两列翻转定向，普通面积仍是 $5$。
推广到 $n$ 维：行加法是剪切，切片的高度和平行截面体积保持；交换坐标轴保持普通体积，
缩放一坐标轴 $c$ 倍使体积乘 $|c|$。按初等操作分解，这些因子相乘得到 $|\det A|$。
奇异变换把像压进低维子空间，$n$ 维体积为零。
对有限体积的可测区域因此有 $\mathrm{Vol}(A\Omega)=|\det A|\mathrm{Vol}(\Omega)$；`det_area` 展示单位正方形的像。
若 A 可逆，该式也允许区域体积为无穷；若 A 奇异且区域体积无穷，像的 n 维体积仍为零，
但右侧会出现未定义的 $0\cdot\infty$，此时应直接使用低维像结论，而不写这个乘积式。

## 与量子力学的连接

对 $N$ 个正交归一单粒子轨道 $\phi_j$，费米子态写成 Slater 行列式
$$ \Psi(x_1,\ldots,x_N)=\frac1{\sqrt{N!}}\det[\phi_j(x_i)]. $$
交换两个粒子坐标交换两行，推导 4.1 给 $\Psi\mapsto-\Psi$；两个粒子占据同一轨道时两列相同，波函数为零。
归一化因子来自展开的 $N!$ 个乘积态：不同置换对应至少一处正交轨道，交叉内积为零；
每项范数为一，总范数平方为 $N!$，故除以 $\sqrt{N!}$。
把反对称性落到两个离散位置：取正交轨道 $\phi_1=(1,0)^T,\phi_2=(0,1)^T$，
位置标号为 0、1。双粒子振幅
$\Psi(x_1,x_2)=[\phi_1(x_1)\phi_2(x_2)-\phi_2(x_1)\phi_1(x_2)]/\sqrt2$。
逐对代入得 $\Psi(0,0)=0$、$\Psi(0,1)=1/\sqrt2$、
$\Psi(1,0)=-1/\sqrt2$、$\Psi(1,1)=0$，概率和为 $0+1/2+1/2+0=1$。
交换位置时振幅变号而模平方不变；若把第二轨道也取为 $(1,0)^T$，
四个行列式都因两列相同成为零，无法再归一化成双费米子态。

坐标变换 $x=x(q)$ 则由局部雅可比 $J_{ij}=\partial x_i/\partial q_j$ 给
$d^nx=|\det J|d^nq$，概率归一化变为 $\int|\psi(x(q))|^2|\det J|d^nq=1$。
三维球坐标的因子 $r^2\sin\theta$ 正来自此式；若希望在平直测度 $d^nq$ 下表示态，
可定义 $\widetilde\psi(q)=|\det J|^{1/2}\psi(x(q))$。这要求在光滑可逆坐标片内使用换元。
"""

REFERENCE = '''def determinants(A):
    A = np.asarray(A)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("需要方阵")
    def expand(M):
        n = M.shape[0]
        if n == 0:
            return 1.0
        total = 0.0
        for j in range(n):
            minor = np.delete(np.delete(M, 0, axis=0), j, axis=1)
            total += (-1)**j * M[0, j] * expand(minor)
        return total
    R = np.array(A, dtype=np.result_type(A, float), copy=True)
    product = 1.0
    for j in range(len(R)):
        pivot = j + int(np.argmax(np.abs(R[j:, j])))
        if abs(R[pivot, j]) <= 1e-10:
            return expand(A), 0.0
        if pivot != j:
            R[[j, pivot]] = R[[pivot, j]]
            product = -product
        value = R[j, j]
        product *= value
        for i in range(j+1, len(R)):
            R[i, j:] -= (R[i, j] / value) * R[j, j:]
    return expand(A), product
'''
CASES = '''def make_cases(rng):
    cases = [(np.zeros((0, 0)),), (np.array([[0., 2.], [3., 4.]]),),
             (np.array([[1., 2.], [2., 4.]]),), (np.array([[2j]]),)]
    for index in range(28):
        n = int(rng.integers(1, 6))
        A = rng.integers(-4, 5, size=(n, n)).astype(float)
        if index % 3 == 0:
            A = A + 1j*rng.integers(-2, 3, size=(n, n))
        if index % 4 == 0 and n > 1:
            A[-1] = A[0]
        cases.append((A,))
    return cases
'''
WRITTEN = [
    {"id": "ch4-w1", "level": "基础", "problem_md": r"从置换定义计算 $\det\begin{pmatrix}3&2\\1&4\end{pmatrix}$，再交换两行。", "solution_md": r"""
1. 依据定义 4.1，$S_2$ 只有列表 $(1,2)$ 与 $(2,1)$，逆序数为 $0,1$。
2. 依据定义 4.2，两项是 $3\cdot4$ 与 $-2\cdot1$，和为 $10$。
3. 依据交换后的条目重新代入定义，行列式为 $1\cdot2-4\cdot3=-10$。
4. 依据推导 4.1 的对换奇偶性，结果必须为原值的负数，与 $-10$ 相符。
""", "pitfalls_md": "二阶反对角项有负号；行列式不是四个元素的乘积。"},
    {"id": "ch4-w2", "level": "基础", "problem_md": r"用按行展开计算 $\det\begin{pmatrix}1&2&0\\0&3&4\\0&0&5\end{pmatrix}$。", "solution_md": r"""
1. 依据式 (4.2)，选第三行展开，前两项的系数为零，只剩 $5C_{33}$。
2. 依据代数余子式定义，$C_{33}=(-1)^6\det\begin{pmatrix}1&2\\0&3\end{pmatrix}$。
3. 依据二阶置换公式，该小行列式为 $1\cdot3-2\cdot0=3$。
4. 依据乘法运算，原值为 $5\cdot3=15$；也符合推导 4.3 的三角主对角线乘积。
""", "pitfalls_md": "按某一行展开要使用那一行的所有代数余子式，不能混用列编号。"},
    {"id": "ch4-w3", "level": "进阶", "problem_md": r"用消元算 $A=\begin{pmatrix}0&2&1\\1&1&0\\2&0&1\end{pmatrix}$ 的行列式。", "solution_md": r"""
1. 依据推导 4.1，交换前两行，得 $B=\begin{pmatrix}1&1&0\\0&2&1\\2&0&1\end{pmatrix}$，且 $\det B=-\det A$。
2. 依据加行倍数保持行列式，$R_3\leftarrow R_3-2R_1$，第三行变为 $(0,-2,1)$。
3. 依据同一规则，$R_3\leftarrow R_3+R_2$，第三行变为 $(0,0,2)$。
4. 依据三角矩阵公式，$\det B=1\cdot2\cdot2=4$；还原交换的符号，$\det A=-4$。
""", "pitfalls_md": "记录换行次数；如果还把主元行归一化，必须另记行缩放因子。"},
    {"id": "ch4-w4", "level": "进阶", "problem_md": r"用克拉默法则解 $2x+y=5,\ x+3y=5$，并说明面积意义。", "solution_md": r"""
1. 依据式 (4.1)，系数行列式为 $2\cdot3-1=5\ne0$，满足克拉默法则条件。
2. 依据式 (4.5)，替换第一列后的值为 $5\cdot3-5=10$，所以 $x=10/5=2$。
3. 依据式 (4.5)，替换第二列后的值为 $2\cdot5-5=5$，所以 $y=5/5=1$。
4. 依据几何小节的底乘高公式，系数矩阵将单位正方形映为面积 $|5|=5$ 的平行四边形，正号表示定向保持。
""", "pitfalls_md": "右端替换的是一列，不能替换一行；普通面积取绝对值。"},
    {"id": "ch4-w5", "level": "考试风格", "problem_md": r"设 $A$ 是三阶可逆矩阵，$\det A=-2$，求 $\det(3A^{-1}A^T)$ 并说明每个因子。", "solution_md": r"""
1. 依据式 (4.4)，$1=\det(AA^{-1})$，所以 $\det A^{-1}=1/(-2)=-1/2$。
2. 依据推导 4.1 中逆置换重编号，$\det A^T=\det A=-2$。
3. 依据逐行线性，三阶矩阵整体乘 $3$ 等于三行各乘 $3$，因子为 $3^3=27$。
4. 依据式 (4.4)，所求为 $27\cdot(-1/2)\cdot(-2)=27$。
""", "pitfalls_md": "整体标量的指数是矩阵阶数；转置不会像交换一对行那样变号。"},
    {"id": "ch4-w6", "level": "考试风格", "problem_md": r"两费米子占据归一且正交的轨道 $\phi,\chi$。写出归一化反对称态，证明同轨道时为零。", "solution_md": r"""
1. 依据二阶行列式定义，$\Psi(x_1,x_2)=c[\phi(x_1)\chi(x_2)-\chi(x_1)\phi(x_2)]$。
2. 依据单粒子正交归一性，两个乘积态各范数为 $1$，交叉内积含 $\langle\phi,\chi\rangle=0$，所以 $\|\Psi\|^2=2|c|^2$。
3. 依据归一化要求并选择正实系数，取 $c=1/\sqrt2$；交换 $x_1,x_2$ 后两项交换，整体变号。
4. 依据列交错性，若 $\chi=\phi$，两列相同，行列式为零；此时不能再归一化为物理态。
""", "pitfalls_md": "同轨道情形给零态，不能继续把它当归一化态；轨道不正交时原归一化因子也需重算。"},
]
QUIZ = [
    {"q_md": "交换两行使行列式如何变化？", "options": ["不变", "变为负值的绝对值", "乘以 -1", "变为零"], "answer": 2, "explain_md": "置换项以位置对换重编号；奇偶性翻转，所以每项符号反转。"},
    {"q_md": r"$R_2\leftarrow R_2+3R_1$ 的行列式因子是？", "options": ["3", "1", "-1", "0"], "answer": 1, "explain_md": "多重线性展开多出一个含重复行的行列式，它为零。"},
    {"q_md": r"$\det A=0$ 可以推出？", "options": ["方程必无解", "方程必有无穷多解", "A 可逆", "A 的核含非零向量"], "answer": 3, "explain_md": "方阵不满秩，核维数为 n-r>0；右端仍决定非齐次方程是否相容。"},
    {"q_md": r"$\det A=-3,\det B=2$，$\det(AB)$ 为？", "options": ["-6", "6", "-1", "5"], "answer": 0, "explain_md": "同阶方阵的乘法公式给 (-3)×2=-6。"},
    {"q_md": "实二维矩阵的普通面积缩放因子是？", "options": [r"$\det A$", r"$|\det A|$", r"$\operatorname{tr}A$", "秩"], "answer": 1, "explain_md": "底乘有向高得到 det A，普通面积丢去定向符号，所以取绝对值。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch4-q{_i}"

CHAPTER = {
    "id": "ch4", "num": 4, "title": "行列式", "depends_on": ["ch2", "ch3"],
    "theory_md": THEORY, "viz": ["det_area"], "written": WRITTEN, "quiz": QUIZ,
    "code_exercises": [{
        "id": "ch4-c1", "title": "展开与消元，两条路径计算行列式",
        "prompt_md": r"实现 `determinants(A)`，返回 `(d_expand, d_eliminate)`。第一项必须从首行余子式递归展开，第二项用部分主元消元及主元乘积，记录换行符号。输入为实或复方阵，阶数 0～5，尺度适中；零阶值为 1。消元以 $10^{-10}$ 判零，遇零主元返回第二项为 0。随机用例避开近奇异阈值边界。不得调用任何 `linalg`（含 det、slogdet）、np.poly/roots/polynomial 或 scipy，禁用 numpy.matrix 系列封装。两条算法在精确算术中由推导 4.2～4.3 相等，浮点中以容差比较。",
        "function_name": "determinants", "starter": "def determinants(A):\n    # TODO：递归展开与带换行符号的主元乘积分别计算\n    return 0.0, 0.0\n",
        "banned": ["np.matrix", "np.asmatrix", "np.mat", "np.bmat", "np.matrixlib", "np.matlib", "np.matrix.I", "np.matrix.getI", "np.poly", "np.roots", "np.polynomial", "np.linalg", "np.linalg.det", "np.linalg.slogdet", "scipy"],
        "reference": REFERENCE, "make_cases": CASES,
    }],
}
