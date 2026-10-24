"""第 7 章：对角化、矩阵指数与线性动力学。"""
THEORY = r"""
## 动机与直觉

对角化把耦合坐标改成独立模式：每个模式只乘一个数。矩阵幂描述离散迭代，矩阵指数描述连续演化。
二者都可沿同一组特征方向计算，但矩阵指数在没有足够特征向量时仍有定义。
本章证明可对角化判据、相似不变量与指数乘法条件，再把它们用于初值问题。

## 定义

**定义 7.1（相似与对角化）**：若 $B=P^{-1}AP$ 且 $P$ 可逆，称 $A,B$ 相似。
若 $B=\Lambda$ 为对角阵，称 $A$ 在指定数域上可对角化。
**定义 7.2（迹）**：$\operatorname{tr}A=\sum_i a_{ii}$；特征多项式继续使用 $p_A(z)=\det(zI-A)$。
**定义 7.3（矩阵指数）**：对有限方阵，定义
$$ e^{tA}=\sum_{k=0}^{\infty}\frac{t^kA^k}{k!},\qquad A^0=I.\tag{7.1} $$
本章动力学取实数时间 $t$、常系数 $A$。幂 $A^k$ 的基础范围为非负整数；负幂需要 $A$ 可逆。

## 核心推导

### 推导 7.1：可对角化当且仅当特征向量张满全空间

动机：特征值的数量按重数总是 $n$，却不保证能选出一套坐标轴；要准确找出缺少的条件。
若存在特征向量基 $v_1,\ldots,v_n$，令 $P=(v_1,\ldots,v_n)$，$\Lambda=\operatorname{diag}(\lambda_1,\ldots,\lambda_n)$。
按列计算，$AP$ 的第 $j$ 列为 $Av_j=\lambda_jv_j$，这也正是 $P\Lambda$ 第 $j$ 列。
因此 $AP=P\Lambda$。因为这些列是基，$P$ 可逆，左乘 $P^{-1}$ 得 $P^{-1}AP=\Lambda$。
若只假定所有特征向量张满，按第 1 章从张成组抽基，即得到这样的 $n$ 个向量。

反向，若 $P^{-1}AP=\Lambda$，左乘 $P$ 得 $AP=P\Lambda$。
逐列比较给 $A P_{:j}=\lambda_jP_{:j}$。可逆矩阵的列无关，特别每列非零，
所以每列都是特征向量，且 $n$ 列张成全空间。双向证明完成：
$$ A=P\Lambda P^{-1}\iff A\text{有一组特征向量基}.\tag{7.2} $$
第 5 章已证明不同本征空间是直和，故在特征多项式分裂的数域中，
可对角化又等价于 $\sum_\lambda\dim E_\lambda=n$。
代数重数总和为 $n$，每个几何重数不超过对应代数重数，所以这也等价于逐个重数相等。

对照 $J=\begin{pmatrix}\lambda&1\\0&\lambda\end{pmatrix}=\lambda I+N$。
$(J-\lambda I)v=0$ 给 $v_2=0$，所有特征向量只张成一条直线。
可以写两列相同方向的 $P$ 使形式上的 $JP=P\Lambda$ 成立，但 $P$ 不可逆，
证明中左乘 $P^{-1}$ 的那一步无法执行。这是缺陷矩阵的具体障碍。
Jordan 形式用广义特征向量补齐：$Ne_1=0,Ne_2=e_1$；第二个基向量的像含第一个方向，
所以留下超对角线的 $1$，不能变成对角阵。此处的二阶链已直接算出，足以说明本章对角化路径失效的机制。

条件与失效：在实数域可对角化要求实特征向量张满；平面旋转可在复数域对角化，却不能在实数域对角化。
重根本身不是障碍，缺少独立本征方向才是；例如 $\lambda I$ 的整个空间都是本征空间。

### 推导 7.2：逐个证明相似不变量

动机：换一套基会改变条目，却不应改变体积因子、特征值和模式总增长率的代数记录。
设 $B=P^{-1}AP$。
由第 4 章乘法公式与 $\det(P^{-1})\det P=1$，
$$ \det B=\det(P^{-1})\det A\det P=\det A. $$
对迹，先从定义证明循环交换：
$$ \operatorname{tr}(XY)=\sum_i\sum_jX_{ij}Y_{ji}
=\sum_j\sum_iY_{ji}X_{ij}=\operatorname{tr}(YX). $$
这里交换的是有限求和和标量乘法，未假设矩阵可交换。
取 $X=P^{-1}A,Y=P$，得 $\operatorname{tr}B=\operatorname{tr}(PP^{-1}A)=\operatorname{tr}A$。
最后 $zI-B=P^{-1}(zI-A)P$，再用乘法公式给
$$ p_B(z)=\det(zI-B)=\det(zI-A)=p_A(z).\tag{7.3} $$
所以根及代数重数保持不变。特别可对角化时 $\det A=\prod_j\lambda_j$，
$\operatorname{tr}A=\sum_j\lambda_j$，因为对角阵的行列式与迹分别为积与和。

条件与失效：相似要求两端使用同一个可逆换基。一般左右独立变换 $Q^{-1}AP$ 不保持这些量，
因为消去 $PP^{-1}$ 的步骤不存在。迹允许循环移动因子，但不允许任意交换中间因子。

### 推导 7.3：指数级数的收敛、求导与对易乘法公式

动机：先确保无穷级数合法，才可按模式计算并用它解微分方程。
使用最大绝对行和范数 $\|A\|_\infty=\max_i\sum_j|a_{ij}|$。
按乘法定义与三角不等式，
$\sum_j|(AB)_{ij}|\le\sum_k|a_{ik}|\sum_j|b_{kj}|\le\sum_k|a_{ik}|\|B\|_\infty$，
再取最大值得 $\|AB\|_\infty\le\|A\|_\infty\|B\|_\infty$。
归纳给 $\|A^k\|_\infty\le\|A\|_\infty^k$，所以式 (7.1) 各项范数之和不超过
$\sum_k |t|^k\|A\|_\infty^k/k!=e^{|t|\|A\|_\infty}$，级数绝对收敛。
在任意 $|t|\le T$，导数项的范数和不超过 $\|A\|_\infty e^{T\|A\|_\infty}$，
原级数与导数级数一致收敛，因此可逐项求导并重编号：
$$ \frac{d}{dt}e^{tA}=\sum_{k=1}^\infty\frac{t^{k-1}A^k}{(k-1)!}=Ae^{tA},\qquad e^{0A}=I.\tag{7.4} $$
同样 $e^{tA}A=Ae^{tA}$，因为每个 $A^k$ 都与 $A$ 交换。

现在设 $AB=BA$。归纳证明二项式：$k=0$ 时成立；若第 $k$ 次展开成立，乘以 $A+B$，
用 $B A^r=A^rB$ 合并相邻两项，系数按
$\binom{k}{j}+\binom{k}{j-1}=\binom{k+1}{j}$ 合成，得到下一次展开。
于是 $(A+B)^k=\sum_{r=0}^k\binom{k}{r}A^rB^{k-r}$。
两指数乘积的双级数绝对收敛，其范数和不超过 $e^{|t|\|A\|}e^{|t|\|B\|}$，
所以可按总次数 $k=r+s$ 分组：
$$ e^{tA}e^{tB}=\sum_{k=0}^\infty t^k\sum_{r=0}^k\frac{A^rB^{k-r}}{r!(k-r)!}
=\sum_{k=0}^\infty\frac{t^k(A+B)^k}{k!}=e^{t(A+B)}.\tag{7.5} $$
取 $B=-A$，得 $e^{tA}e^{-tA}=I$，所以矩阵指数总可逆，即使 $A$ 本身奇异。

若不交换，二项式的合项步骤失效：$(A+B)^2=A^2+AB+BA+B^2$，
而乘积指数的二阶交叉项只有 $AB$，两边相差 $t^2(AB-BA)/2+O(t^3)$。
具体取 $A=\begin{pmatrix}0&1\\0&0\end{pmatrix},B=\begin{pmatrix}0&0\\1&0\end{pmatrix}$。
$A^2=B^2=0$，故 $e^{tA}e^{tB}=(I+tA)(I+tB)=\begin{pmatrix}1+t^2&t\\t&1\end{pmatrix}$。
$(A+B)^2=I$，将级数分奇偶项给
$e^{t(A+B)}=\begin{pmatrix}\cosh t&\sinh t\\\sinh t&\cosh t\end{pmatrix}$；
在 $t=1$ 两个对角元分别为 $(2,1)$ 和 $(\cosh1,\cosh1)$，两矩阵不同。

条件与失效：式 (7.5) 的充分条件是 $AB=BA$；非交换时不能把指数拆成两项乘积。
有限矩阵保证范数估计和级数操作；无界算符还需定义域与收敛性分析，不能直接移植此证明。

### 推导 7.4：对角化计算幂、指数并解常系数初值问题

动机：把指数的无限次矩阵乘法变成每个模式上的一个标量函数。
若 $A=P\Lambda P^{-1}$，由结合律逐个消去相邻的 $P^{-1}P$，
$A^k=P\Lambda^kP^{-1}$；$k=0$ 用 $A^0=I$ 同样成立。
将它代入绝对收敛的式 (7.1)，固定矩阵乘法可与极限交换，得到
$$ A^k=P\Lambda^kP^{-1},\qquad e^{tA}=P\operatorname{diag}(e^{t\lambda_j})P^{-1}.\tag{7.6} $$
为解 $x'=Ax,x(0)=x_0$，先解基坐标 $Pc=x_0$；再设 $y=P^{-1}x$。
$P$ 常数，求导给 $y'=P^{-1}Ax=\Lambda y$，即 $y_j'=\lambda_jy_j$。
标量积分因子给 $\frac{d}{dt}(e^{-\lambda_jt}y_j)=0$，所以 $y_j(t)=c_je^{\lambda_jt}$。
换回原坐标，
$$ x(t)=\sum_jc_je^{\lambda_jt}v_j=e^{tA}x_0.\tag{7.7} $$
一般 $A$ 无法对角化时，式 (7.4) 仍直接证明 $e^{tA}x_0$ 满足方程与初值。
若 $x(t)$ 是任意另一个解，乘积分法则给
$\frac{d}{dt}(e^{-tA}x)=-Ae^{-tA}x+e^{-tA}Ax=0$，
因为 $A$ 与其指数交换。于是 $e^{-tA}x=x_0$，再乘 $e^{tA}$，仍得式 (7.7)，证明唯一性。

**完整手算算例**：$A=\begin{pmatrix}2&1\\0&1\end{pmatrix},x_0=(0,1)^T$。
$p_A(z)=(z-2)(z-1)$；解 $(A-2I)v=0$ 得 $v_1=(1,0)^T$，
解 $(A-I)v=0$ 得 $v_2=(-1,1)^T$。
所以 $P=\begin{pmatrix}1&-1\\0&1\end{pmatrix}$，$P^{-1}=\begin{pmatrix}1&1\\0&1\end{pmatrix}$，
$\Lambda=\operatorname{diag}(2,1)$。逐列 $AP=(2v_1,v_2)=P\Lambda$，且 $\det P=1$。
解 $Pc=x_0$：第二行 $c_2=1$，第一行 $c_1-c_2=0$，故 $c_1=1$。
因此
$$ x(t)=e^{2t}(1,0)^T+e^t(-1,1)^T=(e^{2t}-e^t,e^t)^T. $$
代 $t=0$ 得 $(0,1)^T$；求导为 $(2e^{2t}-e^t,e^t)^T$，
而 $Ax(t)=(2(e^{2t}-e^t)+e^t,e^t)^T$，与导数逐项相等。
同一分解还给 $A^kx_0=(2^k-1,1)^T$，如 $k=3$ 时为 $(7,1)^T$。

条件与失效：对角化路径要求可逆特征向量矩阵；缺陷 Jordan 块不满足。
但对 $J=\lambda I+N,N^2=0$，对易乘法给 $e^{tJ}=e^{\lambda t}(I+tN)$，
额外的 $t$ 因子记录广义特征方向的耦合，说明失败的是纯独立指数模式的展开，不是指数定义本身。
时间依赖 $A(t)$ 时求导与常数换基的步骤需重做；不同时刻不对易时不能直接写 $e^{\int A(t)dt}$。

## 适用条件与失效情形

相似与对角化在指定数域内讨论；矩阵指数却对任意有限实或复方阵存在。
可对角化但 $P$ 近奇异时，基坐标计算会放大浮点误差；谱公式精确成立不代表其数值实现总稳定。
本章练习只实现两个互异特征值的二阶路径，避免把二次根公式冒充一般高阶特征算法。

## 几何直觉

换到本征基后，离散第 $j$ 模式乘 $\lambda_j^k$，连续模式乘 $e^{t\lambda_j}$。
模小于一的离散模式衰减；连续模式的模为 $e^{t\operatorname{Re}\lambda_j}$，衰减判据改成实部为负。
`matrix_power` 用迭代次数滑杆显示同一初始点的离散轨迹，能观察主方向的选择以及增长或收缩。
不能把离散的 $|\lambda|<1$ 判据直接当作微分方程的 $\operatorname{Re}\lambda<0$。

## 与量子力学的连接

对时间无关 Hermitian Hamilton 矩阵 $H$，Schrödinger 方程
$i\hbar\,d|\psi\rangle/dt=H|\psi\rangle$ 等价于 $\psi'=(-iH/\hbar)\psi$，故
$$ U(t)=e^{-iHt/\hbar},\qquad |\psi(t)\rangle=U(t)|\psi(0)\rangle. $$
在正交归一本征基 $H|n\rangle=E_n|n\rangle$ 中，
$U(t)=\sum_n e^{-iE_nt/\hbar}|n\rangle\langle n|$；对角化 H 就是寻找定态与其能量。
对实数 $t,\hbar\ne0$，逐项取伴随给 $U^\dagger=e^{iHt/\hbar}$，
再由对易指数公式得 $U^\dagger U=I$，所以总概率保持。
Heisenberg 绘景把演化移到可观测量：$O_H(t)=U^\dagger O_S U$。
当 $O_S$ 无显式时间依赖，乘积求导给
$\dot O_H=(i/\hbar)U^\dagger(HO_S-O_SH)U=(i/\hbar)[H,O_H]$，
最后一步用了 H 与 U 交换。例如 $H=\hbar\omega\sigma_z/2$ 与 $O_S=\sigma_x$ 不对易，
所以自旋横向分量随时间演化；可逐项算出其演化：令 $\theta=\omega t$，
$U=\operatorname{diag}(e^{-i\theta/2},e^{i\theta/2})$，
$$ U^\dagger\sigma_xU=\begin{pmatrix}0&e^{i\theta}\\e^{-i\theta}&0\end{pmatrix}
=\cos\theta\,\sigma_x-\sin\theta\,\sigma_y,\qquad
\sigma_y=\begin{pmatrix}0&-i\\i&0\end{pmatrix}. $$
初态 $|+x\rangle=(1,1)^T/\sqrt2$ 有 $\langle\sigma_x(t)\rangle=\cos\theta$，
因为夹乘上式给 $(e^{i\theta}+e^{-i\theta})/2$。
相应的 +x 投影为 $(I+\sigma_x)/2$，概率是 $(1+\cos\theta)/2$；
$\theta=\pi/2$ 时为 1/2，$\theta=\pi$ 时为 0，定态相位的相对变化已经成为具体测量概率。
若 $[H,O_S]=0$ 则该量守恒。
"""

REFERENCE = '''def spectral_action(A, x, k, t):
    A = np.asarray(A, dtype=complex)
    x = np.asarray(x, dtype=complex)
    if A.shape != (2, 2) or x.shape != (2,):
        raise ValueError("需要二阶矩阵及二维向量")
    trace = A[0, 0] + A[1, 1]
    determinant = A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]
    # 与 trace**2-4*determinant 恒等，减少接近标量阵时的相消。
    root = np.sqrt((A[0, 0]-A[1, 1])**2 + 4*A[0, 1]*A[1, 0])
    values = np.array([(trace+root)/2, (trace-root)/2])
    if abs(values[0]-values[1]) <= 1e-10:
        raise ValueError("本路径要求两个互异特征值")
    P = np.zeros((2, 2), dtype=complex)
    for j in range(2):
        M = A - values[j]*np.eye(2)
        i = int(np.argmax(np.sum(np.abs(M)**2, axis=1)))
        v = np.array([-M[i, 1], M[i, 0]])
        P[:, j] = v / np.sqrt(np.sum(np.abs(v)**2))
    d = P[0, 0]*P[1, 1] - P[0, 1]*P[1, 0]
    c = np.array([(x[0]*P[1, 1]-P[0, 1]*x[1])/d,
                  (P[0, 0]*x[1]-x[0]*P[1, 0])/d])
    return P @ ((values**k)*c), P @ (np.exp(t*values)*c)
'''
CASES = '''def make_cases(rng):
    cases = [(np.array([[2., 1.], [0., 1.]]), np.array([0., 1.]), 3, .5),
             (np.array([[0., -1.], [1., 0.]]), np.array([1., 0.]), 4, 1.),
             (np.diag([0., 2.]), np.array([1., -2.]), 0, 0.),
             (np.diag([0., -2.]), np.array([1., 2.]), 3, -.2),
             (np.diag([1., 1.+1e-8]), np.array([1., 1.]), 2, .5)]
    for index in range(28):
        angle = rng.uniform(-3., 3.)
        Q = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        P = Q @ np.array([[1., .4], [0., 1.]])
        Pinv = np.array([[1., -.4], [0., 1.]]) @ Q.T
        values = np.array([rng.uniform(-2., -.5), rng.uniform(.5, 2.)])
        if index % 2:
            values = values + 1j*np.array([.4, -.3])
        A = P @ np.diag(values) @ Pinv
        x = rng.normal(size=2) + (1j*rng.normal(size=2) if index % 2 else 0)
        cases.append((A, x, int(rng.integers(0, 9)), float(rng.uniform(-1., 1.))))
    return cases
'''
WRITTEN = [
    {"id": "ch7-w1", "level": "基础", "problem_md": r"对 $A=\begin{pmatrix}2&1\\0&1\end{pmatrix}$ 构造 P、Λ，并说明可对角化。", "solution_md": r"""
1. 依据三角特征多项式，特征值为 $2,1$。
2. 依据核方程，分别取特征向量 $v_1=(1,0)^T,v_2=(-1,1)^T$。
3. 依据列组装，$P=\begin{pmatrix}1&-1\\0&1\end{pmatrix}$，$\det P=1\ne0$，所以两列成基。
4. 依据逐列特征方程，$AP=P\operatorname{diag}(2,1)$，左乘 $P^{-1}$ 得到对角化。
""", "pitfalls_md": "Λ 的特征值排列必须与 P 的列一一对应。"},
    {"id": "ch7-w2", "level": "基础", "problem_md": r"$A=\operatorname{diag}(2,-1)$，求 $A^3$、$e^{tA}$ 及 $x'=Ax,x(0)=(1,2)^T$ 的解。", "solution_md": r"""
1. 依据对角矩阵乘法，$A^3=\operatorname{diag}(8,-1)$。
2. 依据指数级数逐对角条目求和，$e^{tA}=\operatorname{diag}(e^{2t},e^{-t})$。
3. 依据式 (7.7)，$x(t)=(e^{2t},2e^{-t})^T$。
4. 依据求导与代入，$x'=(2e^{2t},-2e^{-t})^T=Ax$，且 $x(0)=(1,2)^T$，核对初值问题。
""", "pitfalls_md": "e^{tA} 不是逐元素取指数；非对角位置为零，指数后仍为零。"},
    {"id": "ch7-w3", "level": "进阶", "problem_md": r"$B=P^{-1}AP$。从条目求和证明迹不变，并说明为何这不要求 AP=PA。", "solution_md": r"""
1. 依据迹与矩阵乘法定义，$\operatorname{tr}(XY)=\sum_{i,j}X_{ij}Y_{ji}$。
2. 依据标量交换律与有限求和换序，它等于 $\sum_{j,i}Y_{ji}X_{ij}=\operatorname{tr}(YX)$。
3. 依据这一恒等式取 $X=P^{-1}A,Y=P$，有 $\operatorname{tr}B=\operatorname{tr}(PP^{-1}A)=\operatorname{tr}A$。
4. 依据证明中的运算对象，交换的是标量条目的顺序，未交换矩阵 A 与 P，因此不需要 AP=PA。
""", "pitfalls_md": "tr(XYZ)=tr(ZXY) 是循环关系，不能据此断言 tr(XYZ)=tr(XZY)。"},
    {"id": "ch7-w4", "level": "进阶", "problem_md": r"$J=\begin{pmatrix}1&1\\0&1\end{pmatrix}$。证明不可对角化，求 $e^{tJ}$。", "solution_md": r"""
1. 依据特征多项式与核方程，唯一特征值为 1，$(J-I)v=0$ 给 $v_2=0$，本征空间只有一维。
2. 依据式 (7.2)，二维空间没有特征向量基，故不可对角化。
3. 依据 $N=J-I$ 的直接乘法，$N^2=0$，指数级数截断成 $e^{tN}=I+tN$。
4. 依据 I 与 N 对易及式 (7.5)，$e^{tJ}=e^t(I+tN)=e^t\begin{pmatrix}1&t\\0&1\end{pmatrix}$，其中 t 项来自广义特征方向。
""", "pitfalls_md": "不可对角化不表示矩阵指数不存在；缺陷块产生多项式乘指数。"},
    {"id": "ch7-w5", "level": "考试风格", "problem_md": r"解 $x'=2x+y,y'=y$，初值 $x(0)=0,y(0)=1$，用对角化并核对结果。", "solution_md": r"""
1. 依据本征方程，$A=\begin{pmatrix}2&1\\0&1\end{pmatrix}$ 的特征向量取 $(1,0)^T,(-1,1)^T$，特征值为 $2,1$。
2. 依据初值的基展开，$(0,1)^T=c_1(1,0)^T+c_2(-1,1)^T$，第二分量给 $c_2=1$，第一分量给 $c_1=1$。
3. 依据独立模式解式 (7.7)，$(x,y)^T=e^{2t}(1,0)^T+e^t(-1,1)^T$，即 $x=e^{2t}-e^t,y=e^t$。
4. 依据逐项求导，$x'=2e^{2t}-e^t=2x+y$，$y'=e^t=y$；代 t=0 得初值，故解已核对。
""", "pitfalls_md": "初值分量不是本征模式系数，先解 Pc=x₀，不能直接把 (0,1) 当 c。"},
    {"id": "ch7-w6", "level": "考试风格", "problem_md": r"$H=\hbar\omega\sigma_z/2$，初态 $(|0\rangle+|1\rangle)/\sqrt2$，求演化态与 $\sigma_x$ 的期望值。", "solution_md": r"""
1. 依据对角 Hamilton 的本征值，$E_0=\hbar\omega/2,E_1=-\hbar\omega/2$。
2. 依据 $U(t)=e^{-iHt/\hbar}$，态变为 $(e^{-i\omega t/2}|0\rangle+e^{i\omega t/2}|1\rangle)/\sqrt2$。
3. 依据 $\sigma_x|0\rangle=|1\rangle,\sigma_x|1\rangle=|0\rangle$，期望值是两交叉振幅之和 $(e^{i\omega t}+e^{-i\omega t})/2$。
4. 依据欧拉公式，该值为 $\cos(\omega t)$；两分量模平方仍各为 $1/2$，所以归一化保持。
""", "pitfalls_md": "相对相位影响干涉期望值；能量本征概率不变不意味着所有可观测量都不变。"},
]
QUIZ = [
    {"q_md": "可对角化的充要条件是？", "options": ["行列式非零", "特征向量张满全空间", "迹非零", "特征值全相同"], "answer": 1, "explain_md": "将特征向量基排成可逆 P，逐列 AP=PΛ；反向由相似对角化读出这套基。"},
    {"q_md": "以下哪项是相似不变量？", "options": ["每个条目", "每列模长", "特征多项式", "最大条目"], "answer": 2, "explain_md": "zI-B=P⁻¹(zI-A)P，行列式乘法公式消去两侧换基因子。"},
    {"q_md": r"何时可保证 $e^{t(A+B)}=e^{tA}e^{tB}$？", "options": [r"$AB=BA$", "A、B 都可逆即可", "A、B 同阶即可", "tr A=tr B 即可"], "answer": 0, "explain_md": "对易使二项式合并合法，再由绝对收敛双级数按总次数求和。"},
    {"q_md": r"$N^2=0$ 时 $e^{tN}$ 为？", "options": ["0", r"$tN$", r"$I$", r"$I+tN$"], "answer": 3, "explain_md": "指数级数的零次项为 I，一次项 tN，从二次起均为零。"},
    {"q_md": "时间无关 Hamilton 的本征态随时间如何变化？", "options": ["能量每次翻倍", "乘以能量决定的相位", "总概率指数衰减", "立刻变为零态"], "answer": 1, "explain_md": "每个本征态乘 e^{-iEt/ħ}，对实能量模为 1，保持概率。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch7-q{_i}"

CHAPTER = {
    "id": "ch7", "num": 7, "title": "对角化与矩阵指数", "depends_on": ["ch4", "ch5", "ch6"],
    "theory_md": THEORY, "viz": ["matrix_power"], "written": WRITTEN, "quiz": QUIZ,
    "code_exercises": [{
        "id": "ch7-c1", "title": "从二次特征方程实现谱分解作用",
        "prompt_md": r"实现 `spectral_action(A, x, k, t)`，返回 `(A^k @ x, exp(t*A) @ x)` 两个形状 (2,) 的数组。第二项 exp 表示矩阵指数。A 为实或复二阶矩阵，保证两个特征值互异、间隔>$10^{-10}$、本征基条件适中，k 为非负整数，t 为实数。必须走对角化路径：手算 tr、det，用二次公式求复特征根，判别式可用恒等形式 $(a-d)^2+4bc$ 计算以减少近标量阵的相消。从 (A−λI)v=0 构造两列 P，再用二阶消元或克拉默公式解 Pc=x，分别缩放模式后乘回 P。可选 A−λI 中模长较大的一行 (a,b)，取核向量 (−b,a) 并归一化。允许 np.sqrt 与作用于特征值的 np.exp，禁止 linalg（含 eig/eigvals/eigh/eigvalsh/matrix_power/solve/inv）、np.poly、np.roots、np.polynomial、scipy；禁用 numpy.matrix 系列封装；不得用循环乘 A 代替指定谱路径。重根与缺陷矩阵不在此函数的输入域，遇重根应抛出 ValueError。",
        "function_name": "spectral_action", "starter": "def spectral_action(A, x, k, t):\n    # TODO：二次根公式 -> 核向量 -> 解 Pc=x -> 缩放本征模式\n    return np.zeros(2, dtype=complex), np.zeros(2, dtype=complex)\n",
        "banned": ["np.matrix", "np.asmatrix", "np.mat", "np.bmat", "np.matrixlib", "np.matlib", "np.matrix.I", "np.matrix.getI", "np.poly", "np.linalg", "np.linalg.eig", "np.linalg.eigvals", "np.linalg.eigh", "np.linalg.eigvalsh", "np.linalg.matrix_power", "np.linalg.solve", "np.linalg.inv", "np.roots", "np.polynomial", "scipy"],
        "reference": REFERENCE, "make_cases": CASES,
    }],
}
