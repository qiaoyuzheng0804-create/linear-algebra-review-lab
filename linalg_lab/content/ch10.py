"""第 10 章：张量积的构造、Kronecker 表示、纠缠与部分迹。"""
THEORY = r"""
## 动机与直觉

两个系统各自能取多种基态，联合系统必须记录每一种“第一个状态与第二个状态”的配对。
系数应对两端分别线性：改变第一个态的叠加再配第二个态，应与先分别配对再相加相同。
这就是双线性的动机。张量积提供一个空间，使任何这样的双线性规则都能通过其上的线性映射表达，
称为泛性质。本章直接用基构造这一空间，再把联合态的相关性落实到可算的矩阵。

## 定义

固定同一数域 $\mathbb F$。设 V、W 维数分别为 m、n，有序基为 $e_1,\ldots,e_m$ 与 $f_1,\ldots,f_n$。
**定义 10.1（构造性张量积）**：建立含 mn 个形式符号 $g_{ij}$ 的自由向量空间 T，
其元素是唯一的有限系数表 $\sum_{ij}c_{ij}g_{ij}$，加法与数乘按系数进行。
若 $v=\sum_i a_ie_i,w=\sum_j b_jf_j$，规定
$$ v\otimes w=\sum_{i,j}a_ib_jg_{ij},\qquad e_i\otimes f_j=g_{ij}.\tag{10.1} $$
T 记作 $V\otimes W$。零维空间的张量积为零空间。
**定义 10.2（联合基与直积态）**：量子系统取复内积空间，基态 $|ij\rangle=|i\rangle_A\otimes|j\rangle_B$，
按 i 为外层、j 为内层排序。非零纯态若能写成 $|a\rangle\otimes|b\rangle$ 称直积态，
否则称纠缠纯态；各局域态可归一化并将剩余相位移入其中一端。
**定义 10.3（密度算符与约化态）**：归一化纯态的密度算符为
$\rho=|\psi\rangle\langle\psi|$；一般密度矩阵满足 $\rho=\rho^\dagger\ge0,\operatorname{tr}\rho=1$。
只观察 A 时使用 $\rho_A=\operatorname{Tr}_B\rho$。部分迹不是把整个矩阵求迹后再分配，
它要保留 A 的行列指标，只对 B 的同一指标求和。

## 核心推导

### 推导 10.1：由基配对得到双线性与维数乘积

动机：先证明定义确实给出所需的空间与运算，避免把维数乘积当作凭空规定的物理口号。
第一步构造 T 的 mn 个坐标槽位，每对 (i,j) 恰对应一个槽位。
第二步令 gᵢⱼ 在自己的槽位取 1、其余取 0；任意系数表等于这些 g 的线性组合，故它们张成 T。
第三步若 $\sum c_{ij}g_{ij}=0$，逐槽位比较给每个 $c_{ij}=0$，故它们线性无关。
所以这是一组基，向量空间的维数定义给
$$ \dim(V\otimes W)=\dim V\,\dim W=mn.\tag{10.2} $$
这一步构造了所有系数表，包含不能写成单个 $a_ib_j$ 的表，而不只收集直积态。

按式 (10.1) 展开，若 v 的坐标为 a、v′ 为 a′，
$(sv+tv')\otimes w=\sum_{ij}(sa_i+ta_i')b_jg_{ij}
=s(v\otimes w)+t(v'\otimes w)$。
同样展开第二端坐标得
$v\otimes(sw+tw')=s(v\otimes w)+t(v\otimes w')$，所以两端分别线性。
但它不是对有序对 (v,w) 的联合线性：$(sv)\otimes(sw)=s^2(v\otimes w)$。

还可检查泛性质解释构造不依赖任意的符号命名。
若 $F:V\times W\to X$ 双线性，先规定线性映射 $L:T\to X$ 的基向量像为
$L(g_{ij})=F(e_i,f_j)$，再线性延拓。
则 $L(v\otimes w)=\sum_{ij}a_ib_jF(e_i,f_j)=F(v,w)$，最后等号逐次使用 F 的双线性。
任何满足该式的 L 在所有 $e_i\otimes f_j$ 上都已被指定，故唯一。
若换一对基作出另一个空间 T′，把同一双线性配对分别通过 T、T′ 分解，
得到互相逆的线性映射：它们的复合在全部直积基上为恒等，从而在全空间为恒等。
所以用基构造选定了坐标表示，而没有使抽象张量积依赖某个特定基。

条件与失效：维数乘积式讨论有限维、同数域空间。直和 $V\oplus W$ 记录一对独立坐标，
基数是 m+n；张量积记录两端基的所有配对，基数是 mn，来自不同构造步骤。
若只允许单个 $v\otimes w$，这些系数表通常对加法不封闭，不能充当整个 T。

### 推导 10.2：算符张量积、Kronecker 分块与复合规则

动机：局域操作 A、B 应分别作用在两个粒子上，并可在联合坐标中执行。
设 $A:V\to V'$、$B:W\to W'$。在每个直积基上规定
$(A\otimes B)(e_i\otimes f_j)=Ae_i\otimes Bf_j$，再线性延拓到整个空间。
对任意直积向量展开式 (10.1)，其像是 $Av\otimes Bw$，与要求一致。
若 $Ae_i=\sum_k A_{ki}e'_k,Bf_j=\sum_l B_{lj}f'_l$，双线性给
$$ (A\otimes B)(e_i\otimes f_j)=\sum_{k,l}A_{ki}B_{lj}e'_k\otimes f'_l. $$
因此联合矩阵的条目为 $(A\otimes B)_{(k,l),(i,j)}=A_{ki}B_{lj}$。
若 B 有 p 行 q 列，按约定的扁平排序就是
$(A\otimes B)_{kp+l,\,iq+j}=A_{ki}B_{lj}$（此式下标从零开始）。
每个 A 条目乘整块 B，这正是 Kronecker 积的分块定义。

接着证明乘法规则。设 C、D 的输出分别是 A、B 的输入，取定义域基 $e_i,f_j$，
$$ (A\otimes B)(C\otimes D)(e_i\otimes f_j)
=(A\otimes B)(Ce_i\otimes Df_j)
=ACe_i\otimes BDf_j
=(AC\otimes BD)(e_i\otimes f_j). $$
第二个等号可以展开 $Ce_i=\sum_k c_ke_k,Df_j=\sum_l d_lf_l$，
两边都是 $\sum_{kl}c_kd_l(Ae_k\otimes Bf_l)$，所以未假设输入一定是某个单独基向量。
两线性算符在全部直积基上一致，对任意 $\sum z_{ij}e_i\otimes f_j$ 逐项作用仍一致，故
$$ (A\otimes B)(C\otimes D)=AC\otimes BD.\tag{10.3} $$

量子空间的联合内积令直积正交归一基满足 $\langle ij|kl\rangle=\delta_{ik}\delta_{jl}$，
双线性展开给
$\langle a\otimes b,c\otimes d\rangle=\langle a,c\rangle\langle b,d\rangle$。
特别 $\|a\otimes b\|=\|a\|\|b\|$。
按上面矩阵元取共轭并交换联合行列，得到
$\overline{A_{ki}B_{lj}}=\bar A_{ki}\bar B_{lj}$，对应
$$ (A\otimes B)^\dagger=A^\dagger\otimes B^\dagger.\tag{10.4} $$
式 (10.3)、(10.4) 进一步给局域酉算符的张量积也酉：
$(U\otimes V)^\dagger(U\otimes V)=U^\dagger U\otimes V^\dagger V=I\otimes I$。

**手算算例：局域翻转。** 取 $X=\begin{pmatrix}0&1\\1&0\end{pmatrix}$，
$X\otimes I=\begin{pmatrix}0&0&1&0\\0&0&0&1\\1&0&0&0\\0&1&0&0\end{pmatrix}$，
其列序是 00、01、10、11。因此 $(X\otimes I)|01\rangle=|11\rangle$；
作用在 $(|00\rangle+|11\rangle)/\sqrt2$ 上，像为 $(|10\rangle+|01\rangle)/\sqrt2$，范数平方仍为 1。

条件与失效：复合规则要求中间两端维数分别相容，且必须按固定的基排序解释条目。
通常 $A\otimes B\ne B\otimes A$，交换两系统要另加交换算符 S，
$S(a\otimes b)=b\otimes a$，不能在同一坐标顺序中直接交换因子。

### 推导 10.3：Bell 态不能分离成两个单比特态

动机：叠加态有时仍是直积态；需要用系数检验证明真正的纠缠。
令 $|\Phi^+\rangle=(|00\rangle+|11\rangle)/\sqrt2$。
假设它能写成 $(a|0\rangle+b|1\rangle)\otimes(c|0\rangle+d|1\rangle)$，其中系数可为复数。
双线性逐项展开得到
$ac|00\rangle+ad|01\rangle+bc|10\rangle+bd|11\rangle$。
四个联合基无关，比较系数必须同时满足
$$ ac=1/\sqrt2,\quad ad=0,\quad bc=0,\quad bd=1/\sqrt2. $$
由 ac 非零知 a、c 都非零；ad=0 与 a≠0 迫使 d=0；
但 bd=1/√2 又要求 d 非零，矛盾，所以不存在这四个系数。
也可把第一、四式相乘得 abcd=1/2，第二、三式相乘得 abcd=0，直接矛盾。
这在复数域仍成立，因为非零复数相乘不会等于零。

条件与失效：判据针对非零纯态。$(|00\rangle+|01\rangle)/\sqrt2$
可以提出 $|0\rangle$，等于 $|0\rangle\otimes(|0\rangle+|1\rangle)/\sqrt2$，
因此有叠加不代表有纠缠。上述比较系数的方法依靠联合基无关，不是依靠态的名称。

### 推导 10.4：密度矩阵、部分迹与完整数字算例

动机：只操作 A 时，所有测量预测应由 A 空间上的一个矩阵给出，不需要保存 B 的可观测信息。
纯态 $\rho=\psi\psi^\dagger$ 满足 $\rho^\dagger=\rho$，
$x^\dagger\rho x=|\psi^\dagger x|^2\ge0$，且
$\operatorname{tr}\rho=\sum_i|\psi_i|^2=1$。
若 $\rho=\sum_a p_a\psi_a\psi_a^\dagger$，$p_a\ge0,\sum p_a=1$，
三条性质逐项求和保留，表示对纯态制备的不确定性。

为了保持任意局域观测量 O 的期望值，要求
$\operatorname{tr}((O\otimes I)\rho)=\operatorname{tr}(O\rho_A)$。
按联合指标展开左边：
$\sum_{i,j,k} O_{ik}\rho_{(k,j),(i,j)}$。
右边为 $\sum_{i,k}O_{ik}(\rho_A)_{ki}$。对任意 O 相等要求
$$ (\rho_A)_{ik}=\sum_j\rho_{(i,j),(k,j)},\quad
\operatorname{Tr}_B\rho=\sum_j(I\otimes\langle j|)\rho(I\otimes|j\rangle).\tag{10.5} $$
前者直接给手写算法，后者说明它在 B 基矢上夹取再求和。
若 B 换正交归一基 $|j'\rangle=\sum_l W_{lj}|l\rangle$，
求和系数变成 $\sum_j\bar W_{lj}W_{sj}=\delta_{ls}$，
所以仍只保留 l=s 的项，证明部分迹不依赖 B 端正交基选择。
约化矩阵依然正：对任意 x，
$x^\dagger\rho_Ax=\sum_j(x\otimes j)^\dagger\rho(x\otimes j)\ge0$；
其迹为 $\sum_{ij}\rho_{(i,j),(i,j)}=1$，自伴性也由条目共轭得到。

**完整手算算例：不等权且含相位的二比特态。**
取 $|\psi\rangle=(\sqrt3|00\rangle+i|11\rangle)/2$，按 00、01、10、11 排列，
列向量为 $(\sqrt3/2,0,0,i/2)^T$，于是
$$ \rho=\begin{pmatrix}
3/4&0&0&-i\sqrt3/4\\0&0&0&0\\0&0&0&0\\i\sqrt3/4&0&0&1/4
\end{pmatrix}. $$
上右角是 $(\sqrt3/2)\overline{i/2}=-i\sqrt3/4$，下左角是它的共轭。
对 B 求迹，逐个条目算：
$(\rho_A)_{00}=\rho_{00,00}+\rho_{01,01}=3/4+0=3/4$；
$(\rho_A)_{01}=\rho_{00,10}+\rho_{01,11}=0+0=0$；
$(\rho_A)_{10}=\rho_{10,00}+\rho_{11,01}=0+0=0$；
$(\rho_A)_{11}=\rho_{10,10}+\rho_{11,11}=0+1/4=1/4$。
所以 $\rho_A=\operatorname{diag}(3/4,1/4)$。同样对 A 求迹得到 $\rho_B=\operatorname{diag}(3/4,1/4)$。
$\rho_A^2=\operatorname{diag}(9/16,1/16)\ne\rho_A$，纯度 $\operatorname{tr}\rho_A^2=5/8$。
局域 $\sigma_z$ 的期望值为 $3/4-1/4=1/2$，局域 $\sigma_x$ 的期望值为 0；
直接在联合空间计算 $\langle\psi|\sigma_z\otimes I|\psi\rangle=3/4-1/4$，结果一致。
在 Bell 态的等权情形，四个部分迹条目变成 $1/2,0,0,1/2$，即 $\rho_A=I_2/2$，纯度为 1/2。

条件与失效：部分迹中只求同一个被消去指标 j 的对角配对，不是累加整块所有条目。
全态的非对角元可能含纠缠相关性，局域矩阵不保留所有这些信息；
例如 Bell 纯态与混合态 $(|00\rangle\langle00|+|11\rangle\langle11|)/2$ 有相同约化态，
所以只知道约化矩阵不能唯一恢复联合态。

### 推导 10.5：Schmidt 分解与约化态判据

动机：把不可分离性转成一个与局域坐标无关、能由 SVD 数出来的秩。
令纯态的系数矩阵为 C，$C_{ij}$ 对应 $|i\rangle_A|j\rangle_B$。
第 9 章 SVD 给 $C=U\Sigma V^\dagger$，逐条目展开并按有限求和换序：
$$ |\psi\rangle=\sum_{ij\alpha}s_\alpha U_{i\alpha}\bar V_{j\alpha}|i\rangle|j\rangle
=\sum_{\alpha=1}^r s_\alpha|u_\alpha\rangle\otimes|w_\alpha\rangle,
\quad w_\alpha=\bar v_\alpha.\tag{10.6} $$
U 的列正交归一；$w_\alpha^\dagger w_\beta=\overline{v_\alpha^\dagger v_\beta}=\delta_{\alpha\beta}$，
所以两端基都正交归一。所有 $s_\alpha>0$，归一化给 $\sum s_\alpha^2=1$，r 称 Schmidt 秩。
这就是 Schmidt 分解；不是对整个 mn 维态向量直接作特征分解，而是先排成 m×n 振幅矩阵再作 SVD。

将该展开代入纯态密度矩阵：
$\rho=\sum_{\alpha\beta}s_\alpha s_\beta
|u_\alpha\rangle\langle u_\beta|\otimes|w_\alpha\rangle\langle w_\beta|$。
部分迹对后因子给
$\sum_j\langle j|w_\alpha\rangle\langle w_\beta|j\rangle
=\langle w_\beta|w_\alpha\rangle=\delta_{\alpha\beta}$，
交叉项因此消去，得到
$$ \rho_A=\sum_\alpha s_\alpha^2|u_\alpha\rangle\langle u_\alpha|=CC^\dagger,
\quad\rho_B=\sum_\alpha s_\alpha^2|w_\alpha\rangle\langle w_\alpha|=C^T\bar C.\tag{10.7} $$
两端非零本征值同为 $s_\alpha^2$，但维数不同时零值数量可不同。
B 端是 $C^T\bar C$ 而非一般情况下的 $C^\dagger C$；二者互为共轭，有相同实谱。

现在证明判据。若 r=1，归一化使 s₁=1，式 (10.6) 就是一个直积态。
反向若 $\psi=a\otimes b$，则 $C_{ij}=a_ib_j$，每列是 a 的倍数，非零 C 的秩为 1，故 Schmidt 秩为 1。
约化态纯度为 $\operatorname{tr}(\rho_A^2)=\sum_\alpha s_\alpha^4$。
利用 $\sum s_\alpha^2=1$，展开平方得
$1=\sum s_\alpha^4+2\sum_{\alpha<\beta}s_\alpha^2s_\beta^2$。
若 r>1，至少一个交叉乘积严格为正，故纯度严格小于 1；若 r=1，纯度等于 1，且 $\rho_A^2=\rho_A$。
因此对联合纯态，约化态纯当且仅当联合态为直积；约化态混合当且仅当联合纯态纠缠。
上例的两个 Schmidt 系数为 $\sqrt3/2,1/2$，相位 i 吸收进一个局域基，r=2，判定纠缠。

条件与失效：这个判据要求联合态已知为纯态。
上述等权经典混合态是两个直积态的概率混合，属于可分混合态，却也给 $I/2$ 的混合约化态；
因此“约化态混合⇒纠缠”在联合混合态上失效，断点在于无法用单个振幅矩阵 C 的 Schmidt 展开表示整个 ρ。

## 适用条件与失效情形

本章区分三个对象：单个直积向量、直积向量的线性叠加、纯态密度矩阵的概率混合。
叠加具有相位与交叉项；概率混合按正权重加密度矩阵。它们的联合测量预测可不同。
部分迹对任意联合密度矩阵都定义良好，但用局域纯度判纠缠只适用于已知联合纯态。
代码中的扁平索引始终按 i·dimB+j，消去 A 和消去 B 的求和式必须分别实现。

## 几何直觉

联合空间的基方向数相乘，N 个比特的完整态空间维数是 $2^N$：
每加入一个二态系统，所有旧基都分别配上 0 与 1，方向数翻倍。
直积态的振幅表被约束为秩一 $ab^T$，只是整个联合空间中的一部分，通常不构成子空间；
相加两个这样的表可能变成秩二。直积是两个系统的配对操作，叠加是在已经构造好的联合空间中作向量加法。
`tensor_dims` 用可调基配对网格显示 mn 个方向；这不是把空间画成 m+n 个局域箭头。

## 与量子力学的连接

两个可区分自旋的联合态属于 $\mathbb C^2\otimes\mathbb C^2$。
只测第一个粒子的 $S_z$，算符是 $(\hbar\sigma_z/2)\otimes I$；
同时测两个粒子的 z 自旋相关性，则使用 $\sigma_z\otimes\sigma_z$。
Bell 态满足 $(\sigma_z\otimes\sigma_z)|\Phi^+\rangle=|\Phi^+\rangle$，
因为 00 与 11 两个分量的特征值乘积都为 +1，故相关期望为 1；
两个局域 $\sigma_z$ 期望却都为 0，因为各自约化态为 I/2。
同样 $\sigma_x\otimes\sigma_x$ 交换 00 与 11，Bell 态仍不变，相关期望也为 1。
经典混合态 $(|00\rangle\langle00|+|11\rangle\langle11|)/2$ 的 x 相关期望为 0：
对 00 作用得到 11，与 00 内积为零；对 11 作用得到 00，与 11 内积为零，再取平均仍为零。
这给出只测局域态不能看到、联合测量可以区分的具体相干信息。
局域酉变换把振幅矩阵变成 $C'=LCR^T$，
因为输出系数为 $\sum_{ij}L_{ki}C_{ij}R_{lj}$；$R^T$ 也是酉阵，
所以由 SVD 的酉不变性，Schmidt 系数不变，局域换基不能消除纠缠纯态的非零 Schmidt 分量。
"""

REFERENCE = '''def tensor_tools(A, B, rho, dims, keep):
    A = np.asarray(A, dtype=complex)
    B = np.asarray(B, dtype=complex)
    rho = np.asarray(rho, dtype=complex)
    da, db = dims
    if A.ndim != 2 or B.ndim != 2 or da < 1 or db < 1 or rho.shape != (da*db, da*db) or keep not in (0, 1):
        raise ValueError("维数或保留系统不合法")
    ar, ac = A.shape
    br, bc = B.shape
    product = np.zeros((ar*br, ac*bc), dtype=complex)
    for i in range(ar):
        for j in range(ac):
            for k in range(br):
                for l in range(bc):
                    product[i*br+k, j*bc+l] = A[i, j]*B[k, l]
    size = da if keep == 0 else db
    reduced = np.zeros((size, size), dtype=complex)
    for i in range(size):
        for k in range(size):
            if keep == 0:
                for j in range(db):
                    reduced[i, k] += rho[i*db+j, k*db+j]
            else:
                for j in range(da):
                    reduced[i, k] += rho[j*db+i, j*db+k]
    return product, reduced
'''
CASES = '''def make_cases(rng):
    cases = []
    for k in range(32):
        A = rng.normal(size=(1+k%3, 1+(k//3)%3)) + 1j*rng.normal(size=(1+k%3, 1+(k//3)%3))
        B = rng.normal(size=(1+(k//2)%3, 1+(k//4)%3)) + 1j*rng.normal(size=(1+(k//2)%3, 1+(k//4)%3))
        da, db = 1+k%3, 1+(k//3)%3
        M = rng.normal(size=(da*db, da*db)) + 1j*rng.normal(size=(da*db, da*db))
        rho = M @ M.conj().T
        rho /= np.trace(rho)
        cases.append((A, B, rho, (da, db), k%2))
    psi = np.array([1., 0., 0., 1.])/np.sqrt(2)
    for keep in (0, 1):
        cases.append((np.eye(2), np.array([[0., 1.], [1., 0.]]), np.outer(psi, psi), (2, 2), keep))
    return cases
'''
WRITTEN = [
    {"id": "ch10-w1", "level": "基础", "problem_md": r"$\dim V=2,\dim W=3$。写出按第一指标为外层排序的张量积基，并证明维数。", "solution_md": r"""
1. 依据配对构造，基候选为 $e_1\otimes f_1,e_1\otimes f_2,e_1\otimes f_3,e_2\otimes f_1,e_2\otimes f_2,e_2\otimes f_3$。
2. 依据自由系数表定义，每个元素 $\sum c_{ij}g_{ij}$ 都是该组的线性组合，故张成。
3. 依据逐坐标比较，零组合在六个槽位都为零，所以每个系数为零，故无关。
4. 依据基的长度定义维数，这六个向量组成基，$\dim(V\otimes W)=6$；不能把直和维数 2+3=5 套过来。
""", "pitfalls_md": "张量积基枚举全部配对；m+n 计算的是直和的基数。"},
    {"id": "ch10-w2", "level": "基础", "problem_md": r"$A=\begin{pmatrix}1&2\\0&-1\end{pmatrix},B=\begin{pmatrix}0&1\\1&0\end{pmatrix}$。求 A⊗B 及其对 |10⟩ 的作用。", "solution_md": r"""
1. 依据分块定义，$A\otimes B=\begin{pmatrix}B&2B\\0&-B\end{pmatrix}$。
2. 依据逐块代入，矩阵为 $\begin{pmatrix}0&1&0&2\\1&0&2&0\\0&0&0&-1\\0&0&-1&0\end{pmatrix}$。
3. 依据基序 00、01、10、11，|10⟩ 是第三个标准列，像为第三列 $(0,2,0,-1)^T$，即 $2|01\rangle-|11\rangle$。
4. 依据局域规则，$A|1\rangle=2|0\rangle-|1\rangle$、$B|0\rangle=|1\rangle$，两端张量积也得同一结果，核对分块索引。
""", "pitfalls_md": "第一系统指标变化对应整块移动；基序交换会改变扁平矩阵表示。"},
    {"id": "ch10-w3", "level": "进阶", "problem_md": r"证明 $(|00\rangle+|11\rangle)/\sqrt2$ 不可写成单个直积，并比较 $(|00\rangle+|01\rangle)/\sqrt2$。", "solution_md": r"""
1. 依据假设直积展开，四个系数为 ac、ad、bc、bd，必须满足 $ac=bd=1/\sqrt2,ad=bc=0$。
2. 依据 ac≠0，a、c 都非零；由 ad=0 及 a≠0 得 d=0。
3. 依据 d=0 必有 bd=0，与 bd=1/√2 矛盾，故第一个态纠缠。
4. 依据双线性，第二个态等于 $|0\rangle\otimes(|0\rangle+|1\rangle)/\sqrt2$，是直积态，说明叠加不自动意味着纠缠。
""", "pitfalls_md": "系数允许复数也不能解除乘积方程的矛盾；不能仅凭出现加号判纠缠。"},
    {"id": "ch10-w4", "level": "进阶", "problem_md": r"$\psi=(|01\rangle+2i|10\rangle)/\sqrt5$。求两端约化密度矩阵及纯度。", "solution_md": r"""
1. 依据外积，ρ 只有 $\rho_{01,01}=1/5,\rho_{10,10}=4/5,\rho_{01,10}=-2i/5,\rho_{10,01}=2i/5$ 非零。
2. 依据对 B 的同指标求和，$\rho_A=\operatorname{diag}(1/5,4/5)$，非对角元各为零，因为相关的联合条目不属于同一个 B 指标。
3. 依据对 A 的同指标求和，$\rho_B=\operatorname{diag}(4/5,1/5)$，两端非零本征值相同但坐标位置不同。
4. 依据平方求迹，两端纯度均为 $(1+16)/25=17/25<1$；已知联合态纯，所以约化态混合证明纠缠。
""", "pitfalls_md": "A、B 的约化矩阵不一定逐元素相同，即使具有相同非零谱。"},
    {"id": "ch10-w5", "level": "考试风格", "problem_md": r"对 $\psi=(|00\rangle+|01\rangle+|10\rangle-|11\rangle)/2$，用振幅矩阵求 Schmidt 系数、约化态和纠缠判定。", "solution_md": r"""
1. 依据联合基系数，$C=\frac12\begin{pmatrix}1&1\\1&-1\end{pmatrix}$。
2. 依据矩阵乘法，$CC^\dagger=I_2/2$，其本征值为 1/2、1/2，所以两个 Schmidt 系数均为 1/√2。
3. 依据列的正交归一化，可取 $u_1=(1,1)^T/\sqrt2,u_2=(1,-1)^T/\sqrt2,w_1=|0\rangle,w_2=|1\rangle$，于是 $\psi=(u_1\otimes w_1+u_2\otimes w_2)/\sqrt2$。
4. 依据部分迹谱，$\rho_A=\rho_B=I_2/2$，纯度为 1/2，Schmidt 秩为 2，所以该归一化联合纯态纠缠。
""", "pitfalls_md": "四个非零系数本身不是纠缠证明；关键是振幅矩阵的秩或 Schmidt 系数。"},
    {"id": "ch10-w6", "level": "考试风格", "problem_md": r"比较 Bell 纯态与 $\rho_c=(|00\rangle\langle00|+|11\rangle\langle11|)/2$ 的局域态和 $\sigma_x\otimes\sigma_x$ 期望值，解释纯态判据的限制。", "solution_md": r"""
1. 依据部分迹同指标求和，两者的两个局域矩阵均为 I₂/2，所以只测单粒子不能区分。
2. 依据 $\sigma_x\otimes\sigma_x$ 交换 00 与 11，Bell 态是其本征值 +1 的态，期望值为 1。
3. 依据混合态的线性期望值，ρc 中两个分量分别给 $\langle00|11\rangle=0$ 与 $\langle11|00\rangle=0$，平均仍为 0。
4. 依据 ρc 已表示成两个直积纯态的正权重混合，它是可分混合态；其约化态虽混合却不证明纠缠，故纯度判据必须附联合态为纯态的前提。
""", "pitfalls_md": "部分迹丢失联合相干信息；相同约化态不意味着相同联合态或相同纠缠性质。"},
]
QUIZ = [
    {"q_md": "3 维与 4 维空间的张量积维数为？", "options": ["7", "12", "4", "3"], "answer": 1, "explain_md": "3×4 对基符号构成无关张成组，维数为乘积。"},
    {"q_md": "哪一个二比特纯态是直积态？", "options": [r"$(|00\rangle+|11\rangle)/\sqrt2$", r"$(|01\rangle+|10\rangle)/\sqrt2$", r"$(|00\rangle+|01\rangle)/\sqrt2$", r"$(|00\rangle-|11\rangle)/\sqrt2$"], "answer": 2, "explain_md": "第三项可以提出 |0⟩，写成 |0⟩⊗(|0⟩+|1⟩)/√2。"},
    {"q_md": r"$\operatorname{Tr}_B\rho$ 的条目是？", "options": [r"$\sum_j\rho_{(i,j),(k,j)}$", r"$\sum_{j,l}\rho_{(i,j),(k,l)}$", r"$\operatorname{tr}\rho$", r"$\rho^2$"], "answer": 0, "explain_md": "保留 i、k，只对同一个 B 指标 j 求对角配对和。"},
    {"q_md": "已知联合纯态，约化态纯当且仅当？", "options": ["联合态是 Bell 态", "局域维数相等", "Schmidt 秩最大", "联合态是直积态"], "answer": 3, "explain_md": "约化本征值为 sα²，纯度为 1 当且仅当只有一个非零 Schmidt 系数。"},
    {"q_md": "C=UΣV† 对应的 B 端 Schmidt 基取什么？", "options": ["U 的行", "V 列的复共轭", "V 的特征值", "Σ 的行和"], "answer": 1, "explain_md": "Cij=Σα U_iα sα overline(V_jα)，因此 B 端系数为 V 列的共轭。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch10-q{_i}"

CHAPTER = {
    "id": "ch10", "num": 10, "title": "张量积与复合量子系统", "depends_on": ["ch2", "ch8", "ch9"],
    "theory_md": THEORY, "viz": ["tensor_dims"], "written": WRITTEN, "quiz": QUIZ,
    "code_exercises": [{
        "id": "ch10-c1", "title": "手写 Kronecker 分块与两个方向的部分迹",
        "prompt_md": r"实现 `tensor_tools(A, B, rho, dims, keep)`，返回 `(A⊗B, reduced)`。A、B 为行列数各 1～3 的复矩阵；rho 是两端维数 dims=(da,db) 的联合密度矩阵，da、db 各 1～3，与 A、B 的维数独立。keep=0 保留 A 系统、对 B 求迹；keep=1 保留 B 系统、对 A 求迹。按 i*db+j 排联合基；Kronecker 积按每个 Aij 乘整块 B 的定义填入。部分迹必须手写同指标求和。禁止 np.kron（含别名与内部重导出）、outer/multiply.outer、einsum、tensordot、trace/diagonal 等代算路径以及外部量子库；可用数组分配、索引、循环、标量乘加。rho 维数不符或 keep 非 0/1 时抛 ValueError。",
        "function_name": "tensor_tools", "starter": "def tensor_tools(A, B, rho, dims, keep):\n    # TODO：按块填积；分别按被消去的系统索引求部分迹\n    return np.zeros((A.shape[0]*B.shape[0], A.shape[1]*B.shape[1])), np.zeros((dims[keep], dims[keep]))\n",
        "banned": ["np.kron", "np.outer", "np.multiply.outer", "np.einsum", "np.tensordot", "np.trace", "np.diagonal", "np.linalg", "np.matrix", "np.asmatrix", "np.mat", "np.bmat", "np.matrixlib", "np.matlib", "scipy", "qutip"],
        "reference": REFERENCE, "make_cases": CASES,
    }],
}
