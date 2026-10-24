"""第 2 章：矩阵与线性变换。"""
THEORY = r"""
## 动机与直觉

矩阵不是数字表的任意排列，而是线性映射在基下的坐标记录。本章解释为什么只需知道基向量的像，
为什么复合映射产生矩阵乘法，以及同一个算符换基后为何出现相似变换。

## 定义

**定义 2.1（线性映射与矩阵表示）**：$T:V\to W$ 满足 $T(au+bv)=aT(u)+bT(v)$。
$B=(b_1,\ldots,b_n)$、$C=(c_1,\ldots,c_m)$ 分别为定义域与陪域的有序基。
表示矩阵记为 $[T]_{C\leftarrow B}$，约定其作用是输入 $B$ 坐标、输出 $C$ 坐标。

**定义 2.2（矩阵乘法）**：$A\in\mathbb F^{m\times n},D\in\mathbb F^{n\times p}$ 时
$(AD)_{ij}=\sum_{k=1}^n A_{ik}D_{kj}$。维数不匹配时此乘法无定义。

**定义 2.3（换基矩阵与相似）**：$P$ 的第 $j$ 列是新基向量 $b_j$ 在旧基 $E$ 下的坐标。
若 $A_B=P^{-1}A_EP$，称 $A_B,A_E$ 相似。
这里同一个空间同时作定义域与陪域；一般 $V\to W$ 的两端换基要使用两套换基矩阵。

## 核心推导

### 推导 2.1：从基向量的像构造整个映射

动机：线性映射有无限多个输入，但基只有有限多个向量；希望用有限数据完整记录它。
任取 $v\in V$，由基的张成与唯一性，$v=\sum_jx_jb_j$，其中 $x=[v]_B$。
把线性性逐项作用于有限和：
$$ T(v)=T\left(\sum_jx_jb_j\right)=\sum_jx_jT(b_j). $$
将每个像在 $C$ 中唯一展开，设 $T(b_j)=\sum_i a_{ij}c_i$。代入并交换有限求和顺序：
$$ T(v)=\sum_jx_j\sum_i a_{ij}c_i
=\sum_i\left(\sum_j a_{ij}x_j\right)c_i. $$
比较 $C$ 坐标，得到 $[T(v)]_C=Ax$，而第 $j$ 列正是 $[T(b_j)]_C$。
若两个线性映射在每个 $b_j$ 上取值相同，上面的展开给出它们在任意 $v$ 上相同，所以映射唯一。
反过来任意指定 $n$ 个像，按 $T(\sum_jx_jb_j)=\sum_jx_jw_j$ 定义映射：
基坐标唯一保证定义没有歧义；分配律保证它线性。因此这些像既足够又可以自由指定。

条件与失效：$B,C$ 必须为基。若用相关生成组，不同展开可能对指定的像产生不同输出；
例如对 $b_2=2b_1$ 却指定 $T(b_2)\ne2T(b_1)$，就违反了展开一致性那一步。

### 推导 2.2：矩阵乘法来自映射复合

动机：依次实施两个线性操作，应能用一个矩阵记录，且要确定作用顺序。
令 $S:U\to V$、$T:V\to W$，在相容基下的矩阵分别为 $D$ 与 $A$。
对 $U$ 的第 $j$ 个基向量 $d_j$，推导 2.1 给
$$ S(d_j)=\sum_kD_{kj}b_k. $$
再作用 $T$ 并用线性性：
$$ T(S(d_j))=\sum_kD_{kj}T(b_k)
=\sum_kD_{kj}\sum_iA_{ik}c_i
=\sum_i\left(\sum_kA_{ik}D_{kj}\right)c_i. $$
因此复合映射在第 $j$ 个基向量上的像，其第 $i$ 坐标就是定义 2.2 的 $(AD)_{ij}$。
每一列都被确定，故 $[T\circ S]=AD$。右边矩阵先作用。
对任意输入 $x$ 也可逐分量核验：
$$ (A(Dx))_i=\sum_kA_{ik}\sum_jD_{kj}x_j
=\sum_j\left(\sum_kA_{ik}D_{kj}\right)x_j=((AD)x)_i. $$

条件与失效：$S$ 的输出基必须就是 $T$ 的输入基；坐标维数相同但基不同也不能直接相乘，必须先换基。
一般 $AD\ne DA$：交换顺序改变中间输入，甚至会造成维数不匹配。

### 推导 2.3：坐标换基与相似关系

动机：物理算符不应因选坐标而改变，但矩阵条目会改变；需要推导精确关系。
在同一空间上设旧基为 $E$、新基为 $B$，令 $P$ 的列为 $[b_j]_E$。
对 $v=\sum_jx_jb_j$，逐项取旧坐标得
$$ [v]_E=\sum_jx_j[b_j]_E=P[v]_B.\tag{2.1} $$
新基既无关又张成，因此 $P$ 是双射，逆矩阵存在。
设 $A_E=[T]_{E\leftarrow E}$、$A_B=[T]_{B\leftarrow B}$。
先将输入换成旧坐标，再施加映射，最后把输出换回新坐标：
$$ [T(v)]_B=P^{-1}[T(v)]_E
=P^{-1}A_E[v]_E
=P^{-1}A_EP[v]_B. $$
对每个 $[v]_B$ 都成立，特别对每个单位坐标列成立，所以逐列有
$$ A_B=P^{-1}A_EP.\tag{2.2} $$
由左乘 $P$ 得 $PA_B=A_EP$，这也给出计算路线：逐列解 $Pz=A_EP_{:j}$，无需显式求逆。
如果定义域换基矩阵为 $P$、陪域换基矩阵为 $Q$，沿同一输入与输出路线得到
$$ [T]_{\text{新陪域}\leftarrow\text{新定义域}}=Q^{-1}AP.\tag{2.3} $$
仅当两端是同一空间并采用同一换基时，才化为相似关系。

例：$A_E=\operatorname{diag}(2,1)$，$b_1=(1,1)^T,b_2=(1,-1)^T$。
$T(b_1)=(2,1)^T=\frac32b_1+\frac12b_2$，
$T(b_2)=(2,-1)^T=\frac12b_1+\frac32b_2$。
按推导 2.1 排列两列，得 $A_B=\begin{pmatrix}3/2&1/2\\1/2&3/2\end{pmatrix}$。

条件与失效：$P$ 必须可逆。相关的“新基”使从旧坐标恢复新坐标不唯一，$P^{-1}$ 步骤失效。
相似变换不要求 $P$ 正交或酉；把 $P^{-1}$ 换成 $P^\dagger$ 需要额外的正交归一条件。

## 适用条件与失效情形

本章坐标表示适用于有限维线性映射。非线性映射如 $x\mapsto x^2$ 不满足推导 2.1 的分配步骤，
不能用一个固定矩阵表达所有输入；微分得到的雅可比矩阵只描述局部线性近似。
对于近奇异换基矩阵，精确逆虽存在，数值求坐标会放大误差；这与“不成基”是不同层面的失效。

## 几何直觉

二维矩阵的两列分别是横向、纵向单位箭头的像。网格上的每个点都按相同系数组合这两列，
因此直线仍映为直线或点，平行方向保持平行或被压缩为零。
换基可以理解为更换刻度和坐标轴，同一向量及同一映射保持不变。

## 与量子力学的连接

算符 $\hat A$ 在正交归一基 $\{|e_i\rangle\}$ 中的条目是
$A_{ij}=\langle e_i|\hat A|e_j\rangle$，即第 $j$ 个基态的像在第 $i$ 个方向上的振幅。
若新基列构成酉矩阵 $U$，则 $[\psi]_{\rm old}=U[\psi]_{\rm new}$ 且
$A_{\rm new}=U^\dagger A_{\rm old}U$，这是式 (2.2) 在酉换基下的形式。
一般非正交基仍可用相似变换，但振幅内积还需同时记录 Gram 矩阵，不能继续直接用坐标的欧氏内积。
"""
MATMUL_REFERENCE = '''def matmul(A, B):
    A, B = np.asarray(A), np.asarray(B)
    m, n = A.shape
    n2, p = B.shape
    if n != n2:
        raise ValueError("维数不匹配")
    C = np.zeros((m, p), dtype=np.result_type(A, B, float))
    for i in range(m):
        for j in range(p):
            for k in range(n):
                C[i, j] += A[i, k] * B[k, j]
    return C
'''
MATMUL_CASES = '''def make_cases(rng):
    cases = [(np.zeros((0, 2)), np.zeros((2, 3))), (np.zeros((2, 0)), np.zeros((0, 3)))]
    for index in range(32):
        m, n, p = [int(x) for x in rng.integers(1, 6, size=3)]
        A, B = rng.normal(size=(m, n)), rng.normal(size=(n, p))
        if index % 4 == 0:
            A = A + 1j*rng.normal(size=A.shape)
            B = B + 1j*rng.normal(size=B.shape)
        cases.append((A, B))
    return cases
'''
BASIS_REFERENCE = '''def in_basis(A, P):
    A, P = np.asarray(A), np.asarray(P)
    n = P.shape[0]
    R = np.zeros((n, 2*n), dtype=np.result_type(A, P, float))
    R[:, :n] = P
    # 用标量乘加构造 A P，避免把解题工作交给矩阵库。
    for i in range(n):
        for j in range(n):
            for k in range(n):
                R[i, n+j] += A[i, k]*P[k, j]
    for col in range(n):
        pivot = col + int(np.argmax(np.abs(R[col:, col])))
        if abs(R[pivot, col]) <= 1e-10:
            raise ValueError("给定列不构成基")
        R[[col, pivot]] = R[[pivot, col]]
        R[col] /= R[col, col]
        for i in range(n):
            if i != col:
                R[i] -= R[i, col]*R[col]
    return R[:, n:]
'''
BASIS_CASES = '''def make_cases(rng):
    cases = [(np.eye(2), np.eye(2)), (np.diag([2., 1.]), np.array([[1., 1.], [1., -1.]]))]
    for index in range(30):
        n = int(rng.integers(1, 6))
        P = rng.normal(size=(n, n))
        A = rng.normal(size=(n, n))
        if index % 4 == 0:
            P = P + .2j*rng.normal(size=P.shape)
            A = A + 1j*rng.normal(size=A.shape)
        # 严格行对角占优，保证可逆，避免测试用例依赖求逆函数。
        for i in range(n):
            P[i, i] = 1 + np.sum(np.abs(P[i]))
        cases.append((A, P))
    return cases
'''
WRITTEN = [
    {"id": "ch2-w1", "level": "基础", "problem_md": r"已知 $T(e_1)=(1,2)^T,T(e_2)=(-1,3)^T$，求标准基下的矩阵及 $T(2,-1)$。", "solution_md": r"""
1. 依据推导 2.1，矩阵的第 $j$ 列为 $T(e_j)$，故 $A=\begin{pmatrix}1&-1\\2&3\end{pmatrix}$。
2. 依据基展开，$(2,-1)^T=2e_1-e_2$。
3. 依据线性性，$T(2e_1-e_2)=2T(e_1)-T(e_2)$。
4. 依据向量加法与标量乘法的逐分量定义，逐分量相减得到 $(2,4)^T-(-1,3)^T=(3,1)^T$，即所求。
""", "pitfalls_md": r"像应排成列；排成行会表示另一个线性映射。"},
    {"id": "ch2-w2", "level": "基础", "problem_md": r"令 $A=\begin{pmatrix}1&1\\0&1\end{pmatrix},D=\operatorname{diag}(2,3)$。计算 $AD,DA$ 并说明作用顺序。", "solution_md": r"""
1. 依据定义 2.2，$AD$ 第一列为 $A(2,0)^T=(2,0)^T$，第二列为 $A(0,3)^T=(3,3)^T$，故 $AD=\begin{pmatrix}2&3\\0&3\end{pmatrix}$。
2. 依据定义 2.2 的矩阵乘法，计算 $DA$：第一列 $D(1,0)^T=(2,0)^T$，第二列 $D(1,1)^T=(2,3)^T$，故 $DA=\begin{pmatrix}2&2\\0&3\end{pmatrix}$。
3. 依据推导 2.2，$AD$ 表示先缩放 $D$ 再剪切 $A$；$DA$ 表示先剪切再缩放。
4. 依据推导 2.1 的基向量像唯一确定映射，两矩阵的第二列不同，即它们对 $e_2$ 的作用不同，所以映射不相同。
""", "pitfalls_md": r"矩阵写在左边不代表先作用；应从右向左追踪输入。"},
    {"id": "ch2-w3", "level": "进阶", "problem_md": r"$T(x,y)=(2x,y)$，$B=((1,1)^T,(1,-1)^T)$。不求逆，直接求 $[T]_B$。", "solution_md": r"""
1. 依据推导 2.1，先求两个新基向量的像：$T(b_1)=(2,1)^T$，$T(b_2)=(2,-1)^T$。
2. 为求 $T(b_1)$ 的 $B$ 坐标，依据 F5 解 $a+b=2,a-b=1$；两式相加得 $a=3/2$，代回得 $b=1/2$。
3. 依据F5 的基坐标定义与标量消元规则，为求 $T(b_2)$ 的 $B$ 坐标，解 $c+d=2,c-d=-1$，得 $c=1/2,d=3/2$。
4. 依据推导 2.1 的列构造规则，把这两组坐标依次排为列，得 $[T]_B=\begin{pmatrix}3/2&1/2\\1/2&3/2\end{pmatrix}$。
""", "pitfalls_md": r"$T(b_j)$ 的标准坐标还不是所求列，必须再在 $B$ 下展开。"},
    {"id": "ch2-w4", "level": "进阶", "problem_md": r"$T:\mathbb R^2\to\mathbb R^3$ 的旧矩阵为 $A$，两端换基矩阵分别为 $P,Q$。推导新矩阵，并说明为什么不能称一般的相似变换。", "solution_md": r"""
1. 依据式 (2.1)，输入旧坐标 $x=P x'$，其中 $P$ 为 $2\times2$。
2. 依据矩阵表示定义，输出旧坐标 $y=Ax=APx'$。
3. 依据式 (2.1) 的输出换基关系及可逆性，输出换基关系为 $y=Qy'$，其中 $Q$ 为 $3\times3$；左乘 $Q^{-1}$ 得 $y'=Q^{-1}APx'$。
4. 依据推导 2.1 的矩阵表示唯一性，对所有 $x'$ 成立，因此新矩阵是 $Q^{-1}AP$。
5. 依据定义 2.3，相似要求同一空间同一换基矩阵的左右逆配对；此处 $P,Q$ 尺寸不同，新旧矩阵还是 $3\times2$，故不是相似变换。
""", "pitfalls_md": r"左侧换的是输出坐标，右侧换的是输入坐标，不能无条件写成 $P^{-1}AP$。"},
    {"id": "ch2-w5", "level": "考试风格", "problem_md": r"在 $\mathbb P_2$（次数不超过 2 的多项式空间）上定义求导算符 $D$。求基 $(1,x,x^2)$ 下矩阵，并计算 $D^2$ 与其核。", "solution_md": r"""
1. 依据求导定义，$D(1)=0,D(x)=1,D(x^2)=2x$。
2. 依据推导 2.1，把三个像的基坐标排成列：$[D]=\begin{pmatrix}0&1&0\\0&0&2\\0&0&0\end{pmatrix}$。
3. 依据复合映射推导 2.2，$[D^2]=[D][D]=\begin{pmatrix}0&0&2\\0&0&0\\0&0&0\end{pmatrix}$。
4. 依据核的定义及矩阵作用公式，对 $p=a+bx+cx^2$，矩阵作用给 $D^2p=2c$，由核的定义要求 $2c=0$，故 $c=0$。
5. 依据定义 1.2 的张成与无关性判据，因此 $\ker D^2=\operatorname{span}(1,x)$，两个基多项式无关，核维数为 $2$。
""", "pitfalls_md": r"$D^2$ 表示求导算符复合两次，不是把每个矩阵元素平方。"},
    {"id": "ch2-w6", "level": "考试风格", "problem_md": r"$A=\begin{pmatrix}0&1\\1&0\end{pmatrix}$，新正交归一基为 $b_1=(1,1)^T/\sqrt2,b_2=(1,-1)^T/\sqrt2$。求新矩阵，并用该表示计算 $A^k$。", "solution_md": r"""
1. 直接代入得 $Ab_1=b_1$、$Ab_2=-b_2$，依据推导 2.1，新矩阵为 $\Lambda=\operatorname{diag}(1,-1)$。
2. 依据正交归一条件与逆矩阵定义，换基矩阵 $U=(b_1,b_2)$ 满足 $U^TU=I$，因为两列长度为 $1$ 且内积为 $0$；因此 $U^{-1}=U^T$。
3. 依据式 (2.2)，$\Lambda=U^TAU$，等价于 $A=U\Lambda U^T$。
4. 依据矩阵乘法结合律与逆矩阵消去关系，对正整数 $k$，连续相乘时相邻的 $U^TU$ 化为 $I$，得到 $A^k=U\Lambda^kU^T$；$k=0$ 时也给 $I$。
5. 依据对角矩阵幂的逐对角元计算规则，因 $\Lambda^k=\operatorname{diag}(1,(-1)^k)$，故偶数 $k$ 时 $A^k=I$，奇数 $k$ 时 $A^k=A$。
""", "pitfalls_md": r"这里能用转置代逆是因为新基正交归一；一般换基不具备该性质。"},
]
QUIZ = [
    {"q_md": "线性变换矩阵的第 2 列记录什么？", "options": ["第 2 个输出坐标轴", "输入第 2 个基向量之像的输出坐标", "任意向量的第 2 坐标", "矩阵的第 2 行"], "answer": 1, "explain_md": "由推导 2.1，令输入坐标为第二个单位列，输出恰为矩阵第二列。"},
    {"q_md": r"$ADx$ 的作用顺序为？", "options": ["先 A 后 D", "顺序无关", "先 D 后 A", "先转置 A"], "answer": 2, "explain_md": r"由推导 2.2，$AD$ 表示复合 $T\circ S$，输入先经过右侧的 $D$。"},
    {"q_md": r"若 $[v]_E=P[v]_B$，同一算符满足？", "options": [r"$A_B=P^{-1}A_EP$", r"$A_B=PA_EP^{-1}$", r"$A_B=PA_EP$", r"$A_B=A_E$ 恒成立"], "answer": 0, "explain_md": "按式 (2.2)，先转旧坐标、作用旧矩阵、再转回新坐标。"},
    {"q_md": r"一般可逆 $P$ 何时可把 $P^{-1}$ 写成 $P^\dagger$？", "options": ["P 是方阵即可", "P 实数即可", "P 对角即可", "P 为酉矩阵"], "answer": 3, "explain_md": r"酉条件 $P^\dagger P=I$ 给出 $P^\dagger=P^{-1}$，一般可逆性不保证。"},
    {"q_md": "把相关向量组当基的主要问题是？", "options": ["坐标不一定唯一", "线性映射一定非线性", "矩阵一定为单位阵", "所有向量都消失"], "answer": 0, "explain_md": "相关组存在非平凡零组合，可给同一向量添加这组系数而不改变向量，展开失去唯一性。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch2-q{_i}"

CHAPTER = {
    "id": "ch2", "num": 2, "title": "矩阵与线性变换", "depends_on": ["ch1"],
    "theory_md": THEORY, "viz": ["grid_transform"],
    "code_exercises": [
        {"id": "ch2-c1", "title": "三重循环矩阵乘法", "prompt_md": r"实现 `matmul(A, B)`，输入为维数相容的二维实或复数组，允许空维度。按定义用标量乘加与循环计算，禁止 `@`、`np.matmul`、`np.dot`、`np.einsum`、`np.tensordot`、`np.inner`、`np.vecdot`、`np.vecmat`、`np.matvec`、`np.convolve`、`np.correlate` 和 `linalg`。返回形状为 $(m,p)$ 的数组。", "function_name": "matmul", "starter": "def matmul(A, B):\n    m, n = A.shape\n    p = B.shape[1]\n    C = np.zeros((m, p), dtype=np.result_type(A, B, float))\n    # TODO：逐元素计算有限求和\n    return C\n", "banned": ["np.linalg", "np.matmul", "np.dot", "np.einsum", "np.tensordot", "np.inner", "np.vecdot", "np.vecmat", "np.matvec", "np.convolve", "np.correlate", "scipy", "@"], "reference": MATMUL_REFERENCE, "make_cases": MATMUL_CASES},
        {"id": "ch2-c2", "title": "在给定基下表示变换", "prompt_md": r"实现 `in_basis(A, P)`：$A$ 为旧基下算符矩阵，$P$ 的列为新基的旧坐标。输入是同阶方阵，保证 $P$ 可逆且尺度适中，可含复数。返回 $P^{-1}AP$，但不许调用 `linalg`；通过逐列构造 $AP$ 并行化简增广矩阵 $[P\mid AP]$ 完成，使用部分主元与 $10^{-10}$ 零阈值。", "function_name": "in_basis", "starter": "def in_basis(A, P):\n    # TODO：构造右端 AP，将 [P | AP] 行化简为 [I | 新矩阵]\n    return np.zeros_like(A)\n", "banned": ["np.linalg", "scipy"], "reference": BASIS_REFERENCE, "make_cases": BASIS_CASES},
    ], "written": WRITTEN, "quiz": QUIZ,
}
