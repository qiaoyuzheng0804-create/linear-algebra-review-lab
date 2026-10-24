"""第 3 章：线性方程组与高斯消元。"""
THEORY = r"""
## 动机与直觉

寻找展开系数、判断态的本征方程、识别约束是否矛盾，都可写为 $Ax=b$。
高斯消元逐次消除变量之间的耦合。本章从解集不变性出发，证明何时有解、何时唯一、
以及自由参数为什么恰好等于未知量数减去秩。

## 定义

**定义 3.1（系统与增广矩阵）**：$A\in\mathbb F^{m\times n},x\in\mathbb F^n,b\in\mathbb F^m$。
$Ax=b$ 为线性方程组，$[A\mid b]$ 为增广矩阵。$b=0$ 称齐次。
**定义 3.2（初等行操作）**：交换两行；一行乘非零标量；把另一行的某标量倍加到本行。
**定义 3.3（阶梯形）**：零行在底部；每个非零行的首非零项称主元，下一行主元严格右移；
主元下方为零。简化阶梯形还要求主元为 $1$ 且主元列其他位置为零。
主元所在列的变量称主元变量，其余称自由变量。秩沿用第 1 章的列空间维数定义。

## 核心推导

### 推导 3.1：行操作为何保持整个解集

动机：消元每一步必须既不引入新解，也不丢掉旧解；只证明旧解仍满足是不够的。
交换第 $i,j$ 个方程仅改变它们的排列，同一个 $x$ 满足前后的条件完全相同；再交换一次恢复原系统。
把第 $i$ 个方程 $a_i x=b_i$ 乘以 $c\ne0$，旧解满足 $ca_ix=cb_i$。
反向乘以 $1/c$ 恢复原等式，所以新解也是旧解。

对 $i\ne j$，把第 $i$ 行改为第 $i$ 行加 $c$ 倍第 $j$ 行：
$$ a_ix=b_i,\quad a_jx=b_j
\quad\longrightarrow\quad (a_i+ca_j)x=b_i+cb_j,\quad a_jx=b_j. $$
旧解将两等式按同样系数组合，满足新系统；新解从新第 $i$ 式减去 $c$ 倍未变的第 $j$ 式，恢复 $a_ix=b_i$。
因此操作逆为加上 $-c$ 倍第 $j$ 行。每个操作都是可逆的，有限次复合仍可逐次倒序撤销，故解集始终相同。
也可记每步左乘一个可逆初等矩阵 $E$，于是 $Ax=b\iff EAx=Eb$。

条件与失效：必须同时操作增广列，且缩放系数不能为零；否则乘以 $1/c$ 的恢复步骤失效。
把本行的 $-1$ 倍加给本行会清零整行，所以第三类操作要求来自另一行。

### 推导 3.2：阶梯形主元数等于秩

动机：希望从消元后看得见的主元数量读出原矩阵的独立信息数。
设可逆行操作总矩阵为 $E$。对任意列系数 $c$，
$$ Ac=0\ \Longrightarrow\ EAc=0,\qquad EAc=0\ \Longrightarrow\ Ac=E^{-1}0=0. $$
所以行变换前后所有列之间的相关关系不变，列空间维数保持。
在阶梯形中设非零行数为 $r$，继续把主元归一并消掉其上方项，可得到相同 $r$ 个主元的简化阶梯形 $R$。
其 $r$ 个主元列在非零行上的坐标构成 $I_r$，因此无关；其他列的底部为零，前 $r$ 个坐标能用主元列展开，
所以所有列均由这些列张成。因此 $\operatorname{rank}R=r$，进而 $\operatorname{rank}A=r$。
这也说明不同合法消元路径虽可能产生不同的中间矩阵，主元总数相同。

条件与失效：用到精确的零判定与可逆消元。部分主元是在当前列选绝对值最大的候选项，以避免除以过小的数；
它改善数值表现，但不能使病态问题变得精确。

### 推导 3.3：Rouché–Capelli 判据与唯一性

动机：把“方程是否互相矛盾”变成一个不依赖消元顺序的秩条件。
$Ax=b$ 等价于 $b=\sum_jx_jA_{:j}$，也就是 $b\in\operatorname{col}A$。
添加列 $b$ 后，如果 $b$ 已在列空间内，张成空间不变，秩不变。
反向，若 $b$ 不在列空间内，取 $A$ 的一组列空间基 $u_1,\ldots,u_r$。
若 $cb+\sum_ja_ju_j=0$ 且 $c\ne0$，则解出 $b=-\sum_j(a_j/c)u_j$，与假设矛盾；
因此 $c=0$，由基无关性所有 $a_j=0$，说明加入 $b$ 会增加一个独立方向。
所以
$$ Ax=b\text{ 有解}\iff \operatorname{rank}[A\mid b]=\operatorname{rank}A.\tag{3.1} $$

消元语言给出同一结论：若出现 $[0\ \cdots\ 0\mid d]$ 且 $d\ne0$，对应 $0=d$，无解，且增广列多出一个主元。
若没有这样的行，设秩为 $r$；任选 $n-r$ 个自由变量，从最后一个主元行向上解出所有主元变量，故存在解。
若 $r=n$，没有自由变量，每个主元变量由回代唯一确定。
若 $r<n$，至少有一个自由变量，取它为 $0$ 或 $1$，回代得到两个不同解，所以非唯一。
综上，
$$ \begin{cases}
\operatorname{rank}A<\operatorname{rank}[A\mid b] &: \text{无解},\\
\operatorname{rank}A=\operatorname{rank}[A\mid b]=n &: \text{唯一解},\\
\operatorname{rank}A=\operatorname{rank}[A\mid b]<n &: \text{无穷多解}.
\end{cases}\tag{3.2} $$

条件与失效：这里在实数或复数域上，存在自由参数时有无穷多个取值；有限域上应改成有限个解。
未知数数目是列数 $n$，不能用方程数 $m$ 替代唯一性判据中的 $n$。

### 推导 3.4：显式构造核的基，证明维数为 n − r

动机：自由变量计数还不足以证明维数，必须展示与每个参数对应的一组无关解向量。
对齐次系统把 $A$ 化为简化阶梯形。设主元列编号为 $p_1,\ldots,p_r$，自由列为 $f_1,\ldots,f_{n-r}$。
第 $i$ 个非零行给出
$$ x_{p_i}+\sum_{j=1}^{n-r}R_{i f_j}x_{f_j}=0,\qquad
x_{p_i}=-\sum_jR_{i f_j}x_{f_j}. $$
令自由参数为 $t_j=x_{f_j}$。构造 $z_j$：其自由坐标满足 $(z_j)_{f_\ell}=\delta_{j\ell}$，
主元坐标为 $(z_j)_{p_i}=-R_{i f_j}$。逐行代入，主元项 $-R_{i f_j}$ 与自由项 $R_{i f_j}$ 相消，
所以每个 $z_j\in\ker A$。
任意解的各自由坐标是 $t_j$，其主元坐标也由同一表达式确定，因此
$$ x=\sum_{j=1}^{n-r}t_j z_j. $$
这证明张成性。若 $\sum_jc_jz_j=0$，查看自由坐标 $f_\ell$，左边恰为 $c_\ell$，
故每个 $c_\ell=0$，证明无关性。因此
$$ \dim\ker A=n-r=n-\operatorname{rank}A.\tag{3.3} $$
$r=n$ 时得到空基，核只有零向量，维数为 $0$，与公式一致。

非齐次有解时取特解 $x_p$，若 $x$ 也是解，线性性给 $A(x-x_p)=b-b=0$。
反向任意 $z\in\ker A$ 满足 $A(x_p+z)=b$。所以解集为 $x_p+\ker A$，仿射维数也为 $n-r$。
它仅在 $b=0$ 时是向量子空间：若 $b\ne0$，零向量不在解集中。

条件与失效：核维数公式适用于所有有限矩阵，不要求方阵或满秩；非齐次解集的参数化要求先有特解。
无解时无法选择 $x_p$，不能给空解集套用这个仿射维数公式。

## 适用条件与失效情形

行操作的精确结论与浮点实现要分开理解。练习用 $10^{-10}$ 作零阈值、输入尺度适中；
接近该阈值的残差可能被当作零，改变主元或矛盾行判断，应结合输入误差讨论。
本章不会把“有 $n$ 条方程”当作“有 $n$ 条独立约束”。

## 几何直觉

系数行非零的实线性方程表示超平面；零系数行若右端为零，不添加约束，若右端非零则使解集为空。
解集是各方程对应集合的交。
独立约束每增加一个，会在有解的前提下减少一个自由方向。
齐次系统的所有超平面经过原点；非齐次交集可能是平移子空间，也可能为空。

## 与量子力学的连接

有限维本征方程 $H\psi=E\psi$ 等价于 $(H-EI)\psi=0$。
只有当 $\operatorname{rank}(H-EI)<n$ 时才有非零态，且该能量的本征空间维数为
$n-\operatorname{rank}(H-EI)$，即几何简并度。
厄米算符可正交对角化，此时它等于该本征值的代数重数；一般矩阵不能先假设二者相等。
零解始终存在，但不能代表归一化量子态，所以求本征态要求非零核向量。
"""
REFERENCE = '''def gaussian(A, b):
    A, b = np.asarray(A), np.asarray(b)
    m, n = A.shape
    R = np.zeros((m, n+1), dtype=np.result_type(A, b, float))
    R[:, :n], R[:, n] = A, b
    pivots = []
    row = 0
    for col in range(n):
        if row == m:
            break
        pivot = row + int(np.argmax(np.abs(R[row:, col])))
        if abs(R[pivot, col]) <= 1e-10:
            continue
        R[[row, pivot]] = R[[pivot, row]]
        R[row] /= R[row, col]
        for i in range(m):
            if i != row:
                R[i] -= R[i, col]*R[row]
        pivots.append(col)
        row += 1
    x = np.zeros(n, dtype=R.dtype)
    if np.any(np.abs(R[row:, n]) > 1e-10):
        return (0, x, np.zeros((n, 0), dtype=R.dtype))
    free = [j for j in range(n) if j not in pivots]
    Z = np.zeros((n, len(free)), dtype=R.dtype)
    for i, col in enumerate(pivots):
        x[col] = R[i, n]
        for j, f in enumerate(free):
            Z[col, j] = -R[i, f]
    for j, f in enumerate(free):
        Z[f, j] = 1
    return (1 if len(pivots) == n else 2, x, Z)
'''
CASES = '''def make_cases(rng):
    cases = [
        (np.zeros((2, 3)), np.zeros(2)),
        (np.zeros((2, 3)), np.array([0., 1.])),
        (np.array([[0., 2.], [1., 3.]]), np.array([4., 7.])),
        (np.zeros((0, 2)), np.zeros(0)),
        (np.zeros((2, 0)), np.zeros(2)),
    ]
    for index in range(32):
        m, n = int(rng.integers(1, 6)), int(rng.integers(1, 6))
        A = rng.integers(-5, 6, size=(m, n)).astype(float)
        x = rng.integers(-3, 4, size=n).astype(float)
        if index % 4 == 0:
            A = A.astype(complex) + 1j*rng.integers(-2, 3, size=A.shape)
        if index % 3 == 0:
            A[:, -1] = 0
        b = np.sum(A*x, axis=1)
        if index % 5 == 0:
            A[-1] = 0
            b[-1] = 1
        cases.append((A, b))
    return cases
'''
WRITTEN = [
    {"id": "ch3-w1", "level": "基础", "problem_md": r"用消元求解 $x+y=3,\ 2x-y=0$。", "solution_md": r"""
1. 依据定义 3.1，增广矩阵为 $\begin{pmatrix}1&1&3\\2&-1&0\end{pmatrix}$。
2. 依据推导 3.1 的可逆行加法，$R_2\leftarrow R_2-2R_1$，得到第二行 $(0,-3,-6)$。
3. 依据非零行缩放，除以 $-3$，得到 $y=2$。
4. 依据回代法与等式代入规则，将 $y=2$ 代回第一式，得 $x=1$。
5. 依据式 (3.2)，两个主元等于未知量数，且无矛盾行，故 $(1,2)^T$ 是唯一解。
""", "pitfalls_md": r"行相减必须同时改增广列；只改左边会改变问题。"},
    {"id": "ch3-w2", "level": "基础", "problem_md": r"讨论 $x+2y=1,\ 2x+4y=3$ 是否有解，并给出两个秩。", "solution_md": r"""
1. 依据定义 3.1，写出 $[A\mid b]=\begin{pmatrix}1&2&1\\2&4&3\end{pmatrix}$。
2. 依据可逆行加法，$R_2\leftarrow R_2-2R_1$，得到 $\begin{pmatrix}1&2&1\\0&0&1\end{pmatrix}$。
3. 依据推导 3.2，左侧系数矩阵只有一个主元，$\operatorname{rank}A=1$；增广矩阵有两个主元，秩为 $2$。
4. 依据式 (3.1)，两秩不同，因此无解；第二行也直接表达矛盾 $0=1$。
""", "pitfalls_md": r"不能看到左侧两行相关就判有无穷多解，还必须检查右端是否满足同样的关系。"},
    {"id": "ch3-w3", "level": "进阶", "problem_md": r"求 $A=\begin{pmatrix}1&2&1\\0&1&1\end{pmatrix}$ 的核的一组基，并验证零度公式。", "solution_md": r"""
1. 依据齐次系统定义，方程为 $x_1+2x_2+x_3=0$、$x_2+x_3=0$。
2. 依据阶梯形主元位置，$x_1,x_2$ 为主元变量，$x_3=t$ 为自由变量。
3. 依据主元方程回代规则，从第二式得 $x_2=-t$，代回第一式得 $x_1=t$，故全部解为 $t(1,-1,1)^T$。
4. 依据推导 3.4，该非零向量张成核且单向量组无关，所以是核的基。
5. 依据推导 3.2 的主元数判秩与式 (3.3)，两个主元给 $\operatorname{rank}A=2$，未知量数 $n=3$；核维数 $1=3-2$，符合式 (3.3)。
""", "pitfalls_md": r"核的向量有三个分量，因为核在定义域 $\mathbb R^3$，不能按两条方程写成二维向量。"},
    {"id": "ch3-w4", "level": "进阶", "problem_md": r"求 $x+y+z=2,\ 2x+2y+2z=4$ 的全部解，指出解集是否为向量空间。", "solution_md": r"""
1. 依据可逆行操作，第二式减去第一式的两倍得 $0=0$，系数与增广秩均为 $1$。
2. 依据式 (3.2)，系统有解，且 $n-r=2$ 个自由变量。取 $y=s,z=t$。
3. 依据回代规则与特解的定义，回代第一式给 $x=2-s-t$，于是 $x_p=(2,0,0)^T$ 是特解。
4. 依据推导 3.4 的特解加核向量表示，将参数分离，全部解为 $(2,0,0)^T+s(-1,1,0)^T+t(-1,0,1)^T$。
5. 两方向在 $y,z$ 坐标组成单位列，故无关；解集是二维仿射平面。零向量不满足原式，依据向量空间零元条件，它不是向量子空间。
""", "pitfalls_md": r"非齐次解集通常不称“解空间”；其方向空间才是 $\ker A$。"},
    {"id": "ch3-w5", "level": "考试风格", "problem_md": r"讨论参数系统 $x+y+z=1,\ x+ty+z=2,\ x+y+tz=2$ 的有解性及解。", "solution_md": r"""
1. 依据推导 3.1，第二、三式分别减第一式，得到 $(t-1)y=1$、$(t-1)z=1$，第一式保留。
2. 当 $t=1$，后两式都是 $0=1$，依据式 (3.1) 无解；此时系数秩为 $1$，增广秩为 $2$。
3. 依据实数域非零标量可除规则，当 $t\ne1$，可除以 $t-1$，得 $y=z=1/(t-1)$。
4. 依据回代规则与分式通分运算，代回第一式，$x=1-2/(t-1)=(t-3)/(t-1)$。
5. 此时阶梯形有三个非零主元，依据式 (3.2) 解唯一，即 $((t-3)/(t-1),1/(t-1),1/(t-1))^T$。
""", "pitfalls_md": r"参数除法前先分情形；接近 $t=1$ 时解分量变大，不能把浮点近似直接等同于 $t=1$。"},
    {"id": "ch3-w6", "level": "考试风格", "problem_md": r"$H=\operatorname{diag}(E_0,E_0,E_1)$，其中 $E_0\ne E_1$。求能量 $E_0$ 的本征空间、简并度与归一化条件。", "solution_md": r"""
1. 依据本征方程，将 $H\psi=E_0\psi$ 改写为 $(H-E_0I)\psi=0$。
2. 依据对角矩阵相减定义与推导 3.2 的主元判秩，相减得 $H-E_0I=\operatorname{diag}(0,0,E_1-E_0)$；因 $E_1-E_0\ne0$，该矩阵秩为 $1$。
3. 依据式 (3.3)，核维数为 $3-1=2$；逐分量方程要求 $\psi_3=0$，前两个分量自由。
4. 依据本征空间的核定义及定义 1.3 的维数，因此本征空间为 $\operatorname{span}(e_1,e_2)$，几何简并度为 $2$。
5. 依据标准复内积的归一化定义，$\psi=(a,b,0)^T$ 必须满足 $|a|^2+|b|^2=1$，这排除了零解。
""", "pitfalls_md": r"自由参数本身不是概率；复振幅的模平方才进入归一化条件。"},
]
QUIZ = [
    {"q_md": "哪种行操作可能改变解集？", "options": ["交换两行", "一行加另一行的 3 倍", "一行乘以 0", "一行乘以 2"], "answer": 2, "explain_md": "推导 3.1 要求操作可逆；乘以 0 会丢掉约束，无法用倒数恢复。"},
    {"q_md": r"$Ax=b$ 有解的充要条件是？", "options": [r"$m=n$", r"$\operatorname{rank}A=\operatorname{rank}[A\mid b]$", r"$b=0$", r"$\operatorname{rank}A=m$ 必须成立"], "answer": 1, "explain_md": "由式 (3.1)，增广列必须位于原列空间内，等价于添加它不增加秩。"},
    {"q_md": r"$A$ 为 $3\times5$，秩为 2，核维数为？", "options": ["1", "2", "5", "3"], "answer": 3, "explain_md": r"式 (3.3) 使用列数：$\dim\ker A=5-2=3$。"},
    {"q_md": "非齐次方程有唯一解需要？", "options": ["两秩相等且等于未知量数", "只有两秩相等", "方程数大于未知量数即可", "只有系数矩阵为方阵"], "answer": 0, "explain_md": "式 (3.2) 要同时排除矛盾与自由变量，故要求两秩相等且为列数。"},
    {"q_md": r"本征方程 $(H-EI)\psi=0$ 有非零解意味着？", "options": ["核为零空间", "系数矩阵满列秩", r"$\operatorname{rank}(H-EI)<n$", r"$\psi$ 必须为实向量"], "answer": 2, "explain_md": "依据式 (3.3)，非零核要求核维数正，即秩小于未知量数。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch3-q{_i}"

CHAPTER = {
    "id": "ch3", "num": 3, "title": "线性方程组与高斯消元", "depends_on": ["ch1", "ch2"],
    "theory_md": THEORY, "viz": [],
    "code_exercises": [{
        "id": "ch3-c1", "title": "带部分主元的通用高斯消元",
        "prompt_md": r"实现 `gaussian(A, b)`。$A$ 为 $(m,n)$ 数组，$b$ 为长度 $m$ 的一维数组，可为实或复数。逐列选绝对值最大的行主元，阈值 $10^{-10}$，化为简化阶梯形。返回 `(status, x, Z)`：无解时 status=0，x 为 n 维零向量，Z 形状 (n,0)；唯一解时 status=1；无穷多解时 status=2。有解时 x 是自由变量全取零的特解，Z 的每一列是依次令某自由变量为 1、其余为 0 得到的核基，按自由列编号升序排列。唯一解时 Z 形状 (n,0)。允许空行或空列。不得调用 `linalg`，不得用伪逆或最小二乘替代。",
        "function_name": "gaussian", "starter": "def gaussian(A, b):\n    # TODO：同时行变换 A 与 b，记录主元，检查矛盾行\n    # TODO：自由变量取零得到特解，再逐个构造核基\n    n = A.shape[1]\n    return 0, np.zeros(n), np.zeros((n, 0))\n",
        "banned": ["np.linalg", "np.linalg.solve", "np.linalg.inv", "scipy"],
        "reference": REFERENCE, "make_cases": CASES,
    }], "written": WRITTEN, "quiz": QUIZ,
}
