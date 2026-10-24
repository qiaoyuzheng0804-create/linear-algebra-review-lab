"""第 6 章：内积、Gram–Schmidt、QR 与最佳逼近。"""
THEORY = r"""
## 动机与直觉

基给坐标，内积才给长度、夹角与概率幅。非正交坐标常把独立信息混在一起，
Gram–Schmidt 逐个扣掉已经解释的分量，留下新方向。
本章由内积正定性推导 Cauchy–Schwarz，再构造正交基、QR 和最佳逼近，形成同一条推理链。

## 定义

**定义 6.1（内积）**：约定 $\langle x,y\rangle$ 第一变量共轭线性、第二变量线性，
$\langle y,x\rangle=\overline{\langle x,y\rangle}$，且 $\langle x,x\rangle\ge0$，等号当且仅当 $x=0$。
实空间去掉共轭。标准复内积为 $\langle x,y\rangle=x^\dagger y$，
$\|x\|=\sqrt{\langle x,x\rangle}$。
**定义 6.2（正交与正交归一）**：$\langle x,y\rangle=0$ 称正交；
$\langle q_i,q_j\rangle=\delta_{ij}$ 称正交归一组。$W^\perp$ 是与子空间 $W$ 所有向量正交的集合。
**定义 6.3（薄 QR 与正交投影）**：满列秩 $A\in\mathbb F^{m\times n}$（$m\ge n$）的薄 QR
为 $A=QR$，$Q$ 是 $m\times n$ 且 $Q^\dagger Q=I_n$，$R$ 是 $n\times n$ 上三角。
若 $p\in W$ 且 $x-p\in W^\perp$，称 $p$ 是 $x$ 到 $W$ 的正交投影。

## 核心推导

### 推导 6.1：由非负二次式完整证明 Cauchy–Schwarz

动机：要把内积理解为有限概率幅，必须证明其大小不能超过两个态的长度之积。
目标是
$$ |\langle x,y\rangle|^2\le\|x\|^2\|y\|^2.\tag{6.1} $$
先单列共线情形：若 $y=0$，两边都是零；若 $y\ne0$ 且 $x=cy$，
第一变量共轭线性给 $\langle x,y\rangle=\bar c\|y\|^2$，
而 $\|x\|^2=|c|^2\|y\|^2$，故两边都为 $|c|^2\|y\|^4$，等号成立。
这也包含 $x=0$。

现在令 $x,y$ 不共线，所以二者非零。实内积下，对每个实数 $t$，
$$ f(t)=\|x+ty\|^2=\|y\|^2t^2+2\langle x,y\rangle t+\|x\|^2\ge0. $$
展开来自内积的两变量分配律，交叉项相等。首项系数 $a=\|y\|^2>0$。
二次式配方给 $f(t)=a(t+b/(2a))^2+c-b^2/(4a)$；
在 $t=-b/(2a)$ 处非负，故 $b^2-4ac\le0$。
代入 $b=2\langle x,y\rangle,c=\|x\|^2$，除以 $4$，就得到式 (6.1)。
若判别式为零，$f(-\langle x,y\rangle/\|y\|^2)=0$，
正定性迫使 $x-\langle x,y\rangle y/\|y\|^2=0$，违反不共线。
因此不共线时严格小于，共线时恰等号。

复内积不能把复数系数直接放进实二次式判别式。记 $z=\langle x,y\rangle$。
若 $z=0$，式 (6.1) 左边为零，不共线时右边为正，结论成立。
若 $z\ne0$，设 $y'=\bar z y/|z|$；第二变量线性给
$\langle x,y'\rangle=\bar z z/|z|=|z|$，且 $\|y'\|=\|y\|$。
于是对实数 $t$，
$$ \|x+ty'\|^2=\|y\|^2t^2+2|z|t+\|x\|^2\ge0. $$
同一配方判别式论证给 $|z|^2\le\|x\|^2\|y\|^2$。
等号时 $x-(|z|/\|y\|^2)y'=0$，说明共线；结合前面的直接计算，等号条件双向成立。

条件与失效：正定性既保证二次式非负，也把零范数推出零向量；不定内积不满足这两步。
复数情形必须先调相位，否则实 $t$ 只约束内积实部，无法控制模长。

### 推导 6.2：Gram–Schmidt 的构造、正交性与张成不变

动机：已有方向上的分量会造成冗余；要从 $v_k$ 减去怎样的组合，才能只保留垂直的新信息？
设已经有正交归一的 $q_1,\ldots,q_{k-1}$。要求
$u_k=v_k-\sum_{j<k}c_jq_j$ 与每个 $q_i$ 正交，逐个取内积得
$$ 0=\langle q_i,u_k\rangle=\langle q_i,v_k\rangle-\sum_{j<k}c_j\delta_{ij}
=\langle q_i,v_k\rangle-c_i. $$
因此系数由条件唯一确定：$c_i=\langle q_i,v_k\rangle$，不是任意选择。
构造为
$$ u_k=v_k-\sum_{j<k}\langle q_j,v_k\rangle q_j,\qquad q_k=u_k/\|u_k\|.\tag{6.2} $$

对输入无关组 $v_1,\ldots,v_n$ 作归纳。
起点 $v_1\ne0$，$q_1=v_1/\|v_1\|$ 是单位向量，且二者张成同一直线。
假设前 $k-1$ 步已正交归一且张成 $W_{k-1}=\operatorname{span}(v_1,\ldots,v_{k-1})$。
上面的逐项内积计算给 $u_k$ 与全部旧 $q_i$ 正交。
若 $u_k=0$，式 (6.2) 会把 $v_k$ 表成旧 $q$ 的组合，进而落入 $W_{k-1}$，与输入无关矛盾；
所以可以除以其非零范数，得到单位 $q_k$，且保持与旧方向正交。
又因为 $q_k$ 是 $v_k$ 与旧 $q$ 的组合，新 $q$ 全部落在 $\operatorname{span}(v_1,\ldots,v_k)$；
反向重排构造式为
$v_k=\sum_{j<k}\langle q_j,v_k\rangle q_j+\|u_k\|q_k$，
加上归纳假设，所有旧 $v$ 与新 $v_k$ 都在新 $q$ 的张成内。
两个包含关系给张成空间相等，归纳完成。

经典法用原始 $v_k$ 一次计算全部投影；修正版每扣去一个方向后，用当前残差计算下一投影。
精确算术中，扣掉 $q_i$ 不改变与其他 $q_j$ 的内积，因为 $\langle q_j,q_i\rangle=0$；
所以两者代数等价。练习采用修正版，减轻浮点中旧投影残留造成的影响。

条件与失效：输入相关时出现 $u_k=0$，归一化那一步不能执行；可以跳过以构造子空间基，
但本章练习的薄 QR 契约要求满列秩，遇到零残差应报告相关，不能伪造一个单位列。

### 推导 6.3：由逐列展开得到 QR 分解

动机：把正交化过程保存在矩阵中，既记录新基，也记录原始坐标的恢复关系。
将 $v_j=A_{:j}$ 代入式 (6.2) 并重排，得
$$ v_j=\sum_{i<j}q_i r_{ij}+q_j r_{jj},\qquad
r_{ij}=\langle q_i,v_j\rangle\ (i<j),\quad r_{jj}=\|u_j\|>0. $$
令 $r_{ij}=0$（$i>j$），$Q=(q_1,\ldots,q_n)$，$R=(r_{ij})$。
$QR$ 的第 $j$ 列按矩阵乘法就是上述展开，因此每列相等，$A=QR$。
正交归一性逐条给 $(Q^\dagger Q)_{ij}=\langle q_i,q_j\rangle=\delta_{ij}$。
$R$ 主对角线非零，上三角回代可唯一解出每个坐标，因而 $R$ 可逆。

条件与失效：薄 $Q$ 一般不是方阵，只能断言 $Q^\dagger Q=I_n$；
$QQ^\dagger$ 是 $m$ 维空间中的投影，$m>n$ 时不能把它当 $I_m$。
数值满秩不能仅凭列数判断，仍须看残差是否落入给定阈值。

### 推导 6.4：误差正交当且仅当最佳逼近

动机：正交投影为何给最短误差？必须证明充分性与必要性，而不只画一条垂线。
设 $W$ 有正交归一基 $q_1,\ldots,q_n$，令
$$ p=\sum_jq_j\langle q_j,x\rangle=QQ^\dagger x,\qquad r=x-p.\tag{6.3} $$
逐项计算 $\langle q_i,r\rangle=\langle q_i,x\rangle-\langle q_i,x\rangle=0$，
由线性组合可知 $r$ 与所有 $w\in W$ 正交。
对任意 $y\in W$，$x-y=r+(p-y)$，交叉内积为零，因此
$$ \|x-y\|^2=\|r\|^2+\|p-y\|^2\ge\|r\|^2.\tag{6.4} $$
正定性给等号当且仅当 $p-y=0$，所以 $p$ 是唯一最佳逼近。

反向，设 $p\in W$ 使误差最小，令 $r=x-p$。对任意 $w\in W$ 与实数 $t$，$p+tw\in W$。
最小性要求
$$ 0\le\|r-tw\|^2-\|r\|^2=t^2\|w\|^2-2t\operatorname{Re}\langle w,r\rangle. $$
若实部非零，取同号且足够小的 $t$，一次项绝对值大于二次项，右端会为负，矛盾。
故实部为零。复空间还可用 $iw\in W$，而
$\langle iw,r\rangle=-i\langle w,r\rangle$ 的实部就是原内积的虚部，故虚部也为零。
所以 $r\perp W$；实空间只需第一部分。双向证明完成。

应用到最小二乘：$W=\operatorname{col}A$，最佳残差 $b-A\hat c$ 与每列正交，逐列就是
$A^\dagger(b-A\hat c)=0$，即 $A^\dagger A\hat c=A^\dagger b$。
满列秩时对任意非零 $c$，$c^\dagger A^\dagger Ac=\|Ac\|^2>0$，所以 $A^\dagger A$ 核为零、可逆。
也可由 $A=QR$ 直接解 $R\hat c=Q^\dagger b$，避免形成正规方程的矩阵积。

条件与失效：正交投影与最佳逼近在有限维子空间上存在且唯一；
若 $A$ 不满列秩，投影点仍唯一，但系数 $\hat c$ 可以沿 $\ker A$ 改变，正规方程不可直接求逆。
反向证明在复空间必须检验 $iw$，否则只证明了实部正交。

### 手算算例：从正交化一直算到投影

令 $v_1=(1,1,0)^T,v_2=(1,0,1)^T$，$A=(v_1,v_2)$。
由定义 $\|v_1\|=\sqrt2$，所以 $q_1=(1,1,0)^T/\sqrt2$，$r_{11}=\sqrt2$。
取内积 $r_{12}=q_1^Tv_2=1/\sqrt2$，扣除投影得
$u_2=v_2-r_{12}q_1=(1/2,-1/2,1)^T$。
其范数平方 $1/4+1/4+1=3/2$，故 $r_{22}=\sqrt6/2$，$q_2=(1,-1,2)^T/\sqrt6$。
逐项核对 $q_1^Tq_2=(1-1+0)/\sqrt{12}=0$，且 $q_2^Tq_2=(1+1+4)/6=1$。
因此
$$ Q=\begin{pmatrix}1/\sqrt2&1/\sqrt6\\1/\sqrt2&-1/\sqrt6\\0&2/\sqrt6\end{pmatrix},\quad
R=\begin{pmatrix}\sqrt2&1/\sqrt2\\0&\sqrt6/2\end{pmatrix}. $$
$QR$ 第一列为 $v_1$，第二列为 $(1/2,1/2,0)^T+(1/2,-1/2,1)^T=v_2$。
对 $b=(1,2,3)^T$，$Q^Tb=(3/\sqrt2,5/\sqrt6)^T$。
所以 $p=Q Q^Tb=(3/2,3/2,0)^T+(5/6,-5/6,10/6)^T=(7/3,2/3,5/3)^T$。
残差 $r=(-4/3,4/3,4/3)^T$，$v_1^Tr=0$、$v_2^Tr=0$，最小误差平方是 $3\cdot16/9=16/3$。
若要系数，先回代 $r_{22}\hat c_2=5/\sqrt6$ 得 $\hat c_2=5/3$，
再由 $\sqrt2\hat c_1+\hat c_2/\sqrt2=3/\sqrt2$ 得 $\hat c_1=2/3$。
于是 $A\hat c=(2/3)v_1+(5/3)v_2=p$，整个流程闭合。

## 适用条件与失效情形

内积的共轭约定贯穿投影系数，不能把 $\langle q,v\rangle$ 写成 $\langle v,q\rangle$ 而不取共轭。
近相关列会产生很小的残差，归一化会放大舍入误差；修正版改善表现，但不会消除病态性。
练习阈值为 $10^{-10}$，输入尺度适中且远离相关边界；精确结论中的非零与数值阈值需分清。

## 几何直觉

GS 每步先沿旧子空间滑动，减去平行部分，再把剩余垂直部分改成单位长度。
`gram_schmidt` 用步骤滑杆显示 $v_1,v_2$ 和扣除投影后的残差。
最佳逼近的直角三角形由式 (6.4) 精确表达：在投影点以外走任何一段，都额外增加一个非负的平方长度。

## 与量子力学的连接

归一化态 $|\psi\rangle$ 在正交归一基 $|n\rangle$ 中写成
$|\psi\rangle=\sum_n c_n|n\rangle$，取内积给 $c_n=\langle n|\psi\rangle$。
Born 规则把 $|c_n|^2$ 解释为该基测量的概率；Cauchy–Schwarz 对单位态给 $|c_n|\le1$，
而展开内积给 $\sum_n|c_n|^2=1$。
测量结果属于一组基态张成的子空间 $W$ 时，投影算符为
$P_W=\sum_{n\in W}|n\rangle\langle n|$。
由基正交性逐项相乘得 $P_W^2=P_W$，取伴随得 $P_W^\dagger=P_W$。
概率为 $\langle\psi|P_W|\psi\rangle=\|P_W\psi\|^2$；若非零，选择该结果后的态是
$P_W|\psi\rangle/\|P_W\psi\|$。
例如 $|\psi\rangle=(|0\rangle+i|1\rangle)/\sqrt2$，测得 $|1\rangle$ 的概率为 $|i/\sqrt2|^2=1/2$。
非正交态不能直接套用这些投影求和式，失效点正是 $\langle q_i,q_j\rangle=\delta_{ij}$ 的消项步骤。
"""

REFERENCE = '''def gram_schmidt_qr(A):
    A = np.asarray(A)
    m, n = A.shape
    if m < n:
        raise ValueError("要求满列秩且 m>=n")
    dtype = np.result_type(A, float)
    Q = np.zeros((m, n), dtype=dtype)
    R = np.zeros((n, n), dtype=dtype)
    for j in range(n):
        v = np.array(A[:, j], dtype=dtype, copy=True)
        for i in range(j):
            R[i, j] = np.sum(np.conj(Q[:, i]) * v)
            v -= R[i, j] * Q[:, i]
        length = np.sqrt(np.sum(np.abs(v)**2))
        if length <= 1e-10:
            raise ValueError("列相关或残差过小")
        R[j, j] = length
        Q[:, j] = v / length
    return Q, R
'''
CASES = '''def make_cases(rng):
    cases = [(np.zeros((3, 0)),), (np.array([[1., 1.], [1., 0.], [0., 1.]]),),
             (np.array([[1j], [1.]]),)]
    for index in range(29):
        n = int(rng.integers(1, 5))
        m = n + int(rng.integers(0, 4))
        A = .2*rng.normal(size=(m, n))
        A[:n] += 3*np.eye(n)
        if index % 2:
            A = A + .2j*rng.normal(size=A.shape)
        cases.append((A,))
    return cases
'''
WRITTEN = [
    {"id": "ch6-w1", "level": "基础", "problem_md": r"$x=(1,i)^T,y=(1,1)^T$，计算内积、范数并核对 Cauchy–Schwarz。", "solution_md": r"""
1. 依据第一变量共轭线性约定，$\langle x,y\rangle=1\cdot1+(-i)\cdot1=1-i$。
2. 依据范数定义，$\|x\|^2=1+|i|^2=2$，$\|y\|^2=1+1=2$。
3. 依据复数模平方，$|1-i|^2=2$，而两范数平方之积为 $4$，所以 $2<4$。
4. 依据式 (6.1) 的等号判据，此处严格小于，二者不共线；也可由第一分量要求系数为 1、第二分量却不相符来核对。
""", "pitfalls_md": "把 i 共轭为 -i；第一变量不取共轭会破坏正定性。"},
    {"id": "ch6-w2", "level": "基础", "problem_md": r"把 $b=(2,0)^T$ 投影到 $W=\operatorname{span}(1,1)^T$，求最小距离。", "solution_md": r"""
1. 依据归一化定义，取 $q=(1,1)^T/\sqrt2$。
2. 依据式 (6.3)，$q^Tb=\sqrt2$，投影 $p=q\sqrt2=(1,1)^T$。
3. 依据向量减法与范数定义，残差 $r=(1,-1)^T$，$q^Tr=0$，距离为 $\sqrt2$。
4. 依据式 (6.4)，任意 $y\in W$ 满足 $\|b-y\|^2=2+\|p-y\|^2$，因此最小距离确为 $\sqrt2$，且只在 $y=p$ 取得。
""", "pitfalls_md": "非单位方向 v 的投影需除以 v 的范数平方；不能直接用 (vᵀb)v。"},
    {"id": "ch6-w3", "level": "进阶", "problem_md": r"对 $A=\begin{pmatrix}1&1\\1&0\end{pmatrix}$ 手算 GS 与薄 QR，规定 R 对角线为正。", "solution_md": r"""
1. 依据式 (6.2)，$r_{11}=\sqrt2$，$q_1=(1,1)^T/\sqrt2$。
2. 依据投影系数公式，$r_{12}=1/\sqrt2$，$u_2=(1,0)^T-(1/2,1/2)^T=(1/2,-1/2)^T$。
3. 依据归一化定义，$r_{22}=1/\sqrt2$，$q_2=(1,-1)^T/\sqrt2$。
4. 依据逐列展开，$Q=\frac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix}$，$R=\begin{pmatrix}\sqrt2&1/\sqrt2\\0&1/\sqrt2\end{pmatrix}$，其乘积两列分别为 $(1,1)^T,(1,0)^T$。
""", "pitfalls_md": "余量归一化后的长度记录在 R 对角线上，不能把它丢掉。"},
    {"id": "ch6-w4", "level": "进阶", "problem_md": r"对 $v_1=(1,2)^T,v_2=(2,4)^T$ 执行 GS，说明哪一步失败及子空间维数。", "solution_md": r"""
1. 依据范数定义，$\|v_1\|=\sqrt5$，故 $q_1=(1,2)^T/\sqrt5$。
2. 依据投影公式，$\langle q_1,v_2\rangle=(2+8)/\sqrt5=2\sqrt5$。
3. 依据式 (6.2)，$u_2=v_2-2\sqrt5q_1=(2,4)^T-(2,4)^T=0$。
4. 依据正定性，零残差没有可除的范数，不能构造第二单位列；又因 $v_2=2v_1$ 且 $v_1\ne0$，张成空间维数为 1。
""", "pitfalls_md": "零残差表示相关，不能随意补一个方向并声称仍是原组的薄 QR。"},
    {"id": "ch6-w5", "level": "考试风格", "problem_md": r"最小化 $\|Ac-b\|$，其中 $A=\begin{pmatrix}1&0\\0&1\\1&1\end{pmatrix},b=(1,1,0)^T$，给出最小二乘解并证明最优。", "solution_md": r"""
1. 依据残差正交的必要条件，$A^TAc=A^Tb$；逐项乘法给 $\begin{pmatrix}2&1\\1&2\end{pmatrix}c=(1,1)^T$。
2. 依据两式相减，$c_1-c_2=0$；代回 $2c_1+c_2=1$ 得 $c_1=c_2=1/3$。
3. 依据矩阵乘法，$p=Ac=(1/3,1/3,2/3)^T$，残差 $r=b-p=(2/3,2/3,-2/3)^T$，与 A 两列点积都为零。
4. 依据式 (6.4)，此点是唯一最佳逼近，最小误差平方为 $4/9+4/9+4/9=4/3$；A 两列无关，故系数也唯一。
""", "pitfalls_md": "正规方程右端是 Aᵀb；最佳点唯一与系数唯一是两个结论，后者需要满列秩。"},
    {"id": "ch6-w6", "level": "考试风格", "problem_md": r"$|\psi\rangle=(1,i,\sqrt2)^T/2$，测量投影 $P=\operatorname{diag}(1,1,0)$。求概率与选择该结果后的归一化态。", "solution_md": r"""
1. 依据态归一化定义，原范数平方为 $(1+1+2)/4=1$。
2. 依据投影作用，$P\psi=(1,i,0)^T/2$，范数平方为 $(1+1)/4=1/2$。
3. 依据 Born 规则，该结果概率为 $\langle\psi|P|\psi\rangle=1/2$。
4. 依据条件态归一化，除以 $\sqrt{1/2}$ 后得到 $(1,i,0)^T/\sqrt2$，其范数平方为 1，且再次作用 P 保持不变。
""", "pitfalls_md": "选择整个二维子空间不会自动选出其中某个基态；归一化除以概率的平方根。"},
]
QUIZ = [
    {"q_md": "标准复内积的共轭放在哪里？", "options": ["第二变量", "第一变量", "两个都取", "都不取"], "answer": 1, "explain_md": "本书采用 x†y，第一变量共轭线性，第二变量线性。"},
    {"q_md": "Cauchy–Schwarz 等号何时成立？", "options": ["向量线性相关，包括零向量", "两个向量正交且均非零", "长度相同即可", "所有实向量"], "answer": 0, "explain_md": "共线直接代入得等号；反向由二次式最小值为零及正定性推出共线。"},
    {"q_md": "GS 的投影系数是？", "options": [r"$\langle v_k,q_j\rangle$", r"$\|q_j\|$", r"$\langle q_j,v_k\rangle$", r"$\|v_k\|$"], "answer": 2, "explain_md": "令 <q_i,u_k>=0 解出该系数；复数情形次序不可倒置。"},
    {"q_md": r"薄 QR 的 $Q\in\mathbb C^{m\times n}$，$m>n$，必有？", "options": [r"$QQ^\dagger=I_m$", r"$Q^{-1}=Q^\dagger$", r"$Q=R$", r"$Q^\dagger Q=I_n$"], "answer": 3, "explain_md": "列正交归一给 Q†Q=I；QQ† 是列空间上的投影。"},
    {"q_md": "有限维子空间最佳逼近的判据是？", "options": ["残差长度为 1", "残差与子空间正交", "系数全为正", "子空间必须是一维"], "answer": 1, "explain_md": "充分性由勾股展开，必要性由沿 w 与 iw 的微小变化证明。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch6-q{_i}"

CHAPTER = {
    "id": "ch6", "num": 6, "title": "内积空间、正交化与 Gram–Schmidt", "depends_on": ["ch1", "ch2", "ch3"],
    "theory_md": THEORY, "viz": ["gram_schmidt"], "written": WRITTEN, "quiz": QUIZ,
    "code_exercises": [{
        "id": "ch6-c1", "title": "修正版 Gram–Schmidt 输出 Q、R",
        "prompt_md": r"实现 `gram_schmidt_qr(A)`，返回 `(Q, R)`，采用修正版 GS：逐列取残差，每减去一个方向后再计算下一投影，复内积用 conj。A 为 m×n、m≥n、满列秩的实或复矩阵，允许 n=0。Q 形状 (m,n)，R 形状 (n,n)，R 对角线为正实数，满足 A≈QR、Q†Q≈I。残差模长≤$10^{-10}$ 时抛出 ValueError；用例尺度适中且远离阈值。不得调用 linalg（含 qr、svd、norm、cholesky）或 scipy，禁用 numpy.matrix 系列封装；模长用平方和开根实现。",
        "function_name": "gram_schmidt_qr", "starter": "def gram_schmidt_qr(A):\n    m, n = A.shape\n    # TODO：当前残差逐次投影，保存 R，归一化得到 Q\n    return np.zeros((m, n)), np.zeros((n, n))\n",
        "banned": ["np.matrix", "np.asmatrix", "np.mat", "np.bmat", "np.matrixlib", "np.matlib", "np.matrix.I", "np.matrix.getI", "np.linalg", "np.linalg.qr", "np.linalg.svd", "np.linalg.norm", "np.linalg.cholesky", "scipy"],
        "reference": REFERENCE, "make_cases": CASES,
    }],
}
