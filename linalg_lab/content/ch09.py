"""第 9 章：由正半定谱分解构造 SVD、低秩最优逼近和最小二乘。"""
THEORY = r"""
## 动机与直觉

特征分解试图给输入与输出使用同一套轴，矩形矩阵连两端维数都不同，不能要求这种结构。
SVD 允许输入先换基、按非负数伸缩、输出再换基，因此能覆盖矩形与亏秩映射。
这些伸缩量同时回答三个问题：保留哪些方向能最好压缩矩阵、怎样拟合无精确解的方程、
以及二粒子纯态有多少不能消除的相关分量。

## 定义

**定义 9.1（奇异值分解）**：$A\in\mathbb C^{m\times n}$ 的完整 SVD 为
$$ A=U\Sigma V^\dagger,\qquad U^\dagger U=I_m,\quad V^\dagger V=I_n.\tag{9.1} $$
Σ 为 m×n 矩形对角阵，对角元 $\sigma_1\ge\cdots\ge\sigma_p\ge0$，$p=\min(m,n)$。
秩 r 是其中正奇异值个数；为写极小极大公式，把 $\sigma_{p+1},\ldots,\sigma_n$ 视为零。
**定义 9.2（两种范数）**：$\|A\|_F^2=\sum_{i,j}|a_{ij}|^2$，
$\|A\|_2=\max_{\|x\|=1}\|Ax\|$。前者累计所有条目，后者记录最大伸缩。
**定义 9.3（伪逆与最小二乘）**：$A^+=V\Sigma^+U^\dagger$，Σ⁺ 为 n×m，
将 Σ 的非零对角元取倒数并转置位置，零仍置零。
最小二乘求使 $\|Ax-b\|^2$ 最小的 x；不要求 $Ax=b$ 有精确解。

## 核心推导

### 推导 9.1：从 A†A 的谱分解逐步构造 SVD

动机：输入空间中先找出彼此正交的伸缩方向，再把它们的像归一化为输出轴。
令 $G=A^\dagger A$。由乘积伴随逆序，$G^\dagger=A^\dagger A=G$，故第 8 章谱定理适用。
若 $Gv=\lambda v,v\ne0$，左乘 $v^\dagger$ 得
$$ \lambda\|v\|^2=v^\dagger A^\dagger Av=(Av)^\dagger(Av)=\|Av\|^2\ge0. $$
范数平方为正，所以 $\lambda\ge0$；这一步是正半定性的来源。
由谱定理选 $G$ 的正交归一本征基 $v_1,\ldots,v_n$，按特征值降序排列。
令 $\sigma_i=\sqrt{\lambda_i}$。当 $\lambda_i>0$ 时定义 $u_i=Av_i/\sigma_i$。
对正奇异值指标 i、j，逐项计算
$$ u_i^\dagger u_j=\frac{v_i^\dagger A^\dagger Av_j}{\sigma_i\sigma_j}
=\frac{\lambda_j v_i^\dagger v_j}{\sigma_i\sigma_j}
=\frac{\sigma_j^2\delta_{ij}}{\sigma_i\sigma_j}=\delta_{ij}. $$
i=j 时商为 1；i≠j 时分子为零。所以这些 u 正交归一，个数 r 不超过 m。

若 $\lambda_i=0$，上面的非负等式给 $\|Av_i\|^2=0$，正定性使 $Av_i=0$。
因此不能在零奇异值处执行除法；这些方向属于核。
任意 $x=\sum_i(v_i^\dagger x)v_i$，代入 A 得
$Ax=\sum_{i=1}^r\sigma_i u_i(v_i^\dagger x)$。
所以 $u_1,\ldots,u_r$ 张成像空间；反向每个 $u_i=Av_i/\sigma_i$ 也在像空间，
故像空间恰为这些 u 的张成，$\operatorname{rank}A=r$。
在其正交补中取基并用 Gram–Schmidt 扩张成 $u_1,\ldots,u_m$ 的完整正交归一基。
设 U、V 分别以 u、v 为列。Σ 的前 r 个对角元为正奇异值，其余为零。
逐列比较，i≤r 时 $Av_i=\sigma_i u_i$；i>r 时两侧均为零，
于是 $AV=U\Sigma$，右乘 $V^\dagger=V^{-1}$ 得式 (9.1)。
实 A 的 G 为实对称，谱基和正交补都可选实，故得到实正交 SVD。

条件与失效：本构造对任意有限矩阵成立，包括 m≠n、r<p 和零矩阵。
零值处改用正交补，正是避免除以零的步骤。正奇异值重数大于一时，谱基不唯一但仍存在。
特征对角化则一开始要求 $Av=\lambda v$ 两边属于同一空间，必须是方阵；
之后还要求本征向量张满以使换基 P 可逆。SVD 求解的是始终为方阵且厄米的 G，
再分别建立两端基，绕过的正是这两处障碍，因而缺陷 Jordan 矩阵同样有 SVD。

### 推导 9.2：范数、极小极大刻画与截断 SVD 的最优性

动机：压缩时丢掉弱方向听起来合理，但需要证明任何其他秩至多 k 的矩阵都不能做得更好。
先从 SVD 推出范数。令 $c=V^\dagger x$，酉性给 $\|c\|=\|x\|$，于是
$\|Ax\|^2=\sum_i\sigma_i^2|c_i|^2\le\sigma_1^2\|x\|^2$；取 x=v₁ 达到等号，故 $\|A\|_2=\sigma_1$。
Frobenius 范数满足 $\|M\|_F^2=\operatorname{tr}(M^\dagger M)$。
左乘酉阵消去 $U^\dagger U$；右乘酉阵用第 7 章迹的循环公式消去 $VV^\dagger$，
所以 $\|A\|_F^2=\|\Sigma\|_F^2=\sum_{i=1}^n\sigma_i^2$。
这里补到 n 的项若超过 p 就为零。

动机：要约束所有低秩候选，需把第 j 个奇异值表达为子空间上的伸缩极值，才能与候选矩阵的核相交。
下一步完整证明要用的极小极大刻画：
$$ \sigma_j(A)=\min_{\dim S=n-j+1}\ \max_{x\in S,\ \|x\|=1}\|Ax\|,\qquad 1\le j\le n.\tag{9.2} $$
取 $S_0=\operatorname{span}(v_j,\ldots,v_n)$，其上
$\|Ax\|^2=\sum_{i\ge j}\sigma_i^2|c_i|^2\le\sigma_j^2$，取 x=vⱼ 达到等号。
因此右边不超过 $\sigma_j$。反向，对任意维数 n−j+1 的 S，
其与 $L=\operatorname{span}(v_1,\ldots,v_j)$ 的交空间至少一维：
$\dim(S\cap L)=\dim S+\dim L-\dim(S+L)\ge n+1-n=1$。
取交中单位 x，则只有前 j 个谱坐标，
$\|Ax\|^2=\sum_{i\le j}\sigma_i^2|c_i|^2\ge\sigma_j^2\sum_{i\le j}|c_i|^2=\sigma_j^2$。
故每个 S 的最大值至少 $\sigma_j$，合并两方向得到式 (9.2)，且最小值由 $S_0$ 实现。

条件与失效：S 是输入空间中指定维数的线性子空间，单位向量必须在该 S 内取；
超过输出维数的奇异值补零，仍可用核方向达到上界。若把奇异值换成一般矩阵的特征值，
按正交坐标写出的伸缩平方和就没有依据，不能沿用这条极值刻画。

令 $0\le k<p$，$A_k=\sum_{i=1}^k\sigma_i u_iv_i^\dagger$。
由同一酉基下相减，$A-A_k=U(\Sigma-\Sigma_k)V^\dagger$，故
$$ \|A-A_k\|_F=\sqrt{\sum_{i>k}\sigma_i^2},\qquad \|A-A_k\|_2=\sigma_{k+1}.\tag{9.3} $$
仅算出截断误差还不够。取任意 $\operatorname{rank}B\le k$，令 E=A−B。
对每个 $1\le j\le n-k$，由 E 的 SVD 取
$S=\operatorname{span}(v_j(E),\ldots,v_n(E))$，其维数为 n−j+1，
在 S 的单位球上 $\|Ex\|\le\sigma_j(E)$。
由秩零度，$\dim\ker B\ge n-k$，所以
$\dim(S\cap\ker B)\ge n-j-k+1$。
从交空间抽取恰为该维数的子空间 T。对其中每个单位 x，Bx=0，因此
$\|Ax\|=\|Ex\|\le\sigma_j(E)$。将 T 代入 A 的极小极大刻画，得到
$$ \sigma_{j+k}(A)\le\max_{x\in T,\|x\|=1}\|Ax\|\le\sigma_j(E). $$
两端非负，可以逐项平方再求和：
$$ \|A-B\|_F^2=\sum_{j=1}^n\sigma_j(E)^2
\ge\sum_{j=1}^{n-k}\sigma_j(E)^2
\ge\sum_{j=1}^{n-k}\sigma_{j+k}(A)^2
=\sum_{i>k}\sigma_i(A)^2. $$
而 j=1 的不等式直接给 $\|A-B\|_2\ge\sigma_{k+1}(A)$。
截断 Aₖ 同时达到两个下界，这就证明 Eckart–Young 最优逼近定理的这两种范数版本。

条件与失效：需要奇异值按降序排列、约束为秩至多 k，并指定 Frobenius 或谱范数。
若 k≥r，截断已经等于 A，误差为零。若截断边界奇异值相等，最优近似可能不唯一；
证明只给全局下界和一个达到者，没有给唯一性。
逐条目绝对值和等其他误差并不享有上面的酉不变范数公式，平方和那一步不能直接复制。

### 推导 9.3：由误差正交性推出正规方程

动机：有噪声时 b 可能不在列空间，先把“最接近”写成不会被任何列方向改善的条件。
设 $r=b-Ax$。对任意 h 与复标量 t，
$$ \|b-A(x+th)\|^2=\|r-tAh\|^2
=\|r\|^2-2\operatorname{Re}(t\,r^\dagger Ah)+|t|^2\|Ah\|^2. $$
若某 h 使 $z=r^\dagger Ah\ne0$，选 $t=\epsilon\bar z/|z|$，
交叉项就是 $-2\epsilon|z|$，二次项为 $\epsilon^2\|Ah\|^2$。
取足够小的正 $\epsilon$（$Ah\ne0$），总误差严格降低，与最小矛盾。
所以最小点必须有 $r^\dagger Ah=0$ 对所有 h 成立，即 r 与列空间正交。
取共轭得到 $h^\dagger A^\dagger r=0$，再逐个取 h=eⱼ 给 $A^\dagger r=0$。
反向若该式成立，交叉项对任意 h 为零，
$\|r-Ah\|^2=\|r\|^2+\|Ah\|^2\ge\|r\|^2$，所以确为全局最小点。
代入 r 的定义并移项得到充要条件
$$ A^\dagger A x=A^\dagger b.\tag{9.4} $$
这也是 b 到列空间正交投影的系数方程，没有假定 A 满秩。

条件与失效：复数情形必须允许扰动方向的相位变化，仅检验实 t 不能消去虚部交叉项。
若列相关，投影 Ax 唯一但表示它的系数 x 可以不唯一；
因此正规方程有解不代表 $A^\dagger A$ 可逆，不能在此直接写逆矩阵。

### 推导 9.4：伪逆、所有最小二乘解与最小范数

动机：把可确定的系数与核中任意分量分开，给亏秩系统一个唯一的最小范数选择。
将 x 写成 $x=Vc$，b 写成 $b=Ud$，酉性给
$$ \|Ax-b\|^2=\|\Sigma c-d\|^2
=\sum_{i=1}^r|\sigma_i c_i-d_i|^2+\sum_{i=r+1}^m|d_i|^2. $$
第一部分每项都非负且可独立归零，故所有最小点满足 $c_i=d_i/\sigma_i$（i≤r）。
其余输入坐标 i>r 不进入误差，可任取。令它们全为零，得到
$$ x_+=A^+b=\sum_{i=1}^r\frac{u_i^\dagger b}{\sigma_i}v_i,
\qquad x=x_++z,\quad z\in\ker A.\tag{9.5} $$
因为 x₊ 只含前 r 个 v，z 只含后面的 v，二者正交，
$\|x_++z\|^2=\|x_+\|^2+\|z\|^2$，所以 x₊ 是唯一的最小范数最小二乘解。
而 $AA^+b=\sum_{i=1}^r u_i u_i^\dagger b$ 正是到列空间的投影，
残差落在正交补，与推导 9.3 完全一致。

若 A 满列秩，r=n，所有 $\sigma_i>0$，
$A^\dagger A=V\operatorname{diag}(\sigma_i^2)V^\dagger$，可逆且逆为
$V\operatorname{diag}(\sigma_i^{-2})V^\dagger$。
右乘 $A^\dagger=V\Sigma^\dagger U^\dagger$ 时，对角因子逐项给
$\sigma_i^{-2}\sigma_i=\sigma_i^{-1}$，所以
$$ x_+=(A^\dagger A)^{-1}A^\dagger b.\tag{9.6} $$
这是列满秩时的唯一系数解；亏秩时式 (9.5) 仍成立，式 (9.6) 的求逆步骤失效。

### 推导 9.5：伪逆的四恒等式为何唯一确定同一个矩阵

动机：最小范数给出了 $A^+b$，还要核对该矩阵与不依赖 b 的 Moore–Penrose 定义相同。
由 $A=U\Sigma V^\dagger$、$A^+=V\Sigma^+U^\dagger$，直接消去中间酉因子得
$AA^+=U(\Sigma\Sigma^+)U^\dagger$ 与 $A^+A=V(\Sigma^+\Sigma)V^\dagger$。
括号内都是前 r 个对角元为 1、其余为零的矩形两端方阵，因而自伴且幂等。
再逐个对角元算 $\sigma_i(1/\sigma_i)\sigma_i=\sigma_i$ 和
$(1/\sigma_i)\sigma_i(1/\sigma_i)=1/\sigma_i$；零值位置三项乘积也为零，所以
$$ AA^+A=A,\quad A^+AA^+=A^+,\quad
(AA^+)^\dagger=AA^+,\quad(A^+A)^\dagger=A^+A. $$
这同时证明两个乘积确是像空间和核正交补上的正交投影。

为证明唯一性，设 X 为任意满足同四式的 n×m 矩阵，令 $Y=V^\dagger XU$，
把输入与输出都按正奇异方向 r 和其余方向分块。写
$\Sigma=\begin{pmatrix}D&0\\0&0\end{pmatrix}$，
$Y=\begin{pmatrix}Y_{11}&Y_{12}\\Y_{21}&Y_{22}\end{pmatrix}$，
其中 $D=\operatorname{diag}(\sigma_1,\ldots,\sigma_r)$ 可逆。
$AXA=A$ 化为 $DY_{11}D=D$，左右乘 $D^{-1}$ 给 $Y_{11}=D^{-1}$。
$AX$ 自伴化为 $\Sigma Y$ 自伴，其右上块 $DY_{12}$ 必为左下零块的伴随，
故 $Y_{12}=0$。$XA$ 自伴同理使 $Y_{21}D=0$，故 $Y_{21}=0$。
最后 $XAX=X$ 化为 $Y\Sigma Y=Y$，左边右下块为零，给 $Y_{22}=0$。
于是 $Y=\Sigma^+$，换回去 $X=V\Sigma^+U^\dagger=A^+$，唯一性证明结束。
r=0 时 A=0，直接由 $XAX=X$ 得 X=0，也与定义一致；满行或满列时空块不添加条件。

条件与失效：四式须同时满足；只要求 $AXA=A$ 不能固定核与残差方向上的作用。
例如 $A=\operatorname{diag}(1,0)$ 时，任意 $X=\operatorname{diag}(1,c)$ 都满足第一式，
第二式却要求 $c=0$，因此“广义逆”这一较弱条件不等于这里唯一的伪逆。

**完整手算算例。** 取 $A=\begin{pmatrix}1&1\\1&-1\\0&0\end{pmatrix}$，$b=(2,0,3)^T$。
$A^TA=\begin{pmatrix}2&0\\0&2\end{pmatrix}$，故 $V=I_2$，$\sigma_1=\sigma_2=\sqrt2$。
$u_1=(1,1,0)^T/\sqrt2$、$u_2=(1,-1,0)^T/\sqrt2$，再补 $u_3=(0,0,1)^T$。
因此 Σ 前两行是 $\sqrt2 I_2$，第三行为零，逐列 $U\Sigma=A$。
$U^Tb=(\sqrt2,\sqrt2,3)^T$，除以前两个奇异值给 $x_+=(1,1)^T$。
$Ax_+=(2,0,0)^T$，残差 $b-Ax_+=(0,0,3)^T$，范数为 3，
$A^T(b-Ax_+)=(0,0)^T$，验证正规方程。
直接算 $A^Tb=(2,2)^T$，乘 $(A^TA)^{-1}=I_2/2$ 也得 $(1,1)^T$。
保留第一奇异方向的秩一近似为
$A_1=\begin{pmatrix}1&0\\1&0\\0&0\end{pmatrix}$，差矩阵第二列为 $(1,-1,0)^T$，
Frobenius 误差 $\sqrt2$，与丢弃的奇异值一致；相同奇异值使其他秩一最优解也可能存在。

条件与失效：接近零的奇异值取倒数会放大噪声；数值截断相当于改变有效秩，必须说明阈值。
显式形成 G 会把非零奇异值的最大最小比平方，因此本章幂迭代练习限定小型、适中尺度与谱间隔，
用来检验构造链，不将其误当作大型病态问题的稳定 SVD 算法。

## 适用条件与失效情形

SVD 的存在性无须满秩；除以奇异值只出现在正值分支。
$A_k$ 是秩至多 k 的最优逼近，$A^+b$ 是最小范数最小二乘解，这两个“最小”约束对象不同。
代码练习内部构造正奇异三元组，对外返回奇异值、重构和伪逆；不比较列向量的符号，
也不把简并空间里不同正交基误判为错误。零方向由理论中的正交补补齐，精简实现只需秩分解即可重构。

## 几何直觉

交互图 `svd_ellipse` 已按三步展示：先 $V^T$ 旋转或反射，再 Σ 沿坐标轴伸缩，最后 U 旋转或反射。
圆本身在第一步仍是圆，因此图还标记同一个输入向量，显示方向确实发生改变。
椭圆半轴长度就是奇异值；奇异值为零时压成线段或点。矩形情形的输出椭球位于像空间，
额外输出正交方向不由任何输入激发，正是最小二乘残差可能留下的地方。

## 与量子力学的连接

双系统纯态写成 $|\psi\rangle=\sum_{i,j}C_{ij}|i\rangle_A|j\rangle_B$，
对振幅矩阵 C 做 $C=U\Sigma V^\dagger$，逐项展开
$C_{ij}=\sum_\alpha\sigma_\alpha U_{i\alpha}\overline{V_{j\alpha}}$，于是
$$ |\psi\rangle=\sum_\alpha s_\alpha|u_\alpha\rangle_A|w_\alpha\rangle_B,
\quad |w_\alpha\rangle=\sum_j\overline{V_{j\alpha}}|j\rangle_B,\quad s_\alpha=\sigma_\alpha.\tag{9.7} $$
B 端必须用 V 列的共轭，这是 $V^\dagger$ 在系数展开中的来源。
两端各自正交归一，态的范数平方就是 $\sum_\alpha s_\alpha^2$；归一化时此和为 1。
这称 Schmidt 分解，第 10 章将从张量积与部分迹证明其纠缠意义。
例如 $|\psi\rangle=(\sqrt3|00\rangle+|11\rangle)/2$，
$C=\operatorname{diag}(\sqrt3/2,1/2)$，其奇异值就是 $s_1=\sqrt3/2,s_2=1/2$。
约化密度矩阵为 $\rho_A=CC^\dagger=\operatorname{diag}(3/4,1/4)$，
特征值是 Schmidt 系数平方，纯度 $\operatorname{tr}(\rho_A^2)=9/16+1/16=5/8<1$。
一般密度矩阵 ρ 厄米且半正定，谱分解中 $\lambda_i\ge0$，
$\rho^\dagger\rho=\rho^2$ 的本征值为 $\lambda_i^2$，故其奇异值就是 $\sqrt{\lambda_i^2}=\lambda_i$。
对一般非厄米振幅矩阵 C 则不能把奇异值和特征值直接等同。
"""

REFERENCE = '''def svd_from_gram(A):
    A = np.asarray(A, dtype=float)
    m, n = A.shape
    if m not in (2, 3) or n not in (2, 3):
        raise ValueError("只处理二阶或三阶两端维数")
    scale = np.sum(A*A)
    singular = np.zeros(min(m, n))
    reconstructed = np.zeros((m, n))
    inverse = np.zeros((n, m))
    if scale == 0:
        return singular, reconstructed, inverse
    G = A.T @ A / scale
    B = G.copy()
    basis = []
    for k in range(min(m, n)):
        best_value = -1.
        best = np.zeros(n)
        # 多起点避免某个固定初向量恰与最大本征空间正交。
        for start in range(n):
            v = np.eye(n)[:, start].copy()
            for q in basis:
                v -= q * (q @ v)
            norm = np.sqrt(np.sum(v*v))
            if norm < 1e-12:
                continue
            v /= norm
            for iteration in range(1500):
                w = B @ v
                for q in basis:
                    w -= q * (q @ w)
                norm = np.sqrt(np.sum(w*w))
                if norm < 1e-12:
                    break
                v = w / norm
                value = v @ B @ v
                if np.sqrt(np.sum((B@v-value*v)**2)) < 1e-13:
                    break
            value = v @ G @ v
            if value > best_value:
                best_value, best = value, v.copy()
        if best_value < 1e-10:
            break
        v = best
        sigma = np.sqrt(best_value * scale)
        u = A @ v / sigma
        singular[k] = sigma
        reconstructed += sigma * u[:, None] * v[None, :]
        inverse += v[:, None] * u[None, :] / sigma
        basis.append(v)
        B -= best_value * v[:, None] * v[None, :]
    return singular, reconstructed, inverse
'''
CASES = '''def make_cases(rng):
    cases = [(np.zeros((2, 3)),), (np.diag([3., 1.]),),
             (np.array([[1., 1.], [1., -1.], [0., 0.]]),),
             (np.array([[0., 2., 0.], [0., 0., 0.]]),)]
    for k in range(28):
        m, n = 2 + k % 2, 2 + (k//2) % 2
        U, _ = np.linalg.qr(rng.normal(size=(m, m)))
        V, _ = np.linalg.qr(rng.normal(size=(n, n)))
        values = np.array([3., 1., .4])[:min(m, n)].copy()
        if k % 3 == 0:
            values[-1] = 0.
        S = np.zeros((m, n))
        for j in range(len(values)):
            S[j, j] = values[j]
        cases.append((U @ S @ V.T,))
    return cases
'''
WRITTEN = [
    {"id": "ch9-w1", "level": "基础", "problem_md": r"对 $A=\begin{pmatrix}0&2\\1&0\end{pmatrix}$ 手算 SVD。", "solution_md": r"""
1. 依据乘法，$A^TA=\operatorname{diag}(1,4)$，降序特征值为 4、1。
2. 依据开方和本征方向，$\sigma_1=2,\sigma_2=1$，$v_1=e_2,v_2=e_1$。
3. 依据 $u_i=Av_i/\sigma_i$，得 $u_1=e_1,u_2=e_2$，因此 U=I，$V=\begin{pmatrix}0&1\\1&0\end{pmatrix}$。
4. 依据重构乘法，$U\operatorname{diag}(2,1)V^T=\begin{pmatrix}0&2\\1&0\end{pmatrix}=A$，两端基均正交归一。
""", "pitfalls_md": "奇异值是 A†A 本征值的平方根，并不直接等于 A 的本征值。"},
    {"id": "ch9-w2", "level": "基础", "problem_md": r"A 为 $\operatorname{diag}(5,3,1)$，求一个最佳秩一近似及两种误差。", "solution_md": r"""
1. 依据 A 为非负对角矩阵，其 SVD 可取 U=V=I，奇异值为 5、3、1。
2. 依据截断定义，只留第一个三元组，$A_1=\operatorname{diag}(5,0,0)$。
3. 依据酉不变范数公式，Frobenius 误差为 $\sqrt{3^2+1^2}=\sqrt{10}$，谱范数误差为 3。
4. 依据极小极大证明的下界，任意秩至多一的 B 均有 $\sigma_1(A-B)\ge3,\sigma_2(A-B)\ge1$，故 F 误差平方至少 10，A₁ 达到下界。
""", "pitfalls_md": "Frobenius 误差为被丢奇异值平方和的平方根，不是直接相加。"},
    {"id": "ch9-w3", "level": "进阶", "problem_md": r"拟合 $A=(1,1)^T,b=(1,3)^T$，由残差正交推正规方程并求解。", "solution_md": r"""
1. 依据列空间，任意 Ax=(x,x)ᵀ；残差为 $(1-x,3-x)^T$。
2. 依据最小点残差正交于列 $(1,1)^T$，有 $(1-x)+(3-x)=0$，即正规方程 2x=4。
3. 依据标量方程，x=2，拟合值 (2,2)ᵀ，残差 (−1,1)ᵀ，误差平方为 2。
4. 依据任意扰动 h 的展开，误差平方 $(1-2-h)^2+(3-2-h)^2=2+2h^2\ge2$，验证全局最小。
""", "pitfalls_md": "最小二乘不把两个互相矛盾的方程强行当作可精确求解的方程组。"},
    {"id": "ch9-w4", "level": "进阶", "problem_md": r"$A=\begin{pmatrix}1&1\\0&0\end{pmatrix},b=(2,1)^T$。求伪逆、所有最小二乘解与其中最小范数解。", "solution_md": r"""
1. 依据秩一 SVD，$\sigma_1=\sqrt2,u_1=e_1,v_1=(1,1)^T/\sqrt2$，所以 $A^+=v_1u_1^T/\sqrt2=\frac12\begin{pmatrix}1&0\\1&0\end{pmatrix}$。
2. 依据伪逆作用，$x_+=A^+b=(1,1)^T$，像为 (2,0)ᵀ，残差 (0,1)ᵀ 与列空间正交。
3. 依据核方程 $z_1+z_2=0$，所有最小二乘解为 $(1+t,1-t)^T$，t 为实数。
4. 依据正交分解，范数平方为 $(1+t)^2+(1-t)^2=2+2t^2$，仅在 t=0 最小，得到唯一最小范数解 (1,1)ᵀ。
""", "pitfalls_md": "A†A 奇异时不能求普通逆；仍有伪逆和最小二乘解。"},
    {"id": "ch9-w5", "level": "考试风格", "problem_md": r"证明任何秩至多 k 的 B 满足 $\|A-B\|_2\ge\sigma_{k+1}(A)$（0≤k<min(m,n)）。", "solution_md": r"""
1. 依据秩零度，$\dim\ker B\ge n-k$；取 $L=\operatorname{span}(v_1,\ldots,v_{k+1})$，维数为 k+1。
2. 依据交空间维数公式，$\dim(L\cap\ker B)\ge1$，可取其中一个单位向量 x。
3. 依据 Bx=0 与 x 只有前 k+1 个右奇异坐标，$\|(A-B)x\|^2=\|Ax\|^2=\sum_{i\le k+1}\sigma_i^2|c_i|^2\ge\sigma_{k+1}^2$。
4. 依据谱范数是所有单位向量上的最大伸缩，$\|A-B\|_2\ge\|(A-B)x\|\ge\sigma_{k+1}$，截断 SVD 达到下界。
""", "pitfalls_md": "x 必须同时位于前 k+1 个右奇异方向与 ker B；任取核向量不保证下界。"},
    {"id": "ch9-w6", "level": "考试风格", "problem_md": r"态 $\psi=(2|00\rangle+i|11\rangle)/\sqrt5$。由振幅矩阵的 SVD 求 Schmidt 系数与约化态本征值。", "solution_md": r"""
1. 依据振幅排布，$C=\operatorname{diag}(2,i)/\sqrt5$，$C^\dagger C=\operatorname{diag}(4,1)/5$。
2. 依据正平方根，Schmidt 系数为 $2/\sqrt5,1/\sqrt5$；可取 $U=\operatorname{diag}(1,i),V=I$。
3. 依据 SVD 系数展开，$\psi=(2/\sqrt5)|0\rangle|0\rangle+(1/\sqrt5)(i|1\rangle)|1\rangle$，局域 i 相位吸收进左基态。
4. 依据 $\rho_A=CC^\dagger$，本征值为 4/5、1/5，纯度为 $16/25+1/25=17/25$；两个非零 Schmidt 系数表明不是直积态。
""", "pitfalls_md": "Schmidt 系数非负，复相位放进局域基；约化态本征值是系数平方。"},
]
QUIZ = [
    {"q_md": "SVD 对哪些矩阵存在？", "options": ["仅实对称方阵", "仅可逆方阵", "任意有限实或复矩阵", "仅满列秩矩阵"], "answer": 2, "explain_md": "A†A 始终厄米半正定；零值分支用正交补，不要求矩形 A 可逆。"},
    {"q_md": r"$A^\dagger A$ 的特征值为 9、4，奇异值是？", "options": ["9、4", "3、2", "81、16", "−3、−2"], "answer": 1, "explain_md": "奇异值按定义是半正定 Gram 矩阵特征值的非负平方根。"},
    {"q_md": "亏秩 A 的伪逆最小二乘解还满足什么？", "options": ["所有解中范数最小", "总能使残差为零", "一定有普通逆", "每个坐标为正"], "answer": 0, "explain_md": "伪逆把核坐标设零；其他最小点添加正交核向量，只会增大范数平方。"},
    {"q_md": "正规方程的几何依据是？", "options": ["残差与 b 平行", "残差长度必须为零", "A 必须方阵", "残差与 A 的列空间正交"], "answer": 3, "explain_md": "沿任意列方向扰动都不能改善误差，交叉项必须为零，即 A†(b−Ax)=0。"},
    {"q_md": "归一化纯态的 Schmidt 系数平方之和为？", "options": ["矩阵秩", "1", "局域维数之和", "可以为负"], "answer": 1, "explain_md": "两端 Schmidt 基正交归一，态范数平方化为系数平方和。"},
]
for _i, _q in enumerate(QUIZ, 1):
    _q["id"] = f"ch9-q{_i}"

CHAPTER = {
    "id": "ch9", "num": 9, "title": "奇异值分解与低秩近似", "depends_on": ["ch5", "ch6", "ch8"],
    "theory_md": THEORY, "viz": ["svd_ellipse"], "written": WRITTEN, "quiz": QUIZ,
    "code_exercises": [{
        "id": "ch9-c1", "title": "幂迭代与收缩构造小型实矩阵 SVD",
        "prompt_md": r"实现 `svd_from_gram(A)`，A 为 m×n 实矩阵，m、n 各为 2 或 3。返回 `(singular, reconstructed, inverse)`：长度 min(m,n) 的降序奇异值（零补尾）、由求得的奇异三元组重构的 m×n 矩阵、由相同三元组构造的 n×m 伪逆。必须形成 G=A.T@A，用幂迭代+秩一收缩求正本征对，再令 σ=√λ、u=Av/σ、累加 σuvᵀ 与 vuᵀ/σ；禁止直接返回输入充当重构。允许 @、sqrt、逐项求和，不用 linalg（含 svd/eig/eigh/eigvalsh/qr/pinv/norm）、roots、poly、scipy 等。输入为零矩阵或非零奇异值在 [0.4,3] 内、相邻不等正值的比值至多 0.8 的矩阵，允许相等奇异值及秩亏。建议先用 Frobenius 范数平方归一化 G，归一化特征值≤1e−10 视为零；用多个坐标起点避免恰好正交，每次收缩后扣除已选基分量。零矩阵单列返回全零；只需正奇异方向的精简分解，无需返回零方向的任意扩张基。",
        "function_name": "svd_from_gram", "starter": "def svd_from_gram(A):\n    # TODO：Gram -> 多起点幂迭代+收缩 -> 正奇异三元组 -> 重构与伪逆\n    m, n = A.shape\n    return np.zeros(min(m,n)), np.zeros((m,n)), np.zeros((n,m))\n",
        "banned": ["np.linalg", "np.linalg.svd", "np.linalg.eig", "np.linalg.eigh", "np.linalg.eigvals", "np.linalg.eigvalsh", "np.linalg.pinv", "np.poly", "np.roots", "np.polynomial", "np.matrix", "np.asmatrix", "np.mat", "np.bmat", "np.matrixlib", "np.matlib", "scipy"],
        "reference": REFERENCE, "make_cases": CASES,
    }],
}
