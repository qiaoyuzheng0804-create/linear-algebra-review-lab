"""公式、诊断、知识依赖与闪卡数据；无 Streamlit 依赖。"""
from collections import Counter

from .content import CHAPTERS, CHAPTER_TITLES

# 每项为（LaTeX，成立条件及对应推导）；与十章正文使用相同的内积和基序约定。
FORMULA_SHEET = {
    "地基": [
        (r"e^{i\theta}=\cos\theta+i\sin\theta", "实数角 θ；由绝对收敛幂级数分组得到，见地基 F2"),
        (r"\operatorname{proj}_v u=\frac{u\cdot v}{v\cdot v}v", "实内积空间，v 非零；见地基 F3"),
        (r"\|u\times v\|=\|u\|\|v\||\sin\theta|", "定向三维欧氏空间，夹角式要求 u、v 非零；见地基 F4"),
        (r"[v]_E=P[v]_B", "P 的列为新基向量的旧坐标；见地基 F5"),
    ],
    "第1章": [
        (r"(-1)v=-v,\quad 0v=0_V", "任意向量空间；推导 1.1"),
        (r"(v_1,\ldots,v_m)\text{相关}\iff\exists c\ne0:\sum_jc_jv_j=0", "固定数域上的有限向量组；推导 1.2"),
        (r"\operatorname{rank}(v_1,\ldots,v_m)=|\text{极大无关组}|", "有限组；交换引理与推导 1.4 保证长度唯一"),
    ],
    "第2章": [
        (r"[T]_{C\leftarrow B}=([T(b_1)]_C,\ldots,[T(b_n)]_C)", "T 线性，B、C 分别为定义域和陪域的基；推导 2.1"),
        (r"(AD)_{ij}=\sum_k A_{ik}D_{kj},\quad [T\circ S]=[T][S]", "维数与中间坐标基相容；推导 2.2"),
        (r"A_B=P^{-1}A_EP", "同一空间同一算符，两端同换基且 P 可逆；式 (2.2)"),
        (r"A_{\rm new}=Q^{-1}A_{\rm old}P", "两端分别用 P、Q 换基，二者可逆；式 (2.3)"),
    ],
    "第3章": [
        (r"Ax=b\text{有解}\iff\operatorname{rank}A=\operatorname{rank}[A\mid b]", "有限矩阵，固定数域；式 (3.1)"),
        (r"Ax=b\text{唯一解}\iff\operatorname{rank}A=\operatorname{rank}[A\mid b]=n", "n 为未知量数；式 (3.2)"),
        (r"\dim\ker A=n-\operatorname{rank}A", "A 有 n 列；式 (3.3)"),
        (r"\{x:Ax=b\}=x_p+\ker A", "已存在特解 x_p；非齐次时是仿射集"),
    ],
    "第4章": [
        (r"\det A=\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)\prod_{i=1}^n a_{i,\sigma(i)}", "A 为 n 阶方阵；Leibniz 定义，式 (4.1)"),
        (r"\det A=\sum_j a_{ij}C_{ij}", "固定第 i 行，C 为代数余子式；式 (4.2)"),
        (r"\det(AB)=\det A\det B", "A、B 为同阶方阵；初等矩阵分解证明，式 (4.4)"),
        (r"A\text{可逆}\iff\det A\ne0", "方阵；伴随恒等式与乘法公式，推导 4.4"),
        (r"x_j=\frac{\det A_j(b)}{\det A}", "Ax=b，det A 非零；A_j(b) 将第 j 列替换为 b，式 (4.5)"),
        (r"A\operatorname{adj}A=\operatorname{adj}A\,A=(\det A)I", "任意有限方阵，包括奇异矩阵；第4章余子式展开推导"),
        (r"\mathrm{Vol}(A\Omega)=|\det A|\,\mathrm{Vol}(\Omega)", "实方阵作用于可测区域，使用普通体积；第4章几何直觉"),
    ],
    "第5章": [
        (r"Av=\lambda v\ (v\ne0)\iff (A-\lambda I)v=0", "有限维方阵；λ 在所选数域中，v 非零；推导 5.1"),
        (r"\det(A-\lambda I)=0", "λ 为特征值的充要条件，A 为方阵；式 (5.2)"),
        (r"1\le\dim\ker(A-\lambda I)\le a_\lambda", "λ 为所选数域内的特征值，a_λ 为代数重数；推导 5.3"),
        (r"p_A(A)=0,\quad p_A(z)=\det(zI-A)", "任意有限方阵；Cayley–Hamilton，式 (5.3)"),
        (r"E_{\lambda_1}+\cdots+E_{\lambda_s}=E_{\lambda_1}\oplus\cdots\oplus E_{\lambda_s}", "互异特征值的本征空间；第5章直和证明，不要求矩阵可逆"),
    ],
    "第6章": [
        (r"\langle u,v\rangle=u^\dagger v", "复欧氏空间，第一变量共轭线性；定义 6.1"),
        (r"|\langle x,y\rangle|^2\le\|x\|^2\|y\|^2", "正定内积，等号当且仅当两向量线性相关；式 (6.1)"),
        (r"u_k=v_k-\sum_{j<k}\langle q_j,v_k\rangle q_j,\quad q_k=u_k/\|u_k\|", "q_j 正交归一，输入组无关确保 u_k 非零；式 (6.2)"),
        (r"A=QR,\quad Q^\dagger Q=I", "A 满列秩的薄 QR，R 上三角且对角线正实；推导 6.3"),
        (r"p=QQ^\dagger x,\quad \|x-y\|^2=\|x-p\|^2+\|p-y\|^2", "Q 的列为 W 的正交归一基，y∈W；式 (6.3)～(6.4)"),
        (r"A^\dagger A\hat c=A^\dagger b", "有限维实或复矩阵，最小二乘正规方程；推导 6.4"),
        (r"\hat c=(A^\dagger A)^{-1}A^\dagger b", "A 满列秩，最小二乘系数唯一；推导 6.4 的正定性论证"),
    ],
    "第7章": [
        (r"A=P\Lambda P^{-1}\iff A\text{有一组特征向量基}", "在所选数域上；重根不自动保证可对角化，式 (7.2)"),
        (r"\det B=\det A,\quad\operatorname{tr}B=\operatorname{tr}A,\quad p_B=p_A", "B=P⁻¹AP，P 可逆；推导 7.2"),
        (r"A^k=P\Lambda^kP^{-1}", "A 可对角化；k 为非负整数，负整数还要求 A 可逆，式 (7.6)"),
        (r"e^{tA}=\sum_{k=0}^{\infty}\frac{t^kA^k}{k!}", "任意有限方阵，有限 t；绝对收敛，式 (7.1)"),
        (r"e^{t(A+B)}=e^{tA}e^{tB}", "充分条件 AB=BA；不对易时二项式合项步骤失效，式 (7.5)"),
        (r"\dot x=Ax,\quad x(t)=e^{tA}x(0)", "常系数有限维系统，初值设在 t=0；式 (7.7)"),
        (r"U(t)=e^{-iHt/\hbar},\quad U(t)^\dagger U(t)=I", "有限维、时间无关厄米 H，实 t、非零实 ħ；第7章 QM 连接"),
    ],
    "第8章": [
        (r"H=H^\dagger\Longrightarrow H=U\Lambda U^\dagger,\quad\Lambda\text{为实对角阵}", "有限维厄米矩阵，U 酉；谱定理"),
        (r"A^\dagger A=AA^\dagger\iff A=U\Lambda U^\dagger", "复方阵；正规矩阵的酉对角化判据"),
        (r"U^\dagger U=I\Longrightarrow\|U\psi\|=\|\psi\|", "有限维酉方阵，标准复内积"),
        (r"H=\sum_\lambda\lambda P_\lambda,\quad P_\lambda P_\mu=\delta_{\lambda\mu}P_\lambda", "有限维厄米矩阵，P_λ 为整个本征空间的正交投影"),
        (r"Hv=\lambda v\Longrightarrow\lambda\in\mathbb R,\quad Hx=\lambda x,\ Hy=\mu y,\ \lambda\ne\mu\Longrightarrow x^\dagger y=0", "H 厄米，本征向量非零；推导 8.1、8.2；实谱不意味着正谱"),
        (r"U^\dagger U=I\iff\langle Ux,Uy\rangle=\langle x,y\rangle\ (\forall x,y)\iff\|Ux\|=\|x\|\ (\forall x)", "U 为复线性方阵；极化恒等式闭合等价链，式 (8.4)"),
        (r"\langle x,y\rangle=\frac14(\|x+y\|^2-\|x-y\|^2-i\|x+iy\|^2+i\|x-iy\|^2)", "第一变量共轭线性的复内积；式 (8.3)，虚部符号与内积约定相关"),
        (r"p(\lambda)=\|P_\lambda\psi\|^2=\langle\psi|P_\lambda|\psi\rangle,\quad\sum_\lambda P_\lambda=I", "归一化纯态、厄米可观测量；简并时投影须覆盖整个本征空间；式 (8.5)"),
    ],
    "第9章": [
        (r"A=U\Sigma V^\dagger", "任意有限实或复 m×n 矩阵；U、V 酉，奇异值非负"),
        (r"A^+=V\Sigma^+U^\dagger", "任意有限矩阵；Σ⁺ 将非零奇异值取倒数，零保留"),
        (r"\|A-A_k\|_2=\sigma_{k+1},\quad \|A-A_k\|_F^2=\sum_{j>k}\sigma_j^2", "奇异值降序，A_k 为截断 SVD，0≤k<min(m,n)；最优秩至多 k 近似"),
        (r"A^\dagger Av_i=\sigma_i^2v_i,\quad u_i=Av_i/\sigma_i", "后式只用于正奇异值；零值处改用正交补扩张；推导 9.1"),
        (r"\sigma_j(A)=\min_{\dim S=n-j+1}\max_{x\in S,\|x\|=1}\|Ax\|", "A 有 n 列，1≤j≤n，超过 min(m,n) 的奇异值补零；式 (9.2)"),
        (r"A^\dagger A x=A^\dagger b\iff b-Ax\perp\operatorname{im}A", "任意有限实或复 A，等价于 x 是最小二乘解；式 (9.4)"),
        (r"x=A^+b+z,\quad z\in\ker A", "所有最小二乘解；A⁺b 为唯一最小范数者；式 (9.5)"),
        (r"A^+b=(A^\dagger A)^{-1}A^\dagger b", "仅当 A 满列秩；亏秩时该求逆步骤失效，式 (9.6)"),
    ],
    "第10章": [
        (r"\dim(V\otimes W)=\dim V\,\dim W", "有限维同数域向量空间"),
        (r"(A\otimes B)(C\otimes D)=AC\otimes BD", "AC 与 BD 的维数分别相容"),
        (r"(A\otimes B)^\dagger=A^\dagger\otimes B^\dagger", "有限实或复矩阵"),
        (r"\operatorname{Tr}_B\rho=\sum_j(I\otimes\langle j|)\rho(I\otimes|j\rangle)", "有限维复合系统，|j⟩ 是 B 的任一正交归一基"),
        (r"(\rho_A)_{ik}=\sum_j\rho_{(i,j),(k,j)}", "联合基序为 i·dimB+j，只消去 B 指标；式 (10.5)"),
        (r"|\psi\rangle=\sum_{\alpha=1}^r s_\alpha|u_\alpha\rangle\otimes|w_\alpha\rangle,\quad s_\alpha>0,\quad\sum_\alpha s_\alpha^2=1", "Schmidt 分解：归一化双系统纯态，两端基各自正交归一；由振幅矩阵 SVD 构造，式 (10.6)"),
        (r"C=U\Sigma V^\dagger,\quad w_\alpha=\bar v_\alpha,\quad\rho_A=CC^\dagger,\quad\rho_B=C^T\bar C", "Cij 是 |i⟩A|j⟩B 的振幅；复系数时 B 端共轭不能省略，式 (10.7)"),
        (r"\operatorname{tr}(\rho_A^2)=\sum_\alpha s_\alpha^4=1\iff r=1\iff\psi\text{为直积态}", "联合态已知为纯态；联合混合态不能由约化纯度单独判纠缠，推导 10.5"),
    ],
}

KNOWLEDGE_MAP = {
    "foundation": [],
    "ch1": ["foundation"],
    "ch2": ["foundation", "ch1"],
    "ch3": ["foundation", "ch1", "ch2"],
    "ch4": ["foundation", "ch2", "ch3"],
    "ch5": ["foundation", "ch3", "ch4"],
    "ch6": ["foundation", "ch1", "ch2", "ch3"],
    "ch7": ["foundation", "ch4", "ch5", "ch6"],
    "ch8": ["foundation", "ch5", "ch6"],
    "ch9": ["foundation", "ch5", "ch6", "ch8"],
    "ch10": ["foundation", "ch2", "ch8", "ch9"],
}


def chapter_label(chapter):
    return "地基" if chapter == "foundation" else f"第{int(chapter[2:])}章：{CHAPTER_TITLES[int(chapter[2:])]}"


# 诊断结果按 id→bool 保存；未答题不计作答错。
DIAGNOSTIC = [
    {"chapter": "foundation", "q_md": r"$i^2$ 等于？", "options": ["1", "0", "-1", "i"], "answer": 2, "explain_md": "依地基 F1，虚数单位按平方等于 -1 定义。"},
    {"chapter": "foundation", "q_md": r"非零实向量正交等价于？", "options": ["点积为零", "叉积为零", "长度相等", "坐标相等"], "answer": 0, "explain_md": r"由 F3，点积为模长乘积乘夹角余弦，非零时为零等价于夹角为 $\pi/2$。"},
    {"chapter": "foundation", "q_md": r"$T(x)=x+2$ 为什么不是线性映射？", "options": ["没有定义域", "不保持零", "不能求逆", "不是满射"], "answer": 1, "explain_md": r"由 F7，线性映射必须 $T(0)=0$，这里 $T(0)=2$。"},
    {"chapter": "ch1", "q_md": "三个向量中一个是其余两个之和，则该组？", "options": ["必无关", "必是基", "秩必为 3", "必相关"], "answer": 3, "explain_md": r"将等式移项得到不全为零的零组合，依据推导 1.2 相关。"},
    {"chapter": "ch2", "q_md": r"若 $T(e_1)=(1,0)^T,T(e_2)=(2,1)^T$，矩阵为？", "options": [r"$\begin{pmatrix}1&2\\0&1\end{pmatrix}$", r"$\begin{pmatrix}1&0\\2&1\end{pmatrix}$", r"$I$", r"$0$"], "answer": 0, "explain_md": "依推导 2.1，把两个像排成列。"},
    {"chapter": "ch3", "q_md": "4 个未知量、秩为 2 的齐次方程组，核维数为？", "options": ["0", "1", "2", "4"], "answer": 2, "explain_md": r"式 (3.3) 给 $4-2=2$。"},
    {"chapter": "ch3", "q_md": "消元出现 0=1 的行意味着？", "options": ["唯一解", "无解", "无穷多解", "可忽略该行"], "answer": 1, "explain_md": "依据推导 3.3，增广列增加主元，秩判据不满足。"},
    {"chapter": "ch4", "q_md": "二维变换的普通面积缩放比为？", "options": [r"$\det A$", r"$\operatorname{tr}A$", r"$|\det A|$", r"$\det A+1$"], "answer": 2, "explain_md": "普通面积非负，因此面积比是行列式的绝对值，符号用于定向。"},
    {"chapter": "ch5", "q_md": r"λ 是 n 阶矩阵 A 的特征值的充要条件是？", "options": [r"$\det A\ne0$", r"$\lambda\ne0$", r"$\operatorname{tr}A=\lambda$", r"$\det(A-\lambda I)=0$"], "answer": 3, "explain_md": "行列式为零等价于 A−λI 有非零核向量，也就是存在非零本征态。"},
    {"chapter": "ch6", "q_md": "Gram–Schmidt 中某步残差为零说明？", "options": ["当前向量在前面向量的张成内", "已得新单位向量", "矩阵必为方阵", "内积失效"], "answer": 0, "explain_md": "扣除已有方向投影后没有剩余分量，因此归一化分母为零，应识别相关性。"},
    {"chapter": "ch9", "q_md": "SVD 对哪些有限矩阵存在？", "options": ["仅可逆方阵", "任意实或复矩阵", "仅厄米矩阵", "仅对角矩阵"], "answer": 1, "explain_md": "推导 9.1 从始终厄米半正定的 A†A 构造；零奇异值处以正交补扩张，左右可使用不同基，所以不要求方阵或满秩。"},
    {"chapter": "ch8", "q_md": "厄米矩阵的特征值必为？", "options": ["正数", "整数", "实数", "纯虚数"], "answer": 2, "explain_md": r"若 $Hv=\lambda v$，则 $\lambda=v^\dagger Hv/(v^\dagger v)$，分子因 H 厄米等于自身共轭，故 λ 实。"},
    {"chapter": "ch7", "q_md": r"常系数系统 $\dot x=Ax$ 的传播矩阵为？", "options": [r"$tA$", r"$A^{-1}$", r"$A+tI$", r"$e^{tA}$"], "answer": 3, "explain_md": r"指数级数逐项求导给 $d e^{tA}/dt=Ae^{tA}$，且在 t=0 等于 I。"},
    {"chapter": "ch10", "q_md": "二维空间与三维空间的张量积维数为？", "options": ["5", "6", "3", "2"], "answer": 1, "explain_md": "推导 10.1 构造 2×3 个独立坐标槽位，直积基无关且张成，所以维数为 6；不是直和的 2+3。"},
]
for _i, _q in enumerate(DIAGNOSTIC, 1):
    _q["id"] = f"diag-{_i}"


def recommend(results):
    """先补地基；否则定位错题最多章的最早出错前置章。"""
    answered = {q["id"]: bool(results[q["id"]]) for q in DIAGNOSTIC if q["id"] in results}
    if not answered:
        return "尚未完成诊断，请先提交作答。"
    wrong = Counter(q["chapter"] for q in DIAGNOSTIC if q["id"] in answered and not answered[q["id"]])
    prefix = f"已答 {len(answered)}/{len(DIAGNOSTIC)} 题。"
    if not wrong:
        return prefix + "已答题全部正确，建议从第1章按学习流程巩固；尚未答题请继续诊断。" if len(answered) < len(DIAGNOSTIC) else prefix + "全部正确，建议从第1章快速复核推导，再完成手写代码练习。"
    if wrong["foundation"]:
        start = "foundation"
    else:
        weakest = min(wrong, key=lambda c: (-wrong[c], int(c[2:])))
        ancestors = set()
        def visit(node):
            for dependency in KNOWLEDGE_MAP[node]:
                if dependency not in ancestors:
                    ancestors.add(dependency)
                    visit(dependency)
        visit(weakest)
        candidates = [c for c in wrong if c == weakest or c in ancestors]
        start = min(candidates, key=lambda c: int(c[2:]))
    distribution = "、".join(f"{chapter_label(c).split('：')[0]} {wrong[c]} 题" for c in KNOWLEDGE_MAP if wrong[c])
    result = prefix + f"错题分布：{distribution}。建议从{chapter_label(start)}开始，先补前置概念。"
    if start not in {"foundation", *(c["id"] for c in CHAPTERS.values())}:
        result += "该章建设中，可先查公式手册并复核已开放的前置章。"
    return result


_CARD_DATA = [
    ("foundation", r"复数模如何定义？", r"$|a+ib|=\sqrt{a^2+b^2}=\sqrt{z\bar z}$，见 F1。"),
    ("foundation", r"欧拉公式如何从级数出现？", r"将 $e^{i\theta}$ 的绝对收敛级数分为偶次和奇次，分别得到 $\cos\theta$ 与 $i\sin\theta$，见 F2。"),
    ("foundation", r"复数相乘如何改变相位？", r"$re^{i\theta}se^{i\phi}=rs e^{i(\theta+\phi)}$，模相乘、相位相加。"),
    ("foundation", r"点积的几何意义是什么？", r"$u\cdot v=\|u\|\|v\|\cos\theta$；比较坐标平方展开与余弦定理，见 F3。"),
    ("foundation", r"投影为什么要除以 $v\cdot v$？", r"令残差 $(u-cv)\cdot v=0$，解得 $c=(u\cdot v)/(v\cdot v)$，要求 $v\ne0$。"),
    ("foundation", r"叉积的长度表达什么？", r"$\|u\times v\|$ 等于两向量张成的平行四边形面积，见 F4。"),
    ("foundation", r"新基列组成 P，坐标如何变？", r"$[v]_E=P[v]_B$，因为 $v=\sum_jx_jb_j$；反向要求 P 可逆。"),
    ("foundation", r"单射与满射分别约束什么？", r"单射要求相同输出只能来自相同输入；满射要求指定陪域的每个元素都能达到，见 F6。"),
    ("foundation", r"保持零能否证明线性？", r"不能。线性必须保持任意线性组合；$x^2$ 保持零却不满足齐次性，见 F7。"),
    ("ch1", r"零向量唯一性的证明关键是什么？", r"若 0、0′ 都是零元，则 $0+0'=0$ 且 $0+0'=0'$，故相等，见推导 1.1。"),
    ("ch1", r"如何证明 $(-1)v=-v$？", r"分配律给 $v+(-1)v=(1-1)v=0$，由加法逆元唯一性得到结论。"),
    ("ch1", r"线性相关的计算判据是什么？", r"存在 $c\ne0$ 使 $Ac=0$，其中 A 的列为待检向量，见推导 1.2。"),
    ("ch1", r"非平凡零组合如何给冗余向量？", r"选非零系数 $c_k$，解出 $v_k=-\sum_{j\ne k}(c_j/c_k)v_j$。"),
    ("ch1", r"极大无关组为什么张成原组？", r"若有向量不在其张成内，加入后仍无关，违反极大性，见推导 1.4。"),
    ("ch1", r"Steinitz 引理每步为什么能替换？", r"新无关向量不在先前 u 的张成中，所以对剩余 w 的展开至少一个系数非零，可解出并替换该 w。"),
    ("ch1", r"基长度为什么唯一？", r"以一组为无关组、另一组为张成组应用交换引理，两方向分别给出长度不等式。"),
    ("ch1", r"极大无关组本身是否唯一？", r"不唯一，长度唯一。不同主元或选取顺序可得到不同基。"),
    ("ch1", r"归一化态集合为什么不构成向量空间？", r"乘以零后不再归一化，违反标量乘法封闭性；完整态向量空间才满足叠加运算。"),
    ("ch2", r"线性映射由什么有限数据唯一确定？", r"定义域的一组基向量的像；对任意输入按基展开后用线性性，见推导 2.1。"),
    ("ch2", r"矩阵第 j 列是什么？", r"第 j 个输入基向量之像在输出基下的坐标。"),
    ("ch2", r"矩阵乘法为什么是行列求和？", r"对输入第 j 个基先 S 后 T，逐层展开系数得到 $\sum_kA_{ik}D_{kj}$，见推导 2.2。"),
    ("ch2", r"AD 中谁先作用？", r"D 先作用，复合映射为 $T\circ S$；中间坐标基须相容。"),
    ("ch2", r"相似关系如何写？", r"$A_B=P^{-1}A_EP$，P 的列是新基的旧坐标；见式 (2.2)。"),
    ("ch2", r"两端独立换基如何写？", r"$A_{\rm new}=Q^{-1}AP$，输入用 P，输出用 Q，见式 (2.3)。"),
    ("ch2", r"何时逆矩阵等于共轭转置？", r"换基矩阵酉时 $U^{-1}=U^\dagger$；一般可逆矩阵没有此性质。"),
    ("ch2", r"算符矩阵元 $A_{ij}$ 的意义？", r"正交归一基中 $A_{ij}=\langle e_i|\hat A|e_j\rangle$，是基态像的展开振幅。"),
    ("ch3", r"行操作保持解集的核心理由？", r"每种操作都可逆：交换再交换、乘非零数后乘倒数、加行倍数后减回，见推导 3.1。"),
    ("ch3", r"为什么不能只操作系数矩阵？", r"等式两边必须同步，否则不再是原方程的等价组合。"),
    ("ch3", r"有解的秩判据是什么？", r"$\operatorname{rank}A=\operatorname{rank}[A\mid b]$，即 b 在 A 的列空间，见式 (3.1)。"),
    ("ch3", r"唯一解还要满足什么？", r"有解且系数秩等于未知量数 n，排除所有自由变量，见式 (3.2)。"),
    ("ch3", r"如何构造零空间基？", r"每次令一个自由变量为 1、其余为 0，按主元方程解出主元变量；这些解在自由坐标上组成单位列。"),
    ("ch3", r"零度为什么为 n−r？", r"上述 n−r 个核向量张成所有齐次解，且查看自由坐标可证明无关，见推导 3.4。"),
    ("ch3", r"非齐次全部解如何表示？", r"$x=x_p+z$，$z\in\ker A$；相减两解得核向量，反向相加可验证，见推导 3.4。"),
    ("ch3", r"部分主元选择什么？", r"当前列尚未使用的行中绝对值最大的元素，减少除以小数造成的误差。"),
    ("ch3", r"本征方程何时存在非零态？", r"$\operatorname{rank}(H-EI)<n$，使核维数正；零解不代表物理态。"),
    ("ch3", r"非齐次解集为什么通常不是子空间？", r"$b\ne0$ 时零向量不满足 Ax=b；有解时它是核的平移，属于仿射集。"),
    ("ch4", r"行列式为何换两行变号？", r"置换求和中交换两个位置使逆序数奇偶翻转，乘积条目重新配对，每项恰变号。"),
    ("ch4", r"把一行加另一行的倍数，行列式为何不变？", r"按行多线性展开为原行列式加两行相同的项；后者由交错性为零。"),
    ("ch4", r"det(AB) 的乘法公式何时成立？", r"A、B 为同阶方阵；用初等矩阵分解和各行操作对 det 的影响，见推导 4.3。"),
    ("ch4", r"奇异矩阵的伴随恒等式还成立吗？", r"成立：$A\operatorname{adj}A=(\det A)I$ 源于余子式展开；只有进一步除以 det A 求逆时要求非零。"),
    ("ch4", r"克拉默法则的分母为何不能为零？", r"推导要由 $x_j\det A=\det A_j(b)$ 除以 det A；为零时不能由此确定唯一系数。"),
    ("ch4", r"二维普通面积因子与定向因子有何区别？", r"普通面积比为 $|\det A|$，有向面积因子为 det A；负号记录翻转。"),
    ("ch5", r"特征多项式的根为何给非零本征向量？", r"根使 det(λI−A)=0，故 λI−A 不可逆、核维数为正，可取非零核向量。"),
    ("ch5", r"代数重数与几何重数分别数什么？", r"前者数特征多项式根的重数，后者数对应核空间维数；$1\le g_\lambda\le a_\lambda$。"),
    ("ch5", r"不同本征值的本征空间为什么是直和？", r"若各空间向量之和为零，逐次施加 $(A-\lambda I)$ 消去某项，以互异差值非零推出每项为零。"),
    ("ch5", r"Cayley–Hamilton 说什么？", r"任意有限方阵满足 $p_A(A)=0$；推导 5.4 先在可对角化情形证明，再用复三角化与对角微扰构造互异谱逼近，最后以矩阵条目的多项式连续性取极限。"),
    ("ch5", r"实旋转矩阵为什么可能没有实特征向量？", r"非零非平角旋转不保留实直线；特征多项式可能只有非实根，改到复数域才可取本征向量。"),
    ("ch5", r"二阶 Jordan 块为何只有一个本征方向？", r"$J=\begin{pmatrix}\lambda&1\\0&\lambda\end{pmatrix}$ 的核方程给 x₂=0，所以本征空间一维而代数重数二。"),
    ("ch6", r"复内积的第一变量有什么约定？", r"第一变量共轭线性：$\langle cx,y\rangle=\bar c\langle x,y\rangle$；第二变量线性。"),
    ("ch6", r"Cauchy–Schwarz 的等号条件是什么？", r"两向量线性相关（含零向量）；正定性把最小二次式为零转成一个向量为零。"),
    ("ch6", r"Gram–Schmidt 残差为何与已有方向垂直？", r"$u=v-\sum q_j\langle q_j,v\rangle$，取与 qᵢ 的内积后利用 δᵢⱼ 恰抵消原分量。"),
    ("ch6", r"QR 的 R 元素从何而来？", r"原列在正交基下的系数 $r_{ij}=q_i^\dagger a_j$；逐列正交化给上三角，满列秩确保对角残差范数为正。"),
    ("ch6", r"正交投影为什么给最佳逼近？", r"$x-y=(x-p)+(p-y)$ 两项正交，误差平方分成固定项与非负项，只有 y=p 时第二项为零。"),
    ("ch6", r"薄 QR 的 Q 是方阵吗？", r"不一定。m×n 的满列秩矩阵给 m×n 的 Q，仅 $Q^\dagger Q=I_n$；m>n 时 $QQ^\dagger$ 为投影而不是 Iₘ。"),
    ("ch7", r"对角化的真正充要条件是什么？", r"所选数域上本征向量张满全空间；排成 P 可逆，逐列 AP=PΛ 才能乘 P⁻¹。"),
    ("ch7", r"重特征值是否必定不能对角化？", r"否。λI 的本征空间是全空间；问题在于每个几何重数是否达到代数重数。"),
    ("ch7", r"矩阵指数为何总存在？", r"次乘性范数给 $\|t^kA^k/k!\|\le |t|^k\|A\|^k/k!$，总和被标量指数控制，故绝对收敛。"),
    ("ch7", r"什么时候可以拆指数 e^(A+B)？", r"AB=BA 是充分条件，使二项式合项合法；一般非交换时二阶交叉项已经不相同。"),
    ("ch7", r"N²=0 的 Jordan 修正如何进入指数？", r"$e^{tN}=I+tN$；若 J=λI+N，二者对易，$e^{tJ}=e^{\lambda t}(I+tN)$。"),
    ("ch7", r"离散与连续稳定判据分别是什么？", r"模式 λ 的离散幅度为 |λ|ᵏ，连续幅度为 $e^{t\operatorname{Re}\lambda}$；衰减分别要求 |λ|<1 与 Re λ<0。"),
    ("ch8", r"厄米矩阵实谱证明的除法发生在哪？", r"先由共轭得 $(\lambda-\bar\lambda)\|x\|^2=0$，再由本征向量非零及正定性除以正范数平方。"),
    ("ch8", r"谱定理归纳为何能限制到正交补？", r"$q^\dagger Hw=(Hq)^\dagger w=\bar\lambda q^\dagger w=0$ 保证不变；限制矩阵元仍满足共轭对称。"),
    ("ch8", r"至少一个复特征向量的存在来自哪里？", r"代数学基本定理给首一特征多项式一个复根，行列式为零再给非零核；不是先假设谱定理。"),
    ("ch8", r"保范数推出保内积用什么？", r"极化恒等式及复线性，把 x±y、x±iy 的范数差恢复为内积实部与虚部。"),
    ("ch8", r"Schur 三角阵为何在正规时对角？", r"首行平方和与首列平方和相等，首列只有对角元，迫使首行其余元为零；对剩余正规块归纳。"),
    ("ch8", r"简并读数的测量概率怎么求？", r"$p(\lambda)=\|P_\lambda\psi\|^2$，Pλ 必须投影整个本征子空间，不能只取其中一个方向。"),
    ("ch8", r"酉演化保持哪些量？", r"保持任意两态内积与态范数；在固定测量基中的各结果概率仍可因相对相位变化而改变。"),
    ("ch8", r"非正规矩阵的谱证明会在哪断掉？", r"正交补未必不变；Schur 路线仍成立，但缺少行列平方和相等的正规性条件，超对角项不能消去。"),
    ("ch9", r"SVD 为何对矩形矩阵存在？", r"A†A 始终为输入空间上的厄米半正定方阵；输出方向另由 Avᵢ/σᵢ 构造，两端不要求同一基。"),
    ("ch9", r"零奇异值对应的 u=Av/σ 怎么办？", r"不能除以零；先由 ‖Av‖²=0 得 Av=0，非零像方向构造完后在输出正交补扩张。"),
    ("ch9", r"截断 SVD 的 F 误差是什么？", r"$\sqrt{\sum_{i>k}\sigma_i^2}$；酉不变性把误差化为丢弃的对角项，极小极大刻画证明其全局最优。"),
    ("ch9", r"正规方程的几何内容是什么？", r"$A^\dagger(b-Ax)=0$，残差与整个列空间正交；沿任何列方向移动都不能减少误差。"),
    ("ch9", r"伪逆为何选出最小范数解？", r"它将所有核方向系数设为零；其他最小二乘解为 x₊+z，二者正交，范数平方额外增加 ‖z‖²。"),
    ("ch9", r"Gram 路径求 SVD 的数值代价是什么？", r"A†A 将非零奇异值的最大最小比平方，可能放大病态性；本章实现限小规模适中谱，检验构造原理。"),
    ("ch10", r"张量积的基为什么有 mn 个？", r"每对局域基占一个独立系数槽位；这些符号张成所有系数表且逐槽比较证明无关。"),
    ("ch10", r"张量积算符的复合规则是什么？", r"$(A\otimes B)(C\otimes D)=AC\otimes BD$；先在 eᵢ⊗fⱼ 上两端分别复合，再按线性性延拓。"),
    ("ch10", r"Bell 态不可分解的系数矛盾是什么？", r"直积要求 ac=bd=1/√2 且 ad=bc=0；前两非零使 abcd=1/2，后两又使 abcd=0。"),
    ("ch10", r"部分迹怎样保留 A 的条目？", r"$(\rho_A)_{ik}=\sum_j\rho_{(i,j),(k,j)}$；保留 i、k，对同一 B 指标 j 求和。"),
    ("ch10", r"Schmidt 分解与 SVD 的复共轭如何对应？", r"C=UΣV† 时，A 端基是 U 列，B 端基是 V 列的共轭，非负 Schmidt 系数是 C 的奇异值。"),
    ("ch10", r"什么时候可以用约化纯度判纠缠？", r"联合态已知为纯态时，约化纯度为 1 等价于 Schmidt 秩一；联合混合态的局域混合性本身不能判纠缠。"),
    ("ch10", r"Bell 态与经典相关混合态能否有相同局域矩阵？", r"可以，均为 I/2；但 Bell 态的 σₓ⊗σₓ 期望为 1，等权 00/11 混合态为 0，联合相干信息不同。"),
]
FLASHCARDS = [{"id": f"card-{i}", "chapter": c, "front": f, "back": b} for i, (c, f, b) in enumerate(_CARD_DATA, 1)]
