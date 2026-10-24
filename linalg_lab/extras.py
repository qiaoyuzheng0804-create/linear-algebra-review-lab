"""公式、诊断、知识依赖与闪卡数据；无 Streamlit 依赖。"""
from collections import Counter

from .content import CHAPTER_TITLES

# 每项为（LaTeX，成立条件）；第 4～10 章为后续内容提供标准公式索引。
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
        (r"\det A=\sum_{\sigma\in S_n}\operatorname{sgn}(\sigma)\prod_{i=1}^n a_{i,\sigma(i)}", "A 为 n 阶方阵；Leibniz 定义"),
        (r"\det(AB)=\det A\det B", "A、B 为同阶方阵"),
        (r"A\text{可逆}\iff\det A\ne0", "方阵；一般矩形矩阵无此行列式判据"),
        (r"\mathrm{Vol}(A\Omega)=|\det A|\,\mathrm{Vol}(\Omega)", "实方阵作用于可测区域，使用普通体积"),
    ],
    "第5章": [
        (r"Av=\lambda v\ (v\ne0)\iff (A-\lambda I)v=0", "有限维方阵；λ 在所选数域中"),
        (r"\det(A-\lambda I)=0", "λ 为特征值的充要条件，A 为方阵"),
        (r"A=P\Lambda P^{-1}\iff A\text{有一组特征向量基}", "在所选数域上；重根不自动保证可对角化"),
        (r"A^k=P\Lambda^kP^{-1}", "A 可对角化；k 为非负整数，负整数还要求 A 可逆"),
    ],
    "第6章": [
        (r"\langle u,v\rangle=u^\dagger v", "复欧氏空间，第一变量共轭线性"),
        (r"u_k=v_k-\sum_{j<k}\langle q_j,v_k\rangle q_j,\quad q_k=u_k/\|u_k\|", "q_j 正交归一，输入组无关确保 u_k 非零"),
        (r"A^\dagger A\hat x=A^\dagger b", "最小二乘的正规方程；有限维实或复矩阵"),
        (r"\hat x=(A^\dagger A)^{-1}A^\dagger b", "A 满列秩，最小二乘解唯一；数值实现宜用 QR 或 SVD"),
    ],
    "第7章": [
        (r"A=U\Sigma V^\dagger", "任意有限实或复 m×n 矩阵；U、V 酉，奇异值非负"),
        (r"A^+=V\Sigma^+U^\dagger", "任意有限矩阵；Σ⁺ 将非零奇异值取倒数，零保留"),
        (r"\|A-A_k\|_2=\sigma_{k+1},\quad \|A-A_k\|_F^2=\sum_{j>k}\sigma_j^2", "奇异值降序，A_k 为截断 SVD，0≤k<min(m,n)；最优秩至多 k 近似"),
    ],
    "第8章": [
        (r"H=H^\dagger\Longrightarrow H=U\Lambda U^\dagger,\quad\Lambda\text{为实对角阵}", "有限维厄米矩阵，U 酉；谱定理"),
        (r"A^\dagger A=AA^\dagger\iff A=U\Lambda U^\dagger", "复方阵；正规矩阵的酉对角化判据"),
        (r"U^\dagger U=I\Longrightarrow\|U\psi\|=\|\psi\|", "有限维酉方阵，标准复内积"),
        (r"H=\sum_\lambda\lambda P_\lambda,\quad P_\lambda P_\mu=\delta_{\lambda\mu}P_\lambda", "有限维厄米矩阵，P_λ 为整个本征空间的正交投影"),
    ],
    "第9章": [
        (r"e^{tA}=\sum_{k=0}^{\infty}\frac{t^kA^k}{k!}", "任意有限方阵，有限 t；级数绝对收敛"),
        (r"\dot x=Ax,\quad x(t)=e^{tA}x(0)", "常系数有限维系统，初值设在 t=0"),
        (r"e^{A+B}=e^Ae^B", "充分条件 AB=BA；一般不对易时不能套用"),
        (r"U(t)=e^{-iHt/\hbar},\quad U(t)^\dagger U(t)=I", "有限维、时间无关厄米 H，实 t、非零实 ħ"),
    ],
    "第10章": [
        (r"\dim(V\otimes W)=\dim V\,\dim W", "有限维同数域向量空间"),
        (r"(A\otimes B)(C\otimes D)=AC\otimes BD", "AC 与 BD 的维数分别相容"),
        (r"(A\otimes B)^\dagger=A^\dagger\otimes B^\dagger", "有限实或复矩阵"),
        (r"\operatorname{Tr}_B\rho=\sum_j(I\otimes\langle j|)\rho(I\otimes|j\rangle)", "有限维复合系统，|j⟩ 是 B 的任一正交归一基"),
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
    "ch7": ["foundation", "ch5", "ch6"],
    "ch8": ["foundation", "ch5", "ch6"],
    "ch9": ["foundation", "ch5", "ch8"],
    "ch10": ["foundation", "ch2", "ch8"],
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
    {"chapter": "ch5", "q_md": "n 阶矩阵可对角化的充要条件是？", "options": ["有 n 个不同条目", "行列式非零", "至少有一个特征值", "有 n 个线性无关的特征向量"], "answer": 3, "explain_md": "将特征向量排成可逆 P，逐列 AP=PΛ，才能得到相似对角化。"},
    {"chapter": "ch6", "q_md": "Gram–Schmidt 中某步残差为零说明？", "options": ["当前向量在前面向量的张成内", "已得新单位向量", "矩阵必为方阵", "内积失效"], "answer": 0, "explain_md": "扣除已有方向投影后没有剩余分量，因此归一化分母为零，应识别相关性。"},
    {"chapter": "ch7", "q_md": "SVD 对哪些有限矩阵存在？", "options": ["仅可逆方阵", "任意实或复矩阵", "仅厄米矩阵", "仅对角矩阵"], "answer": 1, "explain_md": "SVD 不要求方阵、满秩或特征向量张满；左右可使用不同的正交或酉基。"},
    {"chapter": "ch8", "q_md": "厄米矩阵的特征值必为？", "options": ["正数", "整数", "实数", "纯虚数"], "answer": 2, "explain_md": r"若 $Hv=\lambda v$，则 $\lambda=v^\dagger Hv/(v^\dagger v)$，分子因 H 厄米等于自身共轭，故 λ 实。"},
    {"chapter": "ch9", "q_md": r"常系数系统 $\dot x=Ax$ 的传播矩阵为？", "options": [r"$tA$", r"$A^{-1}$", r"$A+tI$", r"$e^{tA}$"], "answer": 3, "explain_md": r"指数级数逐项求导给 $d e^{tA}/dt=Ae^{tA}$，且在 t=0 等于 I。"},
    {"chapter": "ch10", "q_md": "二维空间与三维空间的张量积维数为？", "options": ["5", "6", "3", "2"], "answer": 1, "explain_md": "两套基两两配对得到 2×3=6 个张量积基向量。"},
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
    if start not in {"foundation", "ch1", "ch2", "ch3"}:
        result += "该章建设中，可先查公式手册并复核第1～3章。"
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
]
FLASHCARDS = [{"id": f"card-{i}", "chapter": c, "front": f, "back": b} for i, (c, f, b) in enumerate(_CARD_DATA, 1)]
