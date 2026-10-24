# GitHub 先行检索结论（2026-10-25，Muse 检索）

检索方式：`gh search repos` 多组关键词 + 网页检索。代表性现成项目：

| 项目 | 星标量级 | 覆盖了什么 | 没覆盖什么 |
|---|---|---|---|
| fastai/numerical-linear-algebra | ~11k | Jupyter 笔记本形式的计算线性代数课程（英文），偏数据科学 | 无中文、无自动检查的手写练习、无量子力学衔接、无 Web 应用与进度 |
| kenjihiranabe/The-Art-of-Linear-Algebra | ~21.8k | Strang 教材的图形化笔记（静态图示） | 只有图解笔记，无练习、无检查、无交互应用 |
| QBobWatson/ila | ~0.8k | Interactive Linear Algebra 英文在线教材 | 英文教材形态，无代码手写练习与判题、无物理/量子定向 |
| weijie-chen/Linear-Algebra-With-Python | ~2.6k | 英文讲义 + Python 计算与可视化 | 讲义形态，无判题系统、无进度、无量子衔接 |
| liuyubobobo/Play-with-Linear-Algebra | ~0.4k | 中文 MOOC《玩转线性代数》示例代码 | 课程配套代码，无推导讲解页、无自动检查、无错题/进度系统 |
| zlotus/notes-linear-algebra 等笔记仓 | ~4.3k | 中文线性代数笔记 | 静态笔记，无交互、无练习判题 |
| Ximera《Linear Algebra: An Interactive Introduction》 | — | 英文交互教材（GeoGebra/Octave 嵌入、机器判题） | 英文、通用数学课定位，无手写代码练习、无物理/量子衔接 |

## 结论与本平台差异点

现成项目没有一个同时具备以下组合，故不属原样重做：

1. **中文 Streamlit Web 应用**（与用户已在用的 neural-network-learning-lab 同技术栈、同使用习惯），不是笔记本/静态笔记；
2. **面向物理本科生**：每章设“几何直觉”与“与量子力学的连接”小节，第 8 章（厄米矩阵、酉矩阵与谱定理）与第 10 章（张量积）为量子衔接重点章；
3. **单一学习流程**：先学完整原理与推导 → 再手写代码检验（numpy 从零实现、自动检查且能识别直接调库函数的作弊）→ 看检查结果；不并列多种学法；
4. **完整学习闭环**：地基模块、公式手册、入学诊断测、分层书面习题（基础/进阶/考试风格，含分步解答与常见错误点拨）、错题本、概念闪卡、每日一题、知识地图、本地进度持久化。

（作者：检索与本文 — Muse）
