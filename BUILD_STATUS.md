# 构建状态（接手者只看这一份即可续做）

状态：**IN_PROGRESS**（未标 COMPLETE 前不许对外宣称建成）

## 定位与硬标准（不可降标）
- 中文 Streamlit 线代复习平台，面向物理本科生；单一流程「完整推导 → 手写代码检验 → 自动检查」；每章含：动机与直觉 / 定义 / 核心推导 / 适用条件与失效情形 / 几何直觉 / 与量子力学的连接 六小节 + 可视化（能画的章必配）+ 代码练习（AST 防作弊检查）+ 书面习题 6 道（基础2/进阶2/考试风格2，真分步解答+常见错误点拨）+ 章末自测 5 题（页面内一次一题）。
- 推导标准：完整推导链、逐步展开、禁「易得/显然/不难验证」跳步；推导前有动机、后有适用条件与失效情形。
- 差异化依据见 docs/prior-art.md（Muse 检索，2026-10-25）。
- Codex 调用规矩：`~/.local/bin/codex exec -m gpt-6.1-sol -c model_reasoning_effort="high" --skip-git-repo-check -C <repo> "<brief>" </dev/null`；每轮任务书在 /tmp/la_briefN.txt（若已丢失，按本文件「未完成项」重写即可）；每轮收尾要求 Codex 跑 review-agent 并在报告中交代 skills 使用。
- 质量令（用户 00:52/00:53）：全部建成后先做①独立审读轮（review-agent + 逐章抽查推导完整性/习题解答/防作弊有效性，修完问题）再做②充实检查轮（推导密度、算例、QM 连接是否具体），两轮记录写入本文件后才许标 COMPLETE。撞 Codex 限额（usage limit）时在本文件顶部写 `QUOTA_EXHAUSTED <北京时间>` 与续做起点，commit+push 后结束，勿降标赶工。

## 已完成
- [x] GitHub 先行检索 → docs/prior-art.md（Muse）
- [x] 私有仓库 qiaoyuzheng0804-create/linear-algebra-review-lab 建立、克隆、Codex trust 已加
- [x] 第 1 轮（Codex，tokens 103,276）：基础设施全套——app.py（18 个导航页与通用章节渲染器）、linalg_lab/{progress,checks,viz,extras}.py、地基模块（4 节+自检）、第 1–3 章完整内容（每章 theory 约 3200–3900 字符、书面 6、测验 5、代码练习 1–2）、公式手册数据（42 条）、入学诊断 14 题、闪卡 36 张、知识地图依赖数据、错题本/每日一题页面
  - 验收（Muse）：pytest 独立复跑 55 passed；grep 全仓无「显然/易得/不难验证/类似可证」；通读第 1 章核心推导（Steinitz 交换引理完整证明、双向等价、条件与失效情形齐全）判定合格；review-agent 经 Codex 调用修复 4 项后 No findings
  - 第 4–10 章目前为「建设中」占位页

## 未完成（按顺序续做）
1. [ ] 第 2 轮：第 4–7 章完整内容（规格同第 1–3 章；核心推导指定见各轮任务书要点——第4章行列式：按置换定义→交错多线性→det(AB)=detA·detB 推导链、几何面积意义+viz det_area；第5章特征值：特征多项式、代数/几何重数、Cayley–Hamilton 证明要点完整版、viz eigen_demo；第6章内积与 Gram-Schmidt：Cauchy–Schwarz 由判别式完整推导、GS 正交化构造与 QR、viz gram_schmidt 分步；第7章对角化与矩阵指数：对角化判据推导（特征向量张满）、e^A 定义与微分方程组解、Jordan 形失效情形、viz matrix_power）。每章书面 6（2/2/2 分层）、测验 5、代码练习 ≥1（从零实现、banned 清单齐全）。
2. [ ] 第 3 轮：第 8–10 章（第8章厄米/酉/谱定理为重点章：谱定理对厄米矩阵的完整归纳证明、酉对角化、QM 连接写到可观测量/测量概率层面要具体；第9章 SVD：由 A^TA 谱分解构造 SVD 的完整推导、与特征分解的差别来自哪一步、最小二乘正规方程推导、viz svd_ellipse；第10章张量积：由双线性/多重线性构造、纠缠态不可分解的证明、viz tensor_dims、部分迹初步）+ 公式手册第 4–10 章校对补全 + 闪卡扩充至覆盖全部章节（≥70 张）+ 诊断题与新章节映射复核。
3. [ ] 第 4 轮（打磨轮①独立审读）：review-agent 全仓审查 + 验收方逐章各抽 1 个核心推导通读、各章抽 1 道书面题解答验算、亲手构造作弊代码（别名/间接调用等）验证 checks 拦截；问题修完并复跑全量 pytest。
4. [ ] 第 5 轮（打磨轮②充实检查）：逐章核推导密度（核心结论无只有结论的条目）、每章至少 1 个手算算例、QM 连接具体化（点名具体物理对象/公式，不许口号）；缺失补齐后复跑 pytest。
5. [ ] 最终验收（验收方亲自）：本地 `streamlit run app.py` 实际打开逐页确认渲染与滑杆交互；错题本加入/重练/清除实测；亲手验算 ≥3 道代码练习参考解数值、≥5 道书面题解答；全绿后本文件标 COMPLETE 并 commit+push。
6. [ ] 部署尝试：按 neural-network-learning-lab 路径试部署 Streamlit Cloud；需用户登录授权则如实记录卡点（本地+GitHub 交付为准），不要冒充上线。

## 各轮 Codex tokens
- 第 1 轮：103,276
