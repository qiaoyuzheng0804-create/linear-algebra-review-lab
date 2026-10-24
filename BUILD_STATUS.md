# 构建状态（接手者只看这一份即可续做）

**QUOTA_EXHAUSTED 2026-10-25 01:44（北京时间）**——第 4 轮（独立审读）执行约 75 秒后 Codex 返回 `You've hit your usage limit ... try again at 9:23 PM`（按 UTC 读即北京时间 10-25 05:23 恢复，与 5 小时窗口吻合）。撞限额时工作树干净、第 4 轮未产生任何文件改动（它当时还在审读清单第 4 项的作弊探针阶段）。

**续做起点（额度恢复后从这里接，按序）**：
1. 重跑第 4 轮：任务书全文在 `~/workspace/la-briefs/la_brief4.txt`（持久目录，不在 /tmp）。命令：`~/.local/bin/codex exec -m gpt-6.1-sol -c model_reasoning_effort="high" --skip-git-repo-check -C ~/workspace/repos/linear-algebra-review-lab "$(cat ~/workspace/la-briefs/la_brief4.txt)" </dev/null`。验收方已先行完成其中两项，可直接采信并在 BUILD_STATUS 注明：(a) 防作弊拦截实测——direct/别名导入/getattr/eval/星号导入五类变体对第 3 章练习全部被拦截；(b) 全部 11 道代码练习的 reference 经真实文件探针（~/workspace/la-briefs/probe_checks.py，注意需要 `if __name__ == "__main__"` 守卫否则 spawn 递归）逐一 check_code 通过。第 4 轮还须补：书面题/quiz 答案复核、推导完整性全仓扫描、公式手册一致性、AppTest 全页遍历、review-agent 全仓审查。第 4 轮审读中曾探到 `np.load(..., allow_pickle=True)` 未被 banned 命中（探针输出为空命中清单），续做时评估是否构成绕过面并处理。
2. 第 5 轮（充实检查）：任务书待写，要点已记在本文件未完成项第 4 条；另加明确目标——第 1~3 章 theory 现为 3,885/3,222/3,806 字符，明显薄于第 8~10 章的 8,838/7,649/8,955，充实轮优先把第 1~3 章补到同等密度（增推导与手算算例，不灌水）。
3. 最终验收（本文件未完成项第 5 条）与部署尝试（第 6 条；Streamlit Cloud 需用户 GitHub OAuth 登录，大概率卡住，如实记录即可）。
4. 两轮打磨+最终验收全过后，把本文件状态改为 COMPLETE、删去本 QUOTA 段落或注明恢复时间，commit+push。

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
  - 第 1 轮结束时，第 4–10 章为「建设中」占位页；第 2 轮后仅第 8–10 章仍待建设；第 3 轮已补齐十章，仍待全仓两轮打磨与最终验收


- [x] 第 2 轮（Codex 实现，独立 review-agent 审查，2026-10-24 UTC）：
  - 新增 `linalg_lab/content/ch04.py`～`ch07.py`，注册第 4～7 章；复用通用渲染器与 `det_area`、`eigen_demo`、`gram_schmidt`、`matrix_power`，无单章页面特判。
  - 每章六小节、完整推导与手算算例、具体 QM 公式；每章书面题 6 道（基础/进阶/考试风格各 2）、测验 5 题、代码练习 1 道，共书面 24、测验 20、代码 4。
  - 理论字符数（Python `len(theory_md)`，含 LaTeX/空白）：第4章 5,904；第5章 5,618；第6章 6,282；第7章 6,270。
  - 对齐导航/公式手册/诊断/知识依赖：第7章为对角化与矩阵指数，第9章为 SVD；已开放章不再被诊断建议标成建设中。
  - 新增 `tests/test_round2.py`；适配第1轮注册数量与已开放章断言。参考解经现有检查器执行，还独立验证行列式、特征残差、QR 恒等式、指数级数及近根边界；四章真实页面的可视化参数和测验交互均由 AppTest 验证。
  - 初审修复：NumPy 内部重导出/星号导入可绕过禁用 API；`numpy.matrix.I` 隐藏调用求逆及矩阵幂；二阶 `trace²−4det` 消去使合法近邻特征值误报重根。保留原检查器框架，补 API 导入名、星号导入与封装禁令；判别式改用 `(a-d)²+4bc`，增加真实绕过与近根回归。
  - review-agent 最终结论 **No findings**，独立针对回归 **27 passed**；主代理最终全量 `MPLCONFIGDIR=/tmp/linalg-round2-mpl .venv/bin/python -m pytest -q`：**182 passed in 48.04s**。`git diff --check` 通过。
  - skills 盘点：系统目录含 imagegen、openai-docs、review-agent、skill-creator、skill-installer。实际应用 `review-agent`，按其只读审查流程由独立代理复核；其余四个与本轮内容/算法实现无关，未调用。另按托管环境要求读取 `cloud-environment-runtime`，本会话未提供 environment_status 工具，本轮无需联网或容器操作。
  - 本轮未 git commit/push；整体仍为 **IN_PROGRESS**。后续全仓独立审读轮、充实检查轮与最终人工验收仍须执行，不能将本轮审查当成全项目完成验收。

- [x] 第 3 轮（本轮实现，2026-10-24 UTC）：
  - 新增 `linalg_lab/content/ch08.py`～`ch10.py` 并全部注册，十章内容模块已开放。每章六小节、书面 6 道（基础/进阶/考试风格各 2）、测验 5 题、代码练习 1 道，含手算数字算例和具体 QM 对象与公式。
  - 第8章：厄米实谱与异值正交、代数学基本定理作为首个本征向量来源、正交补不变及限制厄米的完整归纳谱定理、极化恒等式及酉性等价链、Schur 三角化及正规三角阵对角化；代码手写简并谱投影重构与候选谱残差；新增实对称谱图。
  - 第9章：A†A 半正定谱分解构造任意矩形/亏秩 SVD、零值正交补、极小极大刻画与 Eckart–Young 两种范数下界、残差正交推正规方程、伪逆全部解与最小范数；代码多起点幂迭代+收缩构造正奇异三元组并返回奇异值、重构、伪逆；复用已有三步 SVD 图并增加追踪向量。
  - 第10章：自由基系数表构造及维数、双线性与泛性质、Kronecker 分块及复合、Bell 态系数矛盾、完整复数密度矩阵与部分迹、Schmidt/SVD 共轭对应与联合纯态纠缠判据；代码手写 Kronecker 与双向部分迹；张量维数图接入两个交互滑杆。
  - 公式手册第4～10章校对补充，含谱投影测量概率、极化恒等式、SVD/极小极大/伪逆与正规方程、Schmidt 与部分迹公式；闪卡保留旧 ID 并追加至 **81 张**，覆盖全部十章。诊断保持 14 题及旧 ID，复核新章映射、未答题处理和推荐；知识地图第10章增加第9章前置，与 Schmidt 推导一致。
  - 理论字符数（Python `len(theory_md)`，含 LaTeX/空白）：第8章 **8,838**；第9章 **7,649**；第10章 **8,955**。第8章超过第6章 6,282 与第7章 6,270。
  - 新增 `tests/test_round3.py`（99 项）；参考解经过现有 checks 独立进程检查，并用库函数/恒等式独立核对：简并测量概率、候选谱四项残差、SVD 与四条 Penrose 恒等式、局域期望和部分迹、复数 Schmidt 约定。覆盖作弊别名/解包/内部重导出/动态属性/替代高层 API、图形返回与真实页面参数/测验持久化。第一轮测试的注册范围和已开放章推荐断言同步更新。
  - 审查修复后最终全量验证：`MPLCONFIGDIR=/tmp/linalg-round3-mpl .venv/bin/python -m pytest -q`：**293 passed in 37.64s**（本轮新增 `test_round3.py` 99 项）；`git diff --check` 通过。
  - 独立 review-agent 最终结论 **No findings**。初审发现新增 Cayley–Hamilton 闪卡所述证明路线与第5章正文不一致，已改为正文实际采用的可对角化先证→复三角化与对角微扰→多项式连续性取极限，复核通过。审查代理独立运行第3轮非页面测试 **96 passed, 3 deselected**，并以 **3,000** 个合法范围随机矩阵压测 SVD（重复谱、0.8间隔边界、秩亏、0.4～3尺度），奇异值/重构/伪逆最大绝对误差分别约 1.78e−15/2.18e−14/1.77e−12。README 同步十章开放与81张闪卡状态。
  - skills 逐项：系统目录共 imagegen、openai-docs、review-agent、skill-creator、skill-installer。实际应用 `review-agent`，由独立代理按其 SKILL.md 复核全部本轮改动；其余四项与本轮数学内容、NumPy 算法及 Matplotlib 交互无关，未调用。另读取 `cloud-environment-runtime`；未提供 environment_status 工具，本轮未执行网络/容器操作。
  - 未 git commit/push。整体仍为 **IN_PROGRESS**；本轮审查不替代后续全仓独立审读轮、充实检查轮与人工最终验收。
  - 本轮 tokens：本会话执行尚未结束，无法在结束前读取外层 `codex exec` 的最终 `tokens used`；由调用方在结束后补记，不扫日志估算。

## 未完成（按顺序续做）
1. [x] 第 2 轮：第 4–7 章完整内容（规格同第 1–3 章；核心推导指定见各轮任务书要点——第4章行列式：按置换定义→交错多线性→det(AB)=detA·detB 推导链、几何面积意义+viz det_area；第5章特征值：特征多项式、代数/几何重数、Cayley–Hamilton 证明要点完整版、viz eigen_demo；第6章内积与 Gram-Schmidt：Cauchy–Schwarz 由判别式完整推导、GS 正交化构造与 QR、viz gram_schmidt 分步；第7章对角化与矩阵指数：对角化判据推导（特征向量张满）、e^A 定义与微分方程组解、Jordan 形失效情形、viz matrix_power）。每章书面 6（2/2/2 分层）、测验 5、代码练习 ≥1（从零实现、banned 清单齐全）。
2. [x] 第 3 轮：第 8–10 章（第8章厄米/酉/谱定理为重点章：谱定理对厄米矩阵的完整归纳证明、酉对角化、QM 连接写到可观测量/测量概率层面要具体；第9章 SVD：由 A^TA 谱分解构造 SVD 的完整推导、与特征分解的差别来自哪一步、最小二乘正规方程推导、viz svd_ellipse；第10章张量积：由双线性/多重线性构造、纠缠态不可分解的证明、viz tensor_dims、部分迹初步）+ 公式手册第 4–10 章校对补全 + 闪卡扩充至覆盖全部章节（≥70 张）+ 诊断题与新章节映射复核。
3. [ ] 第 4 轮（打磨轮①独立审读）：review-agent 全仓审查 + 验收方逐章各抽 1 个核心推导通读、各章抽 1 道书面题解答验算、亲手构造作弊代码（别名/间接调用等）验证 checks 拦截；问题修完并复跑全量 pytest。
4. [ ] 第 5 轮（打磨轮②充实检查）：逐章核推导密度（核心结论无只有结论的条目）、每章至少 1 个手算算例、QM 连接具体化（点名具体物理对象/公式，不许口号）；缺失补齐后复跑 pytest。
5. [ ] 最终验收（验收方亲自）：本地 `streamlit run app.py` 实际打开逐页确认渲染与滑杆交互；错题本加入/重练/清除实测；亲手验算 ≥3 道代码练习参考解数值、≥5 道书面题解答；全绿后本文件标 COMPLETE 并 commit+push。
6. [ ] 部署尝试：按 neural-network-learning-lab 路径试部署 Streamlit Cloud；需用户登录授权则如实记录卡点（本地+GitHub 交付为准），不要冒充上线。

## 各轮 Codex tokens
- 第 1 轮：103,276

- 第 2 轮：无法逐章拆分。收尾日志快照（主代理截至 2026-10-24 17:15:35 UTC，review-agent 截至 17:14:52 UTC）累计 **5,665,275 tokens**，含缓存输入 5,425,792；非缓存输入＋输出 **239,483 tokens**，已包含 review-agent。此为日志采样时点值，不包含采样后的收尾消息。
  - 验收方（Muse）更正口径：本轮 `codex exec` 结束行 "tokens used" = **138,570**（与第 1 轮 103,276 同口径，即 total − cached input）。上面 5,665,275 的日志快照混入了同夜并行运行的其他 Codex 会话，不作本项目用量。
  - 验收（Muse）：全量 pytest 独立复跑 **182 passed**；grep 第 4–7 章无跳步词；通读第 6 章 Cauchy–Schwarz 推导（判别式法 + 复情形相位旋转处理 + 等号条件）判定正确完整。
- 第 3 轮：`codex exec` 结束行 "tokens used" = **133,624**（同口径）。
  - 验收（Muse）：全量 pytest 独立复跑 **293 passed**；grep 第 8–10 章无跳步词；通读第 8 章推导 8.1/8.2 与第 10 章 Schmidt/部分迹/纯度判据段（含判据对联合混合态失效的断点指认）判定合格；亲手 numpy 复算三处数字全对——√3/2,1/2 例 ρ_A=diag(3/4,1/4)、纯度 5/8；Bell 态 ρ_A=I/2、纯度 1/2；σ_x 对 |+z⟩ 测 +1 概率 1/2。
