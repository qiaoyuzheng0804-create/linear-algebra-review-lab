"""Streamlit 入口；所有章节共用数据驱动的渲染器。"""
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from linalg_lab import progress, viz
from linalg_lab.checks import check_code
from linalg_lab.content import CHAPTERS, CHAPTER_TITLES
from linalg_lab.content.foundation import FOUNDATION
from linalg_lab.extras import (
    DIAGNOSTIC, FLASHCARDS, FORMULA_SHEET, KNOWLEDGE_MAP, chapter_label, recommend,
)


def persist(data):
    try:
        progress.save(data)
    except OSError as exc:
        st.error(f"进度未能写入磁盘：{exc}。当前页面中的记录仍在内存中。")


def show_figure(figure):
    try:
        st.pyplot(figure)
    finally:
        plt.close(figure)


def mistake_control(data, item_id, key):
    def changed():
        progress.set_mistake(data, item_id, st.session_state[key])
        persist(data)

    # 自动收录和手动勾选共享数据源；回调只响应用户操作。
    st.session_state[key] = item_id in data["mistakes"]
    st.checkbox("加入错题本", key=key, on_change=changed)


def render_written(item, chapter, data, prefix):
    st.markdown(item["problem_md"])
    with st.expander("查看分步解答与常见错误"):
        st.markdown(item["solution_md"])
        st.markdown("**常见错误点拨**")
        st.markdown(item["pitfalls_md"])
    state = progress.chapter_progress(data, chapter["id"])
    done = item["id"] in state.get("written", [])
    checked = st.checkbox("已独立完成并核对", value=done, key=f"{prefix}-done-{item['id']}")
    if checked != done:
        progress.mark_written(data, chapter["id"], item["id"], checked)
        persist(data)
    mistake_control(data, item["id"], f"{prefix}-mistake-{item['id']}")


def render_quiz(question, chapter, data, prefix):
    qid = question["id"]
    with st.form(f"{prefix}-form-{qid}"):
        st.markdown(question["q_md"])
        selected = st.radio("选择答案", range(4), format_func=lambda i: question["options"][i], index=None, key=f"{prefix}-choice-{qid}")
        submitted = st.form_submit_button("提交答案")
    if submitted:
        if selected is None:
            st.info("请先选择一个答案。")
        else:
            correct = selected == question["answer"]
            progress.mark_quiz(data, chapter["id"], qid, correct)
            if not correct:
                progress.set_mistake(data, qid)
            persist(data)
            st.session_state[f"{prefix}-result-{qid}"] = correct
    if f"{prefix}-result-{qid}" in st.session_state:
        correct = st.session_state[f"{prefix}-result-{qid}"]
        (st.success if correct else st.error)("回答正确。" if correct else "回答错误，已加入错题本。")
        st.markdown("**正确答案：** " + question["options"][question["answer"]])
        st.markdown(question["explain_md"])
    mistake_control(data, qid, f"{prefix}-mistake-{qid}")


def render_visualization(key, prefix):
    st.markdown("**交互观察**")
    if key in {"grid_transform", "det_area", "eigen_demo", "svd_ellipse", "matrix_power"}:
        st.caption("调整二维矩阵的四个元素，观察同一规则如何作用于所有向量。")
        columns = st.columns(4)
        values = [column.number_input(label, -3., 3., value, .1, key=f"{prefix}-{key}-{label}")
                  for column, label, value in zip(columns, ("a₁₁", "a₁₂", "a₂₁", "a₂₂"), (1., .6, 0., 1.))]
        matrix = np.array(values).reshape(2, 2)
        if key == "matrix_power":
            k = st.slider("迭代次数", 0, 25, 5, key=f"{prefix}-power")
            figure = viz.matrix_power(matrix, k)
        else:
            figure = getattr(viz, key)(matrix)
    elif key == "gram_schmidt":
        step = st.slider("正交化步骤", 1, 2, 1, key=f"{prefix}-gs")
        figure = viz.gram_schmidt([[1, 1], [2, 0]], step)
    elif key == "hermitian_spectrum":
        a = st.slider("厄米矩阵对角元 a", -3., 3., 1., .1, key=f"{prefix}-hermitian-a")
        d = st.slider("厄米矩阵对角元 d", -3., 3., -1., .1, key=f"{prefix}-hermitian-d")
        b = st.slider("实对称耦合 b", -3., 3., .6, .1, key=f"{prefix}-hermitian-b")
        figure = viz.hermitian_spectrum([[a, b], [b, d]])
    elif key == "tensor_dims":
        m = st.slider("A 系统维数", 1, 6, 2, key=f"{prefix}-tensor-a")
        n = st.slider("B 系统维数", 1, 6, 3, key=f"{prefix}-tensor-b")
        figure = viz.tensor_dims(m, n)
    else:
        st.error(f"尚未注册可视化：{key}")
        return
    show_figure(figure)


def render_chapter(chapter, data):
    """仅使用 CHAPTER 协议；不包含单章特判。"""
    cid = chapter["id"]
    st.title(f"第{chapter['num']}章　{chapter['title']}")
    st.caption("前置：" + "、".join(chapter_label(x) for x in chapter["depends_on"]))
    completion_slot = st.empty()
    theory_tab, code_tab, result_tab, written_tab, quiz_tab = st.tabs(["① 理论与推导", "② 手写代码", "③ 检查结果", "分层书面题", "概念测验"])
    with theory_tab:
        st.markdown(chapter["theory_md"])
        for key in chapter["viz"]:
            render_visualization(key, cid)
        done = progress.chapter_progress(data, cid).get("theory_done", False)
        checked = st.checkbox("已理解推导并能复述依据", value=done, key=f"theory-{cid}")
        if checked != done:
            progress.mark_theory(data, cid, checked)
            persist(data)
    with code_tab:
        for exercise in chapter["code_exercises"]:
            eid = exercise["id"]
            st.subheader(exercise["title"])
            st.markdown(exercise["prompt_md"])
            st.caption("已提供 np。输入函数代码，随机检查将比较数值、形状和返回结构。")
            drafts = st.session_state.setdefault("code_drafts", {})
            widget_key = f"source-{eid}"
            if widget_key not in st.session_state:
                st.session_state[widget_key] = drafts.get(eid, exercise["starter"])

            def remember_draft(exercise_id=eid, editor_key=widget_key):
                st.session_state["code_drafts"][exercise_id] = st.session_state[editor_key]

            source = st.text_area("代码", height=320, key=widget_key, on_change=remember_draft)
            if st.button("运行检查", key=f"run-{eid}", type="primary"):
                with st.spinner("正在检查代码与随机用例……"):
                    passed, feedback = check_code(source, exercise)
                st.session_state[f"feedback-{eid}"] = (passed, feedback)
                if passed:
                    progress.mark_code(data, cid, eid)
                    st.success("检查通过，完成记录已保存。可在「检查结果」查看逐条反馈。")
                else:
                    progress.set_mistake(data, eid)
                    st.error("检查未通过，已加入错题本。请在「检查结果」查看原因。")
                persist(data)
    with result_tab:
        for exercise in chapter["code_exercises"]:
            eid = exercise["id"]
            st.subheader(exercise["title"])
            if f"feedback-{eid}" in st.session_state:
                passed, feedback = st.session_state[f"feedback-{eid}"]
                (st.success if passed else st.error)("本次通过" if passed else "本次未通过")
                for line in feedback:
                    st.write(line)
            elif eid in progress.chapter_progress(data, cid).get("code_ex", []):
                st.success("此练习已有通过记录。再次运行可生成本次详细反馈。")
            else:
                st.info("尚无检查结果。先完成手写代码，再运行检查。")
    with written_tab:
        for item in chapter["written"]:
            st.subheader(f"{item['level']} · {item['id']}")
            render_written(item, chapter, data, "chapter")
    with quiz_tab:
        questions = chapter["quiz"]
        question_id = st.selectbox("选择测验题", [q["id"] for q in questions],
                                   format_func=lambda qid: f"第 {next(i for i, q in enumerate(questions, 1) if q['id'] == qid)} 题",
                                   key=f"quiz-select-{cid}")
        question = next(q for q in questions if q["id"] == question_id)
        st.subheader(f"第 {questions.index(question) + 1} 题")
        render_quiz(question, chapter, data, "chapter")
    # 填充顶部占位，在本轮交互更新数据后计算，避免进度滞后一轮。
    completion_slot.progress(progress.completion(data, chapter), text="本章完成度：理论、代码、测验与书面题按项目计数")


def render_overview(data):
    st.title("线性代数复习实验室")
    st.markdown("面向大三物理本科生：从完整推导出发，手写算法检验理解，再连接几何与量子力学。")
    metrics = st.columns(3)
    metrics[0].metric("本轮开放章节", f"{len(CHAPTERS)} / {len(CHAPTER_TITLES)}")
    metrics[1].metric("已掌握闪卡", f"{len(set(data['flashcards_mastered']) & {c['id'] for c in FLASHCARDS})} / {len(FLASHCARDS)}")
    metrics[2].metric("错题本", len(data["mistakes"]))
    st.markdown("**学习流程**：入学诊断定位起点 → 地基补齐语言 → 逐章阅读并复述推导 → 手写代码 → 查看检查结果 → 完成书面题与测验 → 用错题本、闪卡和每日一题复习。")
    for number, title in CHAPTER_TITLES.items():
        st.subheader(f"第{number}章　{title}")
        chapter = CHAPTERS.get(number)
        if chapter is None:
            st.progress(0., text="建设中")
            continue
        state = progress.chapter_progress(data, chapter["id"])
        st.progress(progress.completion(data, chapter))
        st.caption(f"理论：{'完成' if state.get('theory_done') else '待学'} · 代码：{len(state.get('code_ex', []))}/{len(chapter['code_exercises'])} · 测验答对：{state.get('quiz', 0)}/{len(chapter['quiz'])} · 书面题：{len(state.get('written', []))}/{len(chapter['written'])}")


def render_foundation(data):
    st.title("地基")
    st.markdown("先把复振幅、几何、坐标和映射的语言对齐，再进入抽象推导。")
    for section in FOUNDATION:
        st.header(section["title"])
        st.markdown(section["md"])
        for i, check in enumerate(section["selfcheck"], 1):
            st.markdown(f"**自检 {i}**　" + check["q"])
            with st.expander("核对答案"):
                st.markdown(check["answer_md"])
    done = progress.chapter_progress(data, "foundation").get("theory_done", False)
    checked = st.checkbox("四节地基已完成", value=done, key="foundation-done")
    if checked != done:
        progress.mark_theory(data, "foundation", checked)
        persist(data)


def render_formula_sheet():
    st.title("公式手册")
    st.caption("每条公式与成立条件、对应章节的推导配对，可按章节或条件检索。")
    query = st.text_input("筛选章节或条件关键词")
    for label, formulas in FORMULA_SHEET.items():
        selected = [(formula, condition) for formula, condition in formulas if not query or query.lower() in f"{label} {formula} {condition}".lower()]
        if selected:
            st.header(label)
            for formula, condition in selected:
                st.markdown("$$" + formula + "$$")
                st.markdown("成立条件：" + condition)


def render_diagnostic(data):
    st.title("入学诊断")
    st.markdown("14 题覆盖地基与全课程。先凭现有理解作答，再根据错题分布选择起点。")
    selected = {}
    with st.form("diagnostic-form"):
        for i, question in enumerate(DIAGNOSTIC, 1):
            st.markdown(f"**{i} · 考查：{chapter_label(question['chapter'])}**")
            st.markdown(question["q_md"])
            selected[question["id"]] = st.radio("答案", range(4), format_func=lambda j, q=question: q["options"][j], index=None, key=f"diagnostic-{question['id']}")
        submitted = st.form_submit_button("提交诊断", type="primary")
    if submitted:
        if any(value is None for value in selected.values()):
            st.info("请完成全部 14 题后提交；未作答的题目不计作错误。")
        else:
            results = {q["id"]: selected[q["id"]] == q["answer"] for q in DIAGNOSTIC}
            progress.set_diagnostic(data, results)
            for question_id, correct in results.items():
                if not correct:
                    progress.set_mistake(data, question_id)
            persist(data)
    results = data["diagnostic"]
    if results:
        st.success(recommend(results))
        for question in DIAGNOSTIC:
            if question["id"] in results:
                with st.expander(f"{'答对' if results[question['id']] else '答错'} · {question['id']}"):
                    st.markdown(question["q_md"])
                    st.markdown("**正确答案：** " + question["options"][question["answer"]])
                    st.markdown(question["explain_md"])


def exercise_index():
    index = {}
    for chapter in CHAPTERS.values():
        for item in chapter["written"]:
            index[item["id"]] = ("written", item, chapter)
        for item in chapter["quiz"]:
            index[item["id"]] = ("quiz", item, chapter)
        for item in chapter["code_exercises"]:
            index[item["id"]] = ("code", item, chapter)
    for item in DIAGNOSTIC:
        index[item["id"]] = ("diagnostic", item, None)
    return index


def render_mistakes(data):
    st.title("错题本")
    if not data["mistakes"]:
        st.info("错题本为空。测验或代码检查未通过时自动加入，也可手动收录书面题。")
        return
    index = exercise_index()
    for item_id in list(data["mistakes"]):
        if item_id not in index:
            st.write(f"{item_id}：题目数据尚未加载。")
        else:
            kind, item, chapter = index[item_id]
            st.subheader(f"{item_id} · {chapter_label(chapter['id'] if chapter else item['chapter'])}")
            if kind == "written":
                render_written(item, chapter, data, "notebook")
            elif kind == "quiz":
                render_quiz(item, chapter, data, "notebook")
            elif kind == "code":
                st.markdown(item["prompt_md"])
                st.caption(f"在侧边栏进入第{chapter['num']}章的「手写代码」重做。")
            else:
                st.markdown(item["q_md"])
                for option in item["options"]:
                    st.markdown("- " + option)
                with st.expander("答案与依据"):
                    st.markdown(item["options"][item["answer"]])
                    st.markdown(item["explain_md"])
        if st.button("从错题本移除", key=f"remove-{item_id}"):
            progress.set_mistake(data, item_id, False)
            persist(data)
            st.rerun()
        st.divider()


def render_knowledge_map(data):
    st.title("知识地图")
    st.markdown("箭头表示前置知识流向后续章节；每章均以地基为共同前提。")
    # 地基的共同依赖用说明表达，图中保留一条起始边以避免十条跨图长线。
    positions = {"foundation": (0, 0), "ch1": (1, 0), "ch2": (2, 0), "ch3": (3, 0),
                 "ch4": (4, 1), "ch5": (5, 1), "ch6": (4, -1), "ch7": (6, -1),
                 "ch8": (6, 1), "ch9": (7, 1), "ch10": (7, -1)}
    fig, ax = plt.subplots(figsize=(13, 5))
    for node, dependencies in KNOWLEDGE_MAP.items():
        for dep in dependencies:
            if dep == "foundation" and node != "ch1":
                continue
            ax.annotate("", xy=positions[node], xytext=positions[dep], arrowprops={"arrowstyle": "->", "color": "#a0aabe", "shrinkA": 20, "shrinkB": 20, "connectionstyle": "arc3,rad=.08"})
    for node, (x, y) in positions.items():
        label = "地基" if node == "foundation" else f"第{node[2:]}章"
        opened = node in {"foundation", *(c["id"] for c in CHAPTERS.values())}
        ax.text(x, y, label, ha="center", va="center", bbox={"boxstyle": "round,pad=.5", "fc": "#dcecf5" if opened else "#eeeeee", "ec": "#8aa3bc"})
    ax.set(xlim=(-.5, 7.5), ylim=(-1.6, 1.6))
    ax.axis("off")
    ax.set_title("看什么：先补箭头起点的概念，再学习箭头终点的内容。")
    fig.tight_layout()
    show_figure(fig)
    for node, dependencies in KNOWLEDGE_MAP.items():
        label = chapter_label(node)
        st.markdown(f"**{label}**　← " + ("、".join(chapter_label(dep).split("：")[0] for dep in dependencies) if dependencies else "起点"))


def render_flashcards(data):
    st.title("闪卡自测")
    chapter_ids = [chapter_id for chapter_id in KNOWLEDGE_MAP if any(card["chapter"] == chapter_id for card in FLASHCARDS)]
    chapter = st.selectbox("范围", ["全部", *chapter_ids], format_func=lambda x: x if x == "全部" else chapter_label(x))
    unmastered_only = st.checkbox("只看未掌握", value=False)
    cards = [c for c in FLASHCARDS if (chapter == "全部" or c["chapter"] == chapter)
             and (not unmastered_only or c["id"] not in data["flashcards_mastered"])]
    if not cards:
        st.success("当前范围内的闪卡均已掌握。")
        return
    card_id = st.selectbox("选择闪卡", [c["id"] for c in cards], format_func=lambda x: f"{x} · {chapter_label(next(c for c in cards if c['id'] == x)['chapter']).split('：')[0]}")
    card = next(c for c in cards if c["id"] == card_id)
    st.subheader("先在心中作答")
    st.markdown(card["front"])
    if st.toggle("显示背面", key=f"back-{card_id}"):
        st.markdown(card["back"])
    mastered = card_id in data["flashcards_mastered"]
    checked = st.checkbox("已掌握", value=mastered, key=f"mastered-{card_id}")
    if checked != mastered:
        progress.set_flashcard(data, card_id, checked)
        persist(data)
        st.rerun()
    st.caption(f"当前范围 {len(cards)} 张；全部掌握 {len(set(data['flashcards_mastered']) & {c['id'] for c in FLASHCARDS})}/{len(FLASHCARDS)} 张。")


def render_daily(data):
    st.title("每日一题")
    bank = [(chapter, item) for chapter in CHAPTERS.values() for item in chapter["written"]]
    today = date.today()
    chapter, item = bank[today.toordinal() % len(bank)]
    st.caption(f"{today.isoformat()} · 第{chapter['num']}章 · {item['level']} · 每日按日期固定一题")
    render_written(item, chapter, data, f"daily-{today.isoformat()}")


def main():
    st.set_page_config(page_title="线性代数复习实验室", page_icon="∑", layout="wide")
    if "progress" not in st.session_state:
        try:
            st.session_state["progress"] = progress.load()
        except (ValueError, OSError) as exc:
            st.error(str(exc))
            st.stop()
    data = st.session_state["progress"]
    navigation = ["总览", "地基", *[f"第{i}章" for i in CHAPTER_TITLES], "公式手册", "入学诊断", "错题本", "知识地图", "闪卡自测", "每日一题"]
    with st.sidebar:
        st.title("线性代数复习")
        page = st.radio("导航", navigation, key="navigation")
        st.caption("先推导，再手写，最后检查与复习。进度保存在本地 data/progress.json。")
    if page.startswith("第"):
        number = int(page[1:-1])
        if number in CHAPTERS:
            render_chapter(CHAPTERS[number], data)
        else:
            st.title(f"第{number}章　{CHAPTER_TITLES[number]}")
            st.info("建设中")
            st.markdown("核心公式已收录于公式手册。前置章节见知识地图。")
    elif page == "总览":
        render_overview(data)
    elif page == "地基":
        render_foundation(data)
    elif page == "公式手册":
        render_formula_sheet()
    elif page == "入学诊断":
        render_diagnostic(data)
    elif page == "错题本":
        render_mistakes(data)
    elif page == "知识地图":
        render_knowledge_map(data)
    elif page == "闪卡自测":
        render_flashcards(data)
    elif page == "每日一题":
        render_daily(data)


if __name__ == "__main__":
    main()
