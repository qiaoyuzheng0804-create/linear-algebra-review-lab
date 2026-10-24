"""第一轮验收：数据契约、数学例子、判题、持久化与真实页面流程。"""
from collections import Counter
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
import pytest

from linalg_lab import progress, viz
from linalg_lab.checks import banned_hits, check_code, _same
from linalg_lab.content import CHAPTERS, CHAPTER_TITLES
from linalg_lab.content.foundation import FOUNDATION
from linalg_lab.extras import DIAGNOSTIC, FLASHCARDS, FORMULA_SHEET, KNOWLEDGE_MAP, recommend

pytestmark = pytest.mark.filterwarnings("ignore:Glyph .* missing from font")
EXERCISES = [e for c in CHAPTERS.values() for e in c["code_exercises"]]
SECTIONS = ["动机与直觉", "定义", "核心推导", "适用条件与失效情形", "几何直觉", "与量子力学的连接"]


def execute_reference(exercise):
    namespace = {"np": np}
    exec(exercise["reference"], namespace)
    return namespace[exercise["function_name"]]


def test_registry_integrity():
    assert list(CHAPTERS) == [1, 2, 3]
    assert list(CHAPTER_TITLES) == list(range(1, 11))
    ids = []
    for number, chapter in CHAPTERS.items():
        assert chapter["num"] == number
        assert chapter["id"] == f"ch{number}"
        assert chapter["title"] == CHAPTER_TITLES[number]
        assert set(chapter) >= {"id", "num", "title", "depends_on", "theory_md", "viz", "code_exercises", "written", "quiz"}
        assert all(dep in KNOWLEDGE_MAP for dep in chapter["depends_on"])
        for item in chapter["code_exercises"] + chapter["written"] + chapter["quiz"]:
            ids.append(item["id"])
        assert all(callable(getattr(viz, key)) for key in chapter["viz"])
    assert len(ids) == len(set(ids))


@pytest.mark.parametrize("chapter", CHAPTERS.values(), ids=lambda c: c["id"])
def test_theory_sections_and_derivation(chapter):
    headings = [line.removeprefix("## ") for line in chapter["theory_md"].splitlines() if line.startswith("## ")]
    assert headings == SECTIONS
    assert not any(word in chapter["theory_md"] for word in ("易得", "显然", "不难验证", "类似可证"))
    assert chapter["theory_md"].count("动机：") >= 3
    assert chapter["theory_md"].count("条件与失效：") >= 3
    assert len(chapter["theory_md"]) > 2000


@pytest.mark.parametrize("chapter", CHAPTERS.values(), ids=lambda c: c["id"])
def test_written_and_quiz_contract(chapter):
    assert Counter(x["level"] for x in chapter["written"]) == {"基础": 2, "进阶": 2, "考试风格": 2}
    for item in chapter["written"]:
        steps = [line for line in item["solution_md"].splitlines() if line[:1].isdigit()]
        assert len(steps) >= 4
        assert all("依据" in step for step in steps)
        assert item["pitfalls_md"]
        assert not any(word in item["solution_md"] for word in ("易得", "显然", "不难验证", "类似可证"))
    assert len(chapter["quiz"]) == 5
    for question in chapter["quiz"]:
        assert len(question["options"]) == 4
        assert 0 <= question["answer"] < 4
        assert question["explain_md"]


@pytest.mark.parametrize("exercise", EXERCISES, ids=lambda e: e["id"])
def test_reference_passes_random_checker(exercise):
    assert "TODO" in exercise["starter"]
    assert banned_hits(exercise["reference"], exercise["banned"]) == []
    passed, feedback = check_code(exercise["reference"], exercise, seed=1083)
    assert passed, feedback
    assert len(feedback) >= 20


@pytest.mark.parametrize("exercise", EXERCISES, ids=lambda e: e["id"])
def test_direct_library_cheating_rejected(exercise):
    source = f"def {exercise['function_name']}(*args):\n    return np.linalg.solve(args[0], args[1])\n"
    passed, feedback = check_code(source, exercise)
    assert not passed
    assert any("np.linalg.solve" in message for message in feedback)


@pytest.mark.parametrize("source, expected", [
    ("import numpy as other\nf=other.linalg.svd", "np.linalg.svd"),
    ("from numpy.linalg import inv as inverse\nx=inverse(A)", "np.linalg.inv"),
    ("import scipy.linalg as la\nx=la.solve(A,b)", "scipy"),
    ("f=np.linalg.eig\ng=f\nx=g(A)", "np.linalg.eig"),
    ("import numpy.linalg as la\nx=la.det(A)", "np.linalg.det"),
    ("from numpy import linalg as la\nx=la.qr(A)", "np.linalg.qr"),
])
def test_cheating_aliases(source, expected):
    hits = banned_hits(source, ["np.linalg", "scipy"])
    assert any(expected in hit for hit in hits)


@pytest.mark.parametrize("source", ["return A @ B", "return np.matmul(A, B)", "return A.dot(B)"])
def test_matmul_shortcuts_rejected(source):
    exercise = CHAPTERS[2]["code_exercises"][0]
    assert banned_hits("def matmul(A, B):\n    " + source, exercise["banned"])


def test_checker_rejects_wrong_output_and_handles_failure():
    exercise = EXERCISES[0]
    passed, feedback = check_code("def independent_group(A):\n    return False, np.array([], dtype=int)", exercise)
    assert not passed and any("不符" in line for line in feedback)
    assert not check_code("def broken(:", exercise)[0]
    assert not check_code("x = 1", exercise)[0]
    assert not check_code("def independent_group(A):\n    raise ValueError('错误')", exercise)[0]
    timeout_result = check_code("while True:\n    pass", exercise, timeout=1)
    assert not timeout_result[0] and "超过" in timeout_result[1][0]


def test_checker_requires_twenty_cases():
    exercise = dict(EXERCISES[0], make_cases="def make_cases(rng):\n    return [(np.eye(2),)]")
    passed, feedback = check_code(exercise["reference"], exercise)
    assert not passed and "不足 20" in feedback[0]


def test_output_comparison_has_no_broadcasting_false_positive():
    assert not _same(1, np.ones((2, 2)))
    assert not _same([1, 2], (1, 2))
    assert _same((1, np.array([2, 3j])), (1, np.array([2, 3j])))
    assert not _same(np.nan, np.nan)


def test_reference_rank_and_complex_matmul_examples():
    independent = execute_reference(EXERCISES[0])
    related, pivots = independent(np.array([[1., 2., 0.], [2., 4., 1.]]))
    assert related and pivots.tolist() == [0, 2]
    matmul = execute_reference(CHAPTERS[2]["code_exercises"][0])
    A = np.array([[1+1j, 2], [3, -1j]])
    B = np.array([[2, 0], [1j, 1]])
    assert np.allclose(matmul(A, B), A @ B)
    in_basis = execute_reference(CHAPTERS[2]["code_exercises"][1])
    A, P = np.diag([2., 1.]), np.array([[1., 1.], [1., -1.]])
    assert np.allclose(in_basis(A, P), [[1.5, .5], [.5, 1.5]])
    with pytest.raises(ValueError):
        in_basis(A, np.ones((2, 2)))


@pytest.mark.parametrize("seed", [17, 209])
def test_gaussian_solution_certificates(seed):
    exercise = CHAPTERS[3]["code_exercises"][0]
    solve = execute_reference(exercise)
    namespace = {"np": np}
    exec(exercise["make_cases"], namespace)
    for A, b in namespace["make_cases"](np.random.default_rng(seed)):
        status, x, Z = solve(A, b)
        if status:
            assert np.allclose(A @ x, b)
            assert np.allclose(A @ Z, 0)
            assert Z.shape[0] == A.shape[1]
            assert (Z.shape[1] == 0) == (status == 1)
            if Z.shape[1]:
                # 每个自由变量给出单位坐标，独立性可不借助秩库验证。
                assert all(any(np.allclose(row, np.eye(Z.shape[1])[j]) for row in Z) for j in range(Z.shape[1]))
        else:
            assert np.all(x == 0) and Z.shape == (A.shape[1], 0)
            assert not np.allclose(A @ x, b)


def test_progress_roundtrip_and_queries(tmp_path):
    path = tmp_path / "data" / "progress.json"
    data = progress.load(path)
    assert data == progress.empty_progress()
    chapter = CHAPTERS[1]
    assert progress.completion(data, chapter) == 0
    progress.mark_theory(data, "ch1")
    for _ in range(2):
        progress.mark_code(data, "ch1", "ch1-c1")
        progress.mark_written(data, "ch1", "ch1-w1")
        progress.mark_quiz(data, "ch1", "ch1-q1", True)
        progress.set_mistake(data, "ch1-q2")
        progress.set_flashcard(data, "card-1")
    progress.set_diagnostic(data, {"diag-1": False})
    progress.save(data, path)
    restored = progress.load(path)
    assert restored == data
    assert restored["chapters"]["ch1"]["quiz"] == 1
    assert restored["chapters"]["ch1"]["code_ex"] == ["ch1-c1"]
    assert restored["mistakes"] == ["ch1-q2"]
    assert restored["flashcards_mastered"] == ["card-1"]
    assert progress.completion(restored, chapter) == pytest.approx(4 / 13)
    progress.mark_quiz(restored, "ch1", "ch1-q1", False)
    progress.mark_written(restored, "ch1", "ch1-w1", False)
    progress.set_mistake(restored, "ch1-q2", False)
    progress.set_flashcard(restored, "card-1", False)
    assert restored["chapters"]["ch1"]["quiz"] == 0
    assert restored["chapters"]["ch1"]["written"] == []
    assert not restored["mistakes"] and not restored["flashcards_mastered"]
    assert not list(path.parent.glob(".progress-*"))


def test_progress_does_not_overwrite_corrupt_file(tmp_path):
    path = tmp_path / "progress.json"
    path.write_text("{", encoding="utf-8")
    with pytest.raises(ValueError, match="进度文件格式错误"):
        progress.load(path)
    assert path.read_text() == "{"


@pytest.mark.parametrize("function, args", [
    (viz.grid_transform, (np.array([[1., .5], [0., 1.]]),)),
    (viz.det_area, (np.array([[1., 2.], [0., 0.]]),)),
    (viz.eigen_demo, (np.array([[2., 0.], [0., 1.]]),)),
    (viz.eigen_demo, (np.array([[0., -1.], [1., 0.]]),)),
    (viz.gram_schmidt, ([[1, 1], [2, 0]], 2)),
    (viz.gram_schmidt, ([[1, 1], [2, 2]], 2)),
    (viz.svd_ellipse, (np.array([[1., 2.], [0., 0.]]),)),
    (viz.matrix_power, (np.array([[.8, .1], [0., .5]]), 5)),
    (viz.tensor_dims, ()),
])
def test_visualizations_return_figures(function, args):
    figure = function(*args)
    try:
        assert isinstance(figure, Figure)
        figure.canvas.draw()
        titles = [ax.get_title() for ax in figure.axes] + ([figure._suptitle.get_text()] if figure._suptitle else [])
        assert any("看什么" in title for title in titles)
    finally:
        plt.close(figure)


def test_extras_complete_and_acyclic():
    assert len(FOUNDATION) == 4
    assert all(len(section["selfcheck"]) == 2 for section in FOUNDATION)
    assert all(q["q"] and q["answer_md"] for section in FOUNDATION for q in section["selfcheck"])
    assert list(FORMULA_SHEET) == ["地基", *[f"第{i}章" for i in range(1, 11)]]
    assert all(formula and condition for rows in FORMULA_SHEET.values() for formula, condition in rows)
    assert len(DIAGNOSTIC) == 14
    assert {q["chapter"] for q in DIAGNOSTIC} == set(KNOWLEDGE_MAP)
    assert all(len(q["options"]) == 4 and 0 <= q["answer"] < 4 and q["explain_md"] for q in DIAGNOSTIC)
    assert len(FLASHCARDS) >= 30
    assert len({card["id"] for card in FLASHCARDS}) == len(FLASHCARDS)
    assert all(card["front"] and card["back"] and card["chapter"] in KNOWLEDGE_MAP for card in FLASHCARDS)
    def visit(node, ancestors):
        assert node not in ancestors
        for dep in KNOWLEDGE_MAP[node]:
            visit(dep, ancestors | {node})
    for node in KNOWLEDGE_MAP:
        visit(node, set())
    assert all("foundation" in deps for node, deps in KNOWLEDGE_MAP.items() if node != "foundation")


def test_diagnostic_recommendations_follow_error_distribution():
    assert "尚未完成" in recommend({})
    assert "全部正确" in recommend({q["id"]: True for q in DIAGNOSTIC})
    assert "从地基开始" in recommend({"diag-1": False, "diag-6": False})
    assert "从第3章" in recommend({"diag-6": False, "diag-7": False})
    assert "从第1章" in recommend({"diag-4": False, "diag-6": False, "diag-7": False})
    assert "建设中" in recommend({"diag-10": False})
    assert "1/14" in recommend({"diag-1": True})


def test_modules_import_without_starting_app():
    assert callable(importlib.import_module("app").main)


def test_app_navigation_smoke(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(progress, "DEFAULT_PATH", tmp_path / "progress.json")
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run(timeout=20)
    assert not app.exception
    choices = app.sidebar.radio[0].options
    assert len(choices) == 18
    for choice in choices:
        app.sidebar.radio[0].set_value(choice).run(timeout=20)
        assert not app.exception, f"页面 {choice} 异常：{app.exception}"


def test_quiz_wrong_answer_stays_in_notebook(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    path = tmp_path / "progress.json"
    monkeypatch.setattr(progress, "DEFAULT_PATH", path)
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run()
    app.sidebar.radio[0].set_value("第1章").run()
    app.radio(key="chapter-choice-ch1-q1").set_value(1)
    app.button(key="FormSubmitter:chapter-form-ch1-q1-提交答案").click().run()
    assert not app.exception
    assert "ch1-q1" in progress.load(path)["mistakes"]
    app.radio(key="chapter-choice-ch1-q1").set_value(0)
    app.button(key="FormSubmitter:chapter-form-ch1-q1-提交答案").click().run()
    assert progress.load(path)["chapters"]["ch1"]["quiz"] == 1
    app.sidebar.radio[0].set_value("错题本").run()
    app.button(key="remove-ch1-q1").click().run()
    assert "ch1-q1" not in progress.load(path)["mistakes"]


def test_code_submission_and_shared_daily_progress(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    path = tmp_path / "progress.json"
    monkeypatch.setattr(progress, "DEFAULT_PATH", path)
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run()
    app.sidebar.radio[0].set_value("第1章").run()
    exercise = CHAPTERS[1]["code_exercises"][0]
    app.text_area(key="source-ch1-c1").set_value(exercise["reference"])
    app.button(key="run-ch1-c1").click().run(timeout=20)
    assert not app.exception
    assert "ch1-c1" in progress.load(path)["chapters"]["ch1"]["code_ex"]
    app.sidebar.radio[0].set_value("每日一题").run()
    checkbox = next(x for x in app.checkbox if x.label == "已独立完成并核对")
    item_id = checkbox.key.split("-done-")[1]
    checkbox.check().run()
    cid = item_id.split("-")[0]
    assert item_id in progress.load(path)["chapters"][cid]["written"]
    app.sidebar.radio[0].set_value(f"第{cid[2:]}章").run()
    assert app.checkbox(key=f"chapter-done-{item_id}").value


@pytest.mark.parametrize("source", [
    "import numpy as n\ndef matmul(A, B):\n    return n.matmul(A, B)\ndef unused():\n    import numpy.random as n\n",
    "def matmul(A, B):\n    n = np\n    result = n.matmul(A, B)\n    n = np.random\n    return result\n",
    "def matmul(A, B):\n    n: object = np\n    return n.matmul(A, B)\n",
    "def setup():\n    global n\n    import numpy as n\nsetup()\ndef matmul(A, B):\n    return n.matmul(A, B)\n",
    "def helper(lib, A, B):\n    return lib.matmul(A, B)\ndef matmul(A, B):\n    return helper(np, A, B)\n",
    "def matmul(A, B):\n    lib, dummy = np, None\n    return lib.matmul(A, B)\n",
])
def test_banned_bindings_are_not_masked_by_other_scopes_or_later_assignments(source):
    exercise = CHAPTERS[2]["code_exercises"][0]
    passed, feedback = check_code(source, exercise)
    assert not passed
    assert any("np.matmul" in message for message in feedback)


def test_chapter_progress_updates_in_current_interaction(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(progress, "DEFAULT_PATH", tmp_path / "progress.json")
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run()
    app.sidebar.radio[0].set_value("第1章").run()
    assert app.get("progress")[0].value == 0
    app.checkbox(key="theory-ch1").check().run()
    assert app.get("progress")[0].value == int(100 / 13)
    app.checkbox(key="chapter-done-ch1-w1").check().run()
    assert app.get("progress")[0].value == int(200 / 13)


def test_code_draft_survives_page_navigation(tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(progress, "DEFAULT_PATH", tmp_path / "progress.json")
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run()
    app.sidebar.radio[0].set_value("第1章").run()
    draft = "def independent_group(A):\n    # 我的未完成代码草稿\n    R = np.array(A, copy=True)\n"
    app.text_area(key="source-ch1-c1").set_value(draft).run()
    app.sidebar.radio[0].set_value("公式手册").run()
    app.sidebar.radio[0].set_value("第1章").run()
    assert app.text_area(key="source-ch1-c1").value == draft
