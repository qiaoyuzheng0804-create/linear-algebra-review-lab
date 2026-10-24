"""第二轮验收：内容契约、独立数值凭据、反作弊及真实章节交互。"""
from collections import Counter
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
import pytest

from linalg_lab import progress, viz
from linalg_lab.checks import banned_hits, check_code
from linalg_lab.content import CHAPTERS, CHAPTER_TITLES
from linalg_lab.extras import DIAGNOSTIC, FORMULA_SHEET, KNOWLEDGE_MAP, recommend

pytestmark = pytest.mark.filterwarnings("ignore:Glyph .* missing from font")
CHAPTER_NUMBERS = (4, 5, 6, 7)
EXERCISES = [e for n in CHAPTER_NUMBERS for e in CHAPTERS[n]["code_exercises"]]
SECTIONS = ["动机与直觉", "定义", "核心推导", "适用条件与失效情形", "几何直觉", "与量子力学的连接"]


def reference_and_cases(exercise, seed):
    namespace = {"np": np}
    exec(exercise["reference"], namespace)
    exec(exercise["make_cases"], namespace)
    return namespace[exercise["function_name"]], namespace["make_cases"](np.random.default_rng(seed))


@pytest.mark.parametrize("number", CHAPTER_NUMBERS)
def test_round2_content_contract(number):
    c = CHAPTERS[number]
    assert c["id"] == f"ch{number}" and c["num"] == number
    assert c["title"] == CHAPTER_TITLES[number]
    assert set(c["depends_on"]) <= set(KNOWLEDGE_MAP[c["id"]])
    theory = c["theory_md"]
    assert re.findall(r"^## (.+)$", theory, re.M) == SECTIONS
    assert len(theory) >= 3200
    assert theory.count("动机：") >= 4
    assert theory.count("条件与失效：") >= 4
    assert "手算算例" in theory
    assert not any(w in theory for w in ("易得", "显然", "不难验证", "类似可证"))
    assert len(c["written"]) == 6
    assert Counter(w["level"] for w in c["written"]) == {"基础": 2, "进阶": 2, "考试风格": 2}
    for w in c["written"]:
        steps = re.findall(r"^\d+\. (.+)$", w["solution_md"], re.M)
        assert len(steps) >= 4 and all("依据" in step for step in steps)
        assert w["problem_md"] and w["pitfalls_md"]
    assert len(c["quiz"]) == 5
    for q in c["quiz"]:
        assert q["q_md"] and q["explain_md"]
        assert len(q["options"]) == len(set(q["options"])) == 4
        assert isinstance(q["answer"], int) and 0 <= q["answer"] < 4
    assert c["code_exercises"] and c["viz"]
    ids = [x["id"] for key in ("written", "quiz", "code_exercises") for x in c[key]]
    assert len(ids) == len(set(ids))
    for e in c["code_exercises"]:
        assert all(e[key] for key in ("starter", "prompt_md", "banned", "reference", "make_cases", "function_name"))


@pytest.mark.parametrize("number", CHAPTER_NUMBERS)
def test_round2_figures(number):
    arguments = {
        "det_area": (np.array([[2., 1.], [1., 3.]]),),
        "eigen_demo": (np.array([[2., 1.], [0., 1.]]),),
        "gram_schmidt": ([[1., 1.], [1., 0.]], 2),
        "matrix_power": (np.array([[.8, .2], [0., .5]]), 6),
    }
    for key in CHAPTERS[number]["viz"]:
        fig = getattr(viz, key)(*arguments[key])
        try:
            assert isinstance(fig, Figure)
            fig.canvas.draw()
        finally:
            plt.close(fig)


@pytest.mark.parametrize("exercise", EXERCISES, ids=lambda e: e["id"])
def test_round2_reference_passes_checker(exercise):
    assert not banned_hits(exercise["reference"], exercise["banned"])
    passed, feedback = check_code(exercise["reference"], exercise, seed=71027)
    assert passed, feedback
    assert len(feedback) >= 20


CHEAT_VARIANTS = [
    "np.linalg.{api}(A)",
    "import numpy as n\nn.linalg.{api}(A)",
    "import numpy.linalg as la\nla.{api}(A)",
    "from numpy import linalg as la\nla.{api}(A)",
    "from numpy.linalg import {api} as shortcut\nshortcut(A)",
    "f = np.linalg.{api}\ng = f\ng(A)",
    "def helper(lib, A):\n    return lib.{api}(A)\nhelper(np.linalg, A)",
    "lib, extra = np.linalg, None\nlib.{api}(A)",
    "lib: object = np.linalg\nlib.{api}(A)",
    "(lib := np.linalg).{api}(A)",
    "from numpy.linalg import *\n{api}(A)",
    "getattr(np.linalg, '{api}')(A)",
    "np.linalg.__dict__['{api}'](A)",
    "import scipy.linalg as la\nla.{api}(A)",
]


@pytest.mark.parametrize("number,api", [(4, "det"), (5, "eig"), (6, "qr"), (7, "eig")])
@pytest.mark.parametrize("variant", CHEAT_VARIANTS)
def test_round2_cheat_variants_rejected(number, api, variant):
    e = CHAPTERS[number]["code_exercises"][0]
    source = variant.format(api=api)
    assert banned_hits(source, e["banned"])
    passed, feedback = check_code(source, e)
    assert not passed and all("作弊检查失败" in line for line in feedback)


@pytest.mark.parametrize("number,source", [
    (4, "np.linalg.slogdet(A)"),
    (5, "from numpy.linalg import eigvals as roots\nroots(A)"),
    (5, "np.linalg.eigh(A)"),
    (5, "np.linalg.eigvalsh(A)"),
    (6, "np.linalg.cholesky(A)"),
    (6, "np.linalg.svd(A)"),
    (7, "np.linalg.matrix_power(A, 3)"),
    (7, "from scipy.linalg import expm\nexpm(A)"),
    (7, "np.roots([1, 2, 3])"),
    (7, "from numpy.polynomial.polynomial import polyroots\npolyroots([1, 2, 3])"),
    (5, "from numpy.lib._polynomial_impl import eigvals as spectrum\nspectrum(A)"),
    (7, "from numpy.lib._polynomial_impl import roots as shortcut\nshortcut([1, 2, 3])"),
    (7, "from numpy.lib._polynomial_impl import eigvals as spectrum\nspectrum(A)"),
    (7, "from numpy.lib.polynomial import roots as shortcut\nshortcut([1, 2, 3])"),
    (7, "from numpy.lib._polynomial_impl import *\nroots([1, 2, 3])"),
    (4, "from numpy.lib._polynomial_impl import poly as characteristic\ncharacteristic(A)"),
    (5, "np.poly(A)"),
    (7, "np.poly(A)"),
    (7, "np.matrix(A).I"),
    (7, "from numpy import asmatrix as wrap\nwrap(A).I"),
    (7, "from numpy.matrixlib.defmatrix import matrix as wrap\nwrap(A)**3"),
    (7, "from numpy.matlib import eye as factory\nfactory(2).I"),
    (7, "def hidden_inverse(M):\n    return M.getI()"),
])
def test_round2_alternative_library_shortcuts_rejected(number, source):
    e = CHAPTERS[number]["code_exercises"][0]
    assert banned_hits(source, e["banned"])


@pytest.mark.parametrize("import_line,function", [
    ("from numpy.lib._polynomial_impl import roots as shortcut", "shortcut"),
    ("from numpy.lib._polynomial_impl import *", "roots"),
])
def test_reexport_cheat_cannot_pass_numeric_checker(import_line, function):
    e = CHAPTERS[7]["code_exercises"][0]
    source = re.sub(r"    root = .*\n    values = .*\n",
                    f"    values = {function}([1, -trace, determinant])\n", e["reference"])
    source = import_line + "\n" + source
    passed, feedback = check_code(source, e)
    assert not passed and all("作弊检查失败" in line for line in feedback)


def test_matrix_inverse_wrapper_cannot_pass_numeric_checker():
    e = CHAPTERS[7]["code_exercises"][0]
    source = e["reference"]
    start = source.index("    c = np.array(")
    end = source.index("    return P @", start)
    source = source[:start] + "    c = np.asarray(np.matrix(P).I @ x).reshape(2)\n" + source[end:]
    passed, feedback = check_code(source, e)
    assert not passed and all("作弊检查失败" in line for line in feedback)


@pytest.mark.parametrize("seed", [29, 20261024])
def test_determinants_independent_certificates(seed):
    f, cases = reference_and_cases(EXERCISES[0], seed)
    for (A,) in cases:
        before = A.copy()
        expanded, eliminated = f(A)
        target = np.linalg.det(A)
        assert np.allclose(expanded, target)
        assert np.allclose(eliminated, target)
        assert np.array_equal(A, before)
    assert f(np.array([[0., 2.], [3., 4.]])) == (-6., -6.)
    with pytest.raises(ValueError):
        f(np.ones((2, 3)))


@pytest.mark.parametrize("seed", [41, 875])
def test_power_iteration_eigenpair_certificates(seed):
    f, cases = reference_and_cases(EXERCISES[1], seed)
    for A, x0, steps in cases:
        before = x0.copy()
        value, x = f(A, x0, steps)
        assert x.shape == x0.shape
        assert np.isclose(np.vdot(x, x), 1)
        assert np.array_equal(x0, before)
        pivot = x[np.argmax(np.abs(x))]
        assert np.isclose(pivot.imag, 0) and pivot.real > 0
        assert np.allclose(value, np.vdot(x, A @ x))
        if steps >= 15:
            values = np.linalg.eigvalsh(A)
            target = values[np.argmax(np.abs(values))]
            assert np.allclose(value, target)
            assert np.linalg.norm(A @ x-value*x) < 1e-9
    with pytest.raises(ValueError):
        f(np.eye(2), np.zeros(2), 10)
    value, x = f(np.diag([-3., 1.]), np.array([0., 1.]), 80)
    assert value == 1 and np.allclose(x, [0, 1])


@pytest.mark.parametrize("seed", [15, 707])
def test_qr_independent_certificates(seed):
    f, cases = reference_and_cases(EXERCISES[2], seed)
    for (A,) in cases:
        before = A.copy()
        Q, R = f(A)
        m, n = A.shape
        assert Q.shape == (m, n) and R.shape == (n, n)
        assert np.allclose(Q.conj().T @ Q, np.eye(n))
        assert np.allclose(Q @ R, A)
        assert np.allclose(R, np.triu(R))
        assert np.all(np.diag(R).real > 0) and np.allclose(np.diag(R).imag, 0)
        assert np.array_equal(A, before)
    with pytest.raises(ValueError):
        f(np.array([[1., 2.], [2., 4.]]))
    with pytest.raises(ValueError):
        f(np.ones((1, 2)))


def exponential_by_series(A, t):
    # 独立的级数求值，不复用待测的二次根公式或谱分解。
    total = np.eye(len(A), dtype=complex)
    term = total.copy()
    for k in range(1, 80):
        term = (term @ A) * (t / k)
        total += term
    return total


@pytest.mark.parametrize("seed", [13, 906])
def test_spectral_action_independent_certificates(seed):
    f, cases = reference_and_cases(EXERCISES[3], seed)
    for A, x, k, t in cases:
        before_A, before_x = A.copy(), x.copy()
        power, exponential = f(A, x, k, t)
        assert power.shape == exponential.shape == (2,)
        expected = x.copy()
        for _ in range(k):
            expected = A @ expected
        assert np.allclose(power, expected)
        assert np.allclose(exponential, exponential_by_series(A, t) @ x)
        assert np.array_equal(A, before_A) and np.array_equal(x, before_x)
    power, exponential = f(np.array([[0., -1.], [1., 0.]]), np.array([1., 0.]), 4, 1.)
    assert np.allclose(power, [1, 0])
    assert np.allclose(exponential, [np.cos(1), np.sin(1)])
    with pytest.raises(ValueError):
        f(np.array([[1., 1.], [0., 1.]]), np.ones(2), 2, 1.)
    # 特征值合法且本征基条件数为 1，原 trace²−4det 会相消为零。
    diagonal = np.array([1., 1.+1e-8])
    power, exponential = f(np.diag(diagonal), np.ones(2), 2, .5)
    assert np.allclose(power, diagonal**2, rtol=1e-12, atol=1e-14)
    assert np.allclose(exponential, np.exp(.5*diagonal), rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("number,fragments", [
    (4, [r"\det(AB)=\det A\det B", r"x_j=\frac{\det A_j(b)}{\det A}"]),
    (5, [r"p_A(A)=0", r"p_A(z)=\det(zI-A)"]),
    (6, [r"Q^\dagger Q=I", r"QQ^\dagger x"]),
    (7, [r"A^k=P\Lambda^kP^{-1}", r"x(t)=e^{tA}x(0)"]),
])
def test_formula_sheet_matches_theory(number, fragments):
    rows = FORMULA_SHEET[f"第{number}章"]
    theory = CHAPTERS[number]["theory_md"]
    sheet = " ".join(formula for formula, _ in rows)
    compact = lambda s: re.sub(r"\s+", "", s)
    for fragment in fragments:
        assert compact(fragment) in compact(sheet)
        # 初值记号 x0 与 x(0) 表达同一条件；检查传播公式与初值说明。
        if fragment == r"x(t)=e^{tA}x(0)":
            assert r"x(t)=\sum_jc_je^{\lambda_jt}v_j=e^{tA}x_0" in theory
            assert "x(0)=x_0" in theory
        else:
            assert compact(fragment) in compact(theory)
    for _, condition in rows:
        for chapter, tag in re.findall(r"式 \((\d+)\.(\d+)\)", condition):
            assert int(chapter) == number
            assert rf"\tag{{{chapter}.{tag}}}" in theory
    if number == 4:
        assert any("det A 非零" in condition for _, condition in rows)
    if number == 7:
        assert any("AB=BA" in condition for _, condition in rows)
        assert not any("Sigma" in formula for formula, _ in rows)


def test_diagnostic_numbering_and_open_chapter_recommendations():
    assert {q["chapter"] for q in DIAGNOSTIC} == set(KNOWLEDGE_MAP)
    assert all(q["chapter"] == "ch9" for q in DIAGNOSTIC if "SVD" in q["q_md"])
    assert all(q["chapter"] == "ch7" for q in DIAGNOSTIC if "常系数系统" in q["q_md"])
    for q in DIAGNOSTIC:
        if q["chapter"] in {f"ch{n}" for n in CHAPTER_NUMBERS}:
            assert "建设中" not in recommend({q["id"]: False})


@pytest.mark.parametrize("number", CHAPTER_NUMBERS)
def test_new_chapter_page_interaction(number, tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(progress, "DEFAULT_PATH", tmp_path / "progress.json")
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / "app.py")).run(timeout=20)
    app.sidebar.radio[0].set_value(f"第{number}章").run(timeout=20)
    assert not app.exception
    assert app.title[0].value == f"第{number}章　{CHAPTER_TITLES[number]}"
    assert any("## 核心推导" in m.value for m in app.markdown)
    assert not any(m.value == "建设中" for m in app.info)
    assert len(app.text_area) == len(CHAPTERS[number]["code_exercises"])
    assert len(app.radio) == 6  # 导航一项及五道测验
    assert len(app.get("image")) == 1
    if number in (4, 5):
        key = CHAPTERS[number]["viz"][0]
        app.number_input(key=f"ch{number}-{key}-a₁₁").set_value(2.).run()
    elif number == 6:
        app.slider(key="ch6-gs").set_value(2).run()
    else:
        app.slider(key="ch7-power").set_value(9).run()
    assert not app.exception and len(app.get("image")) == 1
    q = CHAPTERS[number]["quiz"][0]
    app.radio(key=f"chapter-choice-{q['id']}").set_value(q["answer"])
    app.button(key=f"FormSubmitter:chapter-form-{q['id']}-提交答案").click().run()
    assert not app.exception
    assert progress.load(tmp_path / "progress.json")["chapters"][f"ch{number}"]["quiz"] == 1
