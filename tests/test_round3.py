"""第三轮验收：十章内容闭环、独立数值凭据、防作弊与真实页面交互。"""
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
from linalg_lab.extras import DIAGNOSTIC, FLASHCARDS, FORMULA_SHEET, KNOWLEDGE_MAP, recommend

pytestmark = pytest.mark.filterwarnings("ignore:Glyph .* missing from font")
NUMBERS = (8, 9, 10)
EXERCISES = [CHAPTERS[n]["code_exercises"][0] for n in NUMBERS]
SECTIONS = ["动机与直觉", "定义", "核心推导", "适用条件与失效情形", "几何直觉", "与量子力学的连接"]


def load_exercise(number, seed=328):
    e = CHAPTERS[number]["code_exercises"][0]
    namespace = {"np": np}
    exec(e["reference"], namespace)
    exec(e["make_cases"], namespace)
    return namespace[e["function_name"]], namespace["make_cases"](np.random.default_rng(seed))


@pytest.mark.parametrize("number", NUMBERS)
def test_round3_chapter_contract(number):
    c = CHAPTERS[number]
    assert c["id"] == f"ch{number}" and c["num"] == number
    assert c["title"] == CHAPTER_TITLES[number]
    assert set(c["depends_on"]) <= set(KNOWLEDGE_MAP[c["id"]])
    theory = c["theory_md"]
    assert re.findall(r"^## (.+)$", theory, re.M) == SECTIONS
    assert len(theory) >= 6000
    if number == 8:
        assert len(theory) >= max(len(CHAPTERS[n]["theory_md"]) for n in (6, 7))
    assert theory.count("动机：") >= 5
    assert theory.count("条件与失效：") >= 5
    assert "完整手算算例" in theory
    assert not re.search("易得|显然|不难验证|类似可证", theory)
    assert len(c["written"]) == 6
    assert Counter(w["level"] for w in c["written"]) == {"基础": 2, "进阶": 2, "考试风格": 2}
    for w in c["written"]:
        steps = re.findall(r"^\d+\. (.+)$", w["solution_md"], re.M)
        assert len(steps) >= 4 and all("依据" in step for step in steps)
        assert w["problem_md"] and w["pitfalls_md"]
        assert not re.search("易得|显然|不难验证|类似可证", w["solution_md"])
    assert len(c["quiz"]) == 5
    for q in c["quiz"]:
        assert len(q["options"]) == len(set(q["options"])) == 4
        assert 0 <= q["answer"] < 4 and q["q_md"] and q["explain_md"]
    assert c["viz"] and c["code_exercises"]
    for e in c["code_exercises"]:
        assert all(e[k] for k in ("starter", "banned", "reference", "make_cases", "prompt_md", "function_name"))


@pytest.mark.parametrize("e", EXERCISES, ids=lambda e: e["id"])
def test_round3_reference_passes_checker(e):
    assert banned_hits(e["reference"], e["banned"]) == []
    passed, feedback = check_code(e["reference"], e, seed=831019)
    assert passed, feedback
    assert len(feedback) >= 20


@pytest.mark.parametrize("e", EXERCISES, ids=lambda e: e["id"])
def test_round3_starter_is_not_a_solution(e):
    passed, feedback = check_code(e["starter"], e)
    assert not passed and any("不符" in line for line in feedback)


CHEATS = [
    "np.{module}.{api}(A)",
    "import numpy as n\nn.{module}.{api}(A)",
    "from numpy.{module} import {api} as f\nf(A)",
    "import numpy.{module} as lib\nlib.{api}(A)",
    "f=np.{module}.{api}\ng=f\ng(A)",
    "def helper(lib):\n    return lib.{api}(A)",
    "lib, _ = np.{module}, 0\nlib.{api}(A)",
    "lib: object = np.{module}\nlib.{api}(A)",
    "(lib := np.{module}).{api}(A)",
    "getattr(np.{module}, '{api}')(A)",
    "np.{module}.__dict__['{api}'](A)",
    "from numpy.{module} import *\n{api}(A)",
]


@pytest.mark.parametrize("number,module,api", [(8, "linalg", "eigh"), (8, "linalg", "eig"),
                                                (9, "linalg", "svd"), (9, "linalg", "eig")])
@pytest.mark.parametrize("variant", CHEATS)
def test_round3_spectral_cheat_variants(number, module, api, variant):
    e = CHAPTERS[number]["code_exercises"][0]
    source = variant.format(module=module, api=api)
    assert banned_hits(source, e["banned"])
    passed, feedback = check_code(source, e)
    assert not passed and all("作弊检查失败" in line for line in feedback)


@pytest.mark.parametrize("number,source", [
    (8, "from numpy.lib._polynomial_impl import eigvals as f\nf(A)"),
    (9, "np.linalg.eigh(A)"),
    (9, "np.linalg.eigvalsh(A)"),
    (9, "np.linalg.pinv(A)"),
    (9, "np.roots([1,2,3])"),
    (9, "np.matrix(A).I"),
    (9, "from scipy.linalg import svd as f\nf(A)"),
    (10, "np.kron(A,B)"),
    (10, "import numpy as n\nf=n.kron\ng=f\ng(A,B)"),
    (10, "from numpy import kron as blocks\nblocks(A,B)"),
    (10, "from numpy.lib._shape_base_impl import kron as f\nf(A,B)"),
    (10, "np.lib.shape_base.kron(A,B)"),
    (10, "def helper(lib):\n    return lib.kron(A,B)"),
    (10, "lib, _ = np, 0\nlib.kron(A,B)"),
    (10, "getattr(np, 'kron')(A,B)"),
    (10, "np.__dict__['kron'](A,B)"),
    (10, "np.outer(A,B)"),
    (10, "np.multiply.outer(A,B)"),
    (10, "np.einsum('ij,kj->ik', A,B)"),
    (10, "np.tensordot(A,B,axes=0)"),
    (10, "rho.reshape(2,2,2,2).trace(axis1=1,axis2=3)"),
    (10, "np.diagonal(rho)"),
    (10, "import qutip\nqutip.ptrace(rho, 0)"),
])
def test_round3_shortcut_surfaces(number, source):
    e = CHAPTERS[number]["code_exercises"][0]
    assert banned_hits(source, e["banned"])
    assert not check_code(source, e)[0]


@pytest.mark.parametrize("seed", [23, 114])
def test_spectral_audit_independent_certificates(seed):
    f, cases = load_exercise(8, seed)
    for H, values, Q, psi in cases:
        originals = [a.copy() for a in (H, values, Q, psi)]
        rebuilt, errors, weights = f(H, values, Q, psi)
        expected_rebuilt = Q @ np.diag(values) @ Q.conj().T
        assert np.allclose(rebuilt, expected_rebuilt)
        assert np.allclose(errors, [np.linalg.norm(H-H.conj().T),
                                   np.linalg.norm(Q.conj().T@Q-np.eye(len(values))),
                                   np.linalg.norm(H@Q-Q@np.diag(values)),
                                   np.linalg.norm(H-expected_rebuilt)])
        expected_weights = []
        for value in np.unique(values):
            columns = Q[:, values == value]
            expected_weights.append(np.linalg.norm(columns @ (columns.conj().T @ psi))**2)
        assert np.allclose(weights, expected_weights)
        if np.max(errors) < 1e-10:
            assert np.isclose(sum(weights), 1)
            assert np.allclose(weights @ np.unique(values), np.vdot(psi, H@psi))
        assert all(np.array_equal(a, b) for a, b in zip((H, values, Q, psi), originals))
    # 简并能量只有一个结果；非单位列不能被重构一项掩盖。
    _, errors, weights = f(np.eye(2), np.ones(2), np.eye(2), np.array([1., 0.]))
    assert weights.shape == (1,) and np.allclose(weights, [1]) and np.allclose(errors, 0)
    _, errors, _ = f(np.zeros((2, 2)), np.zeros(2), np.ones((2, 2)), np.array([1., 0.]))
    assert errors[1] > 0 and errors[3] == 0
    with pytest.raises(ValueError):
        f(np.eye(2), np.ones(3), np.eye(2), np.ones(2))


@pytest.mark.parametrize("seed", [31, 290])
def test_gram_svd_independent_certificates(seed):
    f, cases = load_exercise(9, seed)
    for (A,) in cases:
        before = A.copy()
        s, rebuilt, inverse = f(A)
        assert s.shape == (min(A.shape),) and rebuilt.shape == A.shape and inverse.shape == A.T.shape
        assert np.all(s >= 0) and np.all(np.diff(s) <= 1e-8)
        assert np.allclose(s, np.linalg.svd(A, compute_uv=False), rtol=1e-7, atol=1e-8)
        assert np.allclose(rebuilt, A, atol=1e-8)
        assert np.allclose(inverse, np.linalg.pinv(A), atol=1e-8)
        # 四个 Penrose 恒等式，同时核对亏秩分支。
        assert np.allclose(A@inverse@A, A)
        assert np.allclose(inverse@A@inverse, inverse)
        assert np.allclose(A@inverse, (A@inverse).T)
        assert np.allclose(inverse@A, (inverse@A).T)
        b = np.arange(A.shape[0]) + .7
        x = inverse @ b
        assert np.allclose(A.T@(b-A@x), 0, atol=1e-8)
        assert np.array_equal(A, before)
    # 固定起点与主方向正交、重奇异值、全零、秩一和矩形均需覆盖。
    for A in (np.diag([.4, 3., 1.]), 2*np.eye(3),
              np.array([[0.,0.,0.],[0.,3.,0.]]), np.zeros((3, 2))):
        s, rebuilt, inverse = f(A)
        assert np.allclose(rebuilt, A) and np.allclose(inverse, np.linalg.pinv(A))
        assert np.allclose(s, np.linalg.svd(A, compute_uv=False), atol=1e-8)
    with pytest.raises(ValueError):
        f(np.eye(4))


@pytest.mark.parametrize("seed", [12, 56])
def test_tensor_tools_independent_certificates(seed):
    f, cases = load_exercise(10, seed)
    for A, B, rho, (da, db), keep in cases:
        originals = [a.copy() for a in (A, B, rho)]
        product, reduced = f(A, B, rho, (da, db), keep)
        assert np.allclose(product, np.kron(A, B))
        shaped = rho.reshape(da, db, da, db)
        expected = np.trace(shaped, axis1=1, axis2=3) if keep == 0 else np.trace(shaped, axis1=0, axis2=2)
        assert np.allclose(reduced, expected)
        assert np.isclose(np.trace(reduced), 1) and np.allclose(reduced, reduced.conj().T)
        assert np.linalg.eigvalsh(reduced).min() >= -1e-12
        # 局域任意复厄米可观测量的期望值与联合空间计算一致。
        d = da if keep == 0 else db
        O = np.arange(d*d).reshape(d, d) + 1j*np.eye(d)
        O += O.conj().T
        joint = np.kron(O, np.eye(db)) if keep == 0 else np.kron(np.eye(da), O)
        assert np.allclose(np.trace(O@reduced), np.trace(joint@rho))
        assert all(np.array_equal(a, b) for a, b in zip((A, B, rho), originals))
    for keep in (0, 1):
        psi = np.array([np.sqrt(3), 0, 0, 1j])/2
        _, reduced = f(np.eye(2), np.eye(2), np.outer(psi, psi.conj()), (2, 2), keep)
        assert np.allclose(reduced, np.diag([.75, .25]))
        assert np.isclose(np.trace(reduced@reduced), .625)
    with pytest.raises(ValueError):
        f(np.eye(2), np.eye(2), np.eye(3), (2,2), 0)
    with pytest.raises(ValueError):
        f(np.eye(2), np.eye(2), np.eye(4), (2,2), 2)


def test_quantum_complex_schmidt_and_spin_examples():
    # 独立核对正文中 B 端复共轭约定；选非实右奇异基。
    C = np.array([[1, 1j], [2j, 1]], dtype=complex)/np.sqrt(7)
    U, s, Vh = np.linalg.svd(C)
    W = Vh.T
    psi = sum(s[j]*np.kron(U[:,j], W[:,j]) for j in range(2))
    assert np.allclose(psi, C.ravel())
    f, _ = load_exercise(10)
    _, rho_b = f(np.eye(2), np.eye(2), np.outer(psi, psi.conj()), (2,2), 1)
    assert np.allclose(rho_b, C.T@C.conj())
    assert np.allclose(np.linalg.eigvalsh(rho_b), np.sort(s*s))
    psi = np.array([np.sqrt(3), 1])/2
    plus = np.array([1,1])/np.sqrt(2)
    minus = np.array([1,-1])/np.sqrt(2)
    assert np.allclose([abs(plus@psi)**2, abs(minus@psi)**2], [(2+np.sqrt(3))/4, (2-np.sqrt(3))/4])
    for t in (0., .4, 1.7):
        psi = np.array([np.exp(-1j*t/2), 1j*np.exp(1j*t/2)])/np.sqrt(2)
        assert np.isclose(abs(plus@psi)**2, (1-np.sin(t))/2)


@pytest.mark.parametrize("function,args", [(viz.hermitian_spectrum, ([[2.,1.],[1.,-1.]],)),
                                            (viz.hermitian_spectrum, (np.eye(2),)),
                                            (viz.svd_ellipse, ([[1.,2.],[0.,0.]],)),
                                            (viz.tensor_dims, (2,3)), (viz.tensor_dims, (3,2)),
                                            (viz.tensor_dims, (1,1)), (viz.tensor_dims, (6,6))])
def test_round3_visualizations(function, args):
    figure = function(*args)
    try:
        assert isinstance(figure, Figure)
        figure.canvas.draw()
        if function == viz.tensor_dims:
            assert len(figure.axes[0].patches) == args[0]*args[1]
        if function == viz.svd_ellipse:
            assert len(figure.axes) == 3
            assert all(len(ax.texts) >= 2 for ax in figure.axes)
    finally:
        plt.close(figure)


def test_round3_flashcards_and_formula_manual():
    assert len(FLASHCARDS) >= 70
    assert {f"ch{n}" for n in range(1,11)} <= {c["chapter"] for c in FLASHCARDS}
    assert len({c["id"] for c in FLASHCARDS}) == len(FLASHCARDS)
    assert all(c["front"] and c["back"] for c in FLASHCARDS)
    for n in range(4, 11):
        theory = CHAPTERS[n]["theory_md"]
        for _, condition in FORMULA_SHEET[f"第{n}章"]:
            for chapter, tag in re.findall(r"式 \((\d+)\.(\d+)\)", condition):
                assert int(chapter) == n
                assert rf"\tag{{{chapter}.{tag}}}" in theory
    sheet8 = " ".join(f for f,_ in FORMULA_SHEET["第8章"])
    sheet9 = " ".join(f for f,_ in FORMULA_SHEET["第9章"])
    sheet10 = " ".join(f for f,_ in FORMULA_SHEET["第10章"])
    assert r"H=\sum_\lambda\lambda P_\lambda" in sheet8 and r"p(\lambda)=\|P_\lambda\psi\|^2" in sheet8
    assert r"A=U\Sigma V^\dagger" in sheet9 and r"A^+=V\Sigma^+U^\dagger" in sheet9
    assert r"\operatorname{Tr}_B\rho" in sheet10 and r"s_\alpha|u_\alpha\rangle\otimes|w_\alpha\rangle" in sheet10
    assert any("Schmidt" in c for _,c in FORMULA_SHEET["第10章"])
    assert "ch9" in KNOWLEDGE_MAP["ch10"]


def test_round3_diagnostic_mapping_and_recommendations():
    assert len(DIAGNOSTIC) == 14
    assert {q["chapter"] for q in DIAGNOSTIC} == set(KNOWLEDGE_MAP)
    assert [q["id"] for q in DIAGNOSTIC] == [f"diag-{j}" for j in range(1,15)]
    for q in DIAGNOSTIC:
        assert len(q["options"]) == 4 and 0 <= q["answer"] < 4
        if q["chapter"] in {"ch8", "ch9", "ch10"}:
            result = recommend({q["id"]:False})
            assert f"建议从第{q['chapter'][2:]}章" in result
            assert "建设中" not in result and "已答 1/14" in result
    by_chapter = {q["chapter"]:q["id"] for q in DIAGNOSTIC}
    assert "建议从第8章" in recommend({by_chapter[c]:False for c in ("ch8","ch9","ch10")})
    assert "建议从第6章" in recommend({by_chapter[c]:False for c in ("ch6","ch10")})
    assert "建议从地基" in recommend({by_chapter[c]:False for c in ("foundation","ch10")})
    assert "尚未完成" in recommend({"unknown":False})
    assert "已答 1/14" in recommend({by_chapter["ch10"]:True})
    assert "错题分布" not in recommend({q["id"]:True for q in DIAGNOSTIC})


@pytest.mark.parametrize("number", NUMBERS)
def test_round3_chapter_page_interaction(number, tmp_path, monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.setattr(progress, "DEFAULT_PATH", tmp_path / "progress.json")
    app = AppTest.from_file(str(Path(__file__).resolve().parents[1]/"app.py")).run(timeout=20)
    app.sidebar.radio[0].set_value(f"第{number}章").run(timeout=20)
    assert not app.exception
    assert app.title[0].value == f"第{number}章　{CHAPTER_TITLES[number]}"
    assert len(app.get("image")) == 1 and len(app.text_area) == 1 and len(app.radio) == 6
    assert any("## 核心推导" in x.value for x in app.markdown)
    if number == 8:
        app.slider(key="ch8-hermitian-b").set_value(1.2).run()
    elif number == 9:
        app.number_input(key="ch9-svd_ellipse-a₁₂").set_value(-1.4).run()
    else:
        app.slider(key="ch10-tensor-a").set_value(4).run()
        app.slider(key="ch10-tensor-b").set_value(2).run()
    assert not app.exception and len(app.get("image")) == 1
    q = CHAPTERS[number]["quiz"][0]
    app.radio(key=f"chapter-choice-{q['id']}").set_value(q["answer"])
    app.button(key=f"FormSubmitter:chapter-form-{q['id']}-提交答案").click().run()
    assert not app.exception
    state = progress.load(tmp_path / "progress.json")["chapters"][f"ch{number}"]
    assert state["quiz"] == 1 and q["id"] in state["quiz_correct_ids"]
