"""全仓独立审读后的实质回归：绕过、合法等价算法与跨页持久化。"""
import ast
from datetime import date
from pathlib import Path
import re

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

from linalg_lab import progress
from linalg_lab.checks import banned_hits, check_code, _same
from linalg_lab.content import CHAPTERS
from linalg_lab.extras import DIAGNOSTIC, FORMULA_SHEET

pytestmark = pytest.mark.filterwarnings("ignore:Glyph .* missing from font")
EXERCISES = [e for c in CHAPTERS.values() for e in c['code_exercises']]
APP = str(Path(__file__).resolve().parents[1] / 'app.py')


@pytest.mark.parametrize('exercise', EXERCISES, ids=lambda e:e['id'])
@pytest.mark.parametrize('variant', ['alias', 'getattr', 'matrix', 'default', 'lambda', 'exec', 'eval', 'pickle'])
def test_all_chapters_reject_indirect_shortcuts(exercise, variant):
    # 注入到完整参考解前，防止仅因缺少目标函数而假称拦截成功。
    probes = {
        'alias': 'from numpy.linalg import inv as reciprocal\nf = reciprocal\n',
        'getattr': "f = getattr(np, 'lin' + 'alg')\n",
        'matrix': 'from numpy import asmatrix as wrap\nf = lambda A: wrap(A).I\n',
        'default': 'def hidden(A, f=np.linalg.svd):\n    return f(A)\n',
        'lambda': 'hidden = lambda A: np.linalg.inv(A)\n',
        'exec': "exec('f=np.' + 'linalg.inv')\n",
        'eval': "f = eval('np.' + 'linalg.inv')\n",
        'pickle': "def hidden():\n    return np.load('answer.npy', allow_pickle=True)\n",
    }
    source = probes[variant] + exercise['reference']
    assert banned_hits(source, exercise['banned'])
    passed, feedback = check_code(source, exercise)
    assert not passed and all('作弊检查失败' in line for line in feedback), feedback


def equivalent_source(exercise):
    """改变分配API、消元行遍历顺序；保持主元列/奇异值等有序契约。"""
    class Rewrite(ast.NodeTransformer):
        def visit_Call(self, node):
            self.generic_visit(node)
            if isinstance(node.func, ast.Attribute) and node.func.attr == 'zeros':
                node.func.attr = 'full'
                # zeros 未指定dtype时为float；full若用整数0会改变数值契约。
                node.keywords.append(ast.keyword(arg='fill_value', value=ast.Constant(0.0)))
            return node
        def visit_For(self, node):
            self.generic_visit(node)
            # 消元同列中各非主元行操作互不依赖；QR的j顺序保持。
            if (exercise['id'] not in {'ch8-c1','ch9-c1'}
                    and isinstance(node.target, ast.Name) and node.target.id == 'i'
                    and isinstance(node.iter, ast.Call) and isinstance(node.iter.func, ast.Name)
                    and node.iter.func.id == 'range' and len(node.iter.args) == 1):
                node.iter = ast.Call(func=ast.Name(id='reversed',ctx=ast.Load()), args=[node.iter], keywords=[])
            return node
    source = ast.unparse(ast.fix_missing_locations(Rewrite().visit(ast.parse(exercise['reference']))))
    if exercise['id'] == 'ch2-c1':
        # 独立换循环嵌套 i,j,k -> k,j,i：乘积累计顺序仍逐k。
        source = '''def matmul(A, B):
    m, n = A.shape
    p = B.shape[1]
    C = np.full((m,p), 0, dtype=np.result_type(A,B,float))
    for k in range(n):
        for j in reversed(range(p)):
            for i in range(m):
                C[i,j] += A[i,k]*B[k,j]
    return C
'''
    if exercise['id'] == 'ch3-c1':
        source = source.replace('R = np.full((m, n + 1), dtype=np.result_type(A, b, float), fill_value=0.0)',
                                'R = np.hstack((A, b.reshape(m, 1))).astype(np.result_type(A, b, float))')
        assert 'np.hstack' in source
    # 合法导入和别名不能被matrix/库写保护误杀。
    return 'import numpy as numeric\n' + source.replace('np.', 'numeric.')


@pytest.mark.parametrize('exercise', EXERCISES, ids=lambda e:e['id'])
def test_all_chapters_accept_equivalent_implementations(exercise):
    source = equivalent_source(exercise)
    assert not banned_hits(source, exercise['banned'])
    passed, feedback = check_code(source, exercise, seed=4719)
    assert passed, feedback


@pytest.mark.parametrize('api', ['inner', 'vecdot', 'vecmat', 'matvec', 'convolve', 'correlate'])
def test_multiplication_alternative_kernels_blocked(api):
    exercise = CHAPTERS[2]['code_exercises'][0]
    passed, feedback = check_code(f'def matmul(A,B):\n    return np.{api}(A,B.T)',exercise)
    assert not passed and '作弊检查失败' in feedback[0]


@pytest.mark.parametrize('prefix', [
    'np.allclose = lambda *args, **kwargs: True\n',
    'import numpy as lib\nlib.allclose = lambda *args, **kwargs: True\n',
    'lib = np\nlib.allclose = lambda *args, **kwargs: True\n',
    'def tamper(lib):\n    lib.allclose = lambda *args, **kwargs: True\ntamper(np)\n',
    'def tamper(lib):\n    lib._module.allclose = lambda *args, **kwargs: True\ntamper(np)\n',
    'def tamper(lib):\n    del lib.allclose\ntamper(np)\n',
    "np.__dict__['allclose'] = lambda *args, **kwargs: True\n",
    'import numpy as lib\ndef tamper(x):\n    x.allclose = lambda *a, **k: True\ntamper(lib)\n',
])
def test_wrong_solution_cannot_rewrite_verifier(prefix):
    exercise = CHAPTERS[2]['code_exercises'][0]
    source = prefix+'def matmul(A,B):\n    return np.zeros((A.shape[0], B.shape[1]))\n'
    passed, feedback = check_code(source,exercise)
    assert not passed and ('不允许' in feedback[0]), feedback
    assert not np.allclose([0],[1])  # 子进程之外也保持正常


def test_indices_require_integer_dtype():
    e = CHAPTERS[1]['code_exercises'][0]
    source = e['reference'].replace('np.array(pivots, dtype=int)','np.array(pivots, dtype=float)')
    passed, feedback = check_code(source,e)
    assert not passed and any('不符' in x for x in feedback)
    assert not _same(np.array([1.,2.]),np.array([1,2]))
    assert _same(np.array([1,2],dtype=np.int32),np.array([1,2],dtype=np.int64))


def test_numpy_imports_preserve_legal_public_apis():
    e = CHAPTERS[3]['code_exercises'][0]
    source = equivalent_source(e).replace('import numpy as numeric','from numpy import hstack\nimport numpy as numeric').replace('numeric.hstack','hstack')
    passed, feedback = check_code(source,e)
    assert passed, feedback


@pytest.mark.parametrize('source', [
    "np.lib.format.open_memmap('answer.npy')",
    "from numpy.lib.format import read_array\nread_array(f)",
    "np.fromregex('answer.txt', '.*', float)",
    "from numpy import DataSource\nDataSource().open('answer.txt')",
    "np.lib.npyio.DataSource().open('answer.txt')",
])
def test_file_loading_alternatives_blocked(source):
    e = EXERCISES[0]
    passed, feedback = check_code(source+'\n'+e['reference'], e)
    assert not passed and '作弊检查失败' in feedback[0]


@pytest.mark.parametrize('chapter',CHAPTERS.values(),ids=lambda c:c['id'])
def test_quiz_latex_commands_have_one_backslash(chapter):
    for q in chapter['quiz']:
        for text in [q['q_md'],q['explain_md'],*q['options']]:
            assert not re.search(r'(?<!\\)\\\\[A-Za-z|]',text),(q['id'],text)


def test_volume_formula_covers_singular_infinite_boundary():
    condition = next(condition for f,condition in FORMULA_SHEET['第4章'] if 'Vol' in f)
    assert '有限体积' in condition and '0·∞' in condition
    assert '0\\cdot\\infty' in CHAPTERS[4]['theory_md']


def fresh_app(monkeypatch,tmp_path):
    monkeypatch.setattr(progress,'DEFAULT_PATH',tmp_path/'progress.json')
    app=AppTest.from_file(APP).run(timeout=30)
    assert not app.exception
    return app


def navigate(app,page):
    app.sidebar.radio[0].set_value(page).run(timeout=30)
    assert not app.exception,(page,app.exception)


def test_every_page_visualization_redraw_and_persistence(tmp_path,monkeypatch):
    app=fresh_app(monkeypatch,tmp_path)
    pages=list(app.sidebar.radio[0].options)
    assert len(pages)==18
    for page in pages:
        navigate(app,page)
        if not page.startswith('第'):continue
        n=int(page[1:-1]);c=CHAPTERS[n]
        app.checkbox(key=f'theory-ch{n}').check().run()
        app.checkbox(key=f'chapter-done-ch{n}-w1').check().run()
        assert not app.exception
        # 无图的第1、3章明确跳过；其余各图修改真实控件后比较图像的内容地址。
        for viz in c['viz']:
            before=app.get('imgs')[0].proto.SerializeToString() if app.get('imgs') else app.get('image')[0].proto.SerializeToString()
            if viz in {'grid_transform','det_area','eigen_demo','svd_ellipse','matrix_power'}:
                app.number_input(key=f'ch{n}-{viz}-a₁₂').set_value(-1.4).run()
            elif viz=='gram_schmidt':app.slider(key='ch6-gs').set_value(2).run()
            elif viz=='hermitian_spectrum':app.slider(key='ch8-hermitian-b').set_value(1.4).run()
            elif viz=='tensor_dims':app.slider(key='ch10-tensor-a').set_value(4).run()
            assert not app.exception
            after=app.get('imgs')[0].proto.SerializeToString() if app.get('imgs') else app.get('image')[0].proto.SerializeToString()
            assert before!=after,(n,viz,'未重绘')
            if viz=='matrix_power':app.slider(key='ch7-power').set_value(8).run()
            if viz=='tensor_dims':app.slider(key='ch10-tensor-b').set_value(5).run()
            assert not app.exception
    # 完成记录经过所有18页切换、重新载入磁盘，且新会话可恢复。
    saved=progress.load(tmp_path/'progress.json')
    for n in CHAPTERS:
        assert saved['chapters'][f'ch{n}']['theory_done']
        assert f'ch{n}-w1' in saved['chapters'][f'ch{n}']['written']
        navigate(app,f'第{n}章')
        assert app.checkbox(key=f'theory-ch{n}').value
        assert app.checkbox(key=f'chapter-done-ch{n}-w1').value
    restarted=AppTest.from_file(APP).run(timeout=30)
    navigate(restarted,'第1章')
    assert restarted.checkbox(key='theory-ch1').value
    assert restarted.checkbox(key='chapter-done-ch1-w1').value


def test_quiz_single_question_navigation_all_answers_and_notebook(tmp_path,monkeypatch):
    app=fresh_app(monkeypatch,tmp_path)
    for n,c in CHAPTERS.items():
        navigate(app,f'第{n}章')
        for q in c['quiz']:
            app.selectbox(key=f'quiz-select-ch{n}').set_value(q['id']).run()
            assert len([x for x in app.radio if x.label=='选择答案'])==1
            app.radio(key=f"chapter-choice-{q['id']}").set_value(q['answer'])
            app.button(key=f"FormSubmitter:chapter-form-{q['id']}-提交答案").click().run()
            assert not app.exception and app.success
        assert progress.load(tmp_path/'progress.json')['chapters'][f'ch{n}']['quiz']==5
    # 真正从错题本重练，而非仅回原章答对。
    navigate(app,'第1章')
    app.selectbox(key='quiz-select-ch1').set_value('ch1-q1').run()
    app.radio(key='chapter-choice-ch1-q1').set_value(1)
    app.button(key='FormSubmitter:chapter-form-ch1-q1-提交答案').click().run()
    app.checkbox(key='chapter-mistake-ch1-w1').check().run()
    navigate(app,'错题本')
    app.radio(key='notebook-choice-ch1-q1').set_value(0)
    app.button(key='FormSubmitter:notebook-form-ch1-q1-提交答案').click().run()
    assert not app.exception
    assert progress.load(tmp_path/'progress.json')['chapters']['ch1']['quiz']==5
    app.checkbox(key='notebook-done-ch1-w1').check().run()
    app.button(key='remove-ch1-q1').click().run()
    app.button(key='remove-ch1-w1').click().run()
    assert not app.exception and not progress.load(tmp_path/'progress.json')['mistakes']
    navigate(app,'总览');navigate(app,'错题本')
    assert any('错题本为空' in x.value for x in app.info)


def test_formula_diagnostic_flashcard_ui(tmp_path,monkeypatch):
    app=fresh_app(monkeypatch,tmp_path)
    navigate(app,'公式手册')
    count=sum(map(len,FORMULA_SHEET.values()))
    assert sum(x.value.startswith('成立条件：') for x in app.markdown)==count
    app.text_input[0].set_value('极化').run()
    assert sum(x.value.startswith('成立条件：') for x in app.markdown)==1
    navigate(app,'入学诊断')
    app.button(key='FormSubmitter:diagnostic-form-提交诊断').click().run()
    assert not progress.load(tmp_path/'progress.json')['diagnostic']
    for q in DIAGNOSTIC:app.radio(key=f"diagnostic-{q['id']}").set_value(q['answer'])
    app.button(key='FormSubmitter:diagnostic-form-提交诊断').click().run()
    assert all(progress.load(tmp_path/'progress.json')['diagnostic'].values())
    navigate(app,'总览');navigate(app,'入学诊断')
    assert len(progress.load(tmp_path/'progress.json')['diagnostic'])==14
    navigate(app,'闪卡自测')
    from linalg_lab.extras import FLASHCARDS
    card=FLASHCARDS[0]
    assert not any(x.value==card['back'] for x in app.markdown)
    app.toggle(key=f"back-{card['id']}").set_value(True).run()
    assert any(x.value==card['back'] for x in app.markdown)
    app.checkbox(key=f"mastered-{card['id']}").check().run()
    navigate(app,'总览');navigate(app,'闪卡自测')
    assert app.checkbox(key=f"mastered-{card['id']}").value
    assert card['id'] in progress.load(tmp_path/'progress.json')['flashcards_mastered']
    assert not app.exception


def test_daily_question_follows_date_and_preserves_shared_progress(tmp_path,monkeypatch):
    import datetime
    class FixedDate(date):
        current=date(2026,10,25)
        @classmethod
        def today(cls):return cls.current
    monkeypatch.setattr(datetime,'date',FixedDate)
    app=fresh_app(monkeypatch,tmp_path)
    bank=[(c,w) for c in CHAPTERS.values() for w in c['written']]
    for day in [date(2026,10,25),date(2026,10,26),date(2026,10,25)]:
        FixedDate.current=day
        navigate(app,'每日一题')
        c,w=bank[day.toordinal()%len(bank)]
        assert any(day.isoformat() in x.value for x in app.caption)
        assert any(x.value==w['problem_md'] for x in app.markdown)
        app.checkbox(key=f"daily-{day.isoformat()}-done-{w['id']}").check().run()
        navigate(app,f"第{c['num']}章")
        assert app.checkbox(key=f"chapter-done-{w['id']}").value
        assert w['id'] in progress.load(tmp_path/'progress.json')['chapters'][c['id']]['written']
    assert not app.exception
