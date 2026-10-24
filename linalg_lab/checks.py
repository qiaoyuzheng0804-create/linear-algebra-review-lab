"""AST 反直接调库检查与独立进程随机验证；不在主应用进程执行学生代码。"""
from __future__ import annotations

import ast
import multiprocessing as mp

import numpy as np


def _canonical(name):
    return "np" + name[5:] if name == "numpy" or name.startswith("numpy.") else name


def banned_hits(source, banned):
    tree = ast.parse(source)
    blocked = {_canonical(x) for x in banned}
    # NumPy 可能在其他模块重导出 linalg 的函数。禁用整个 linalg 时，
    # 同时保留其公开 API 名，与下面的属性/导入符号检查使用同一策略。
    if "np.linalg" in blocked:
        blocked.update(f"np.linalg.{name}" for name, value in vars(np.linalg).items()
                       if not name.startswith("_") and callable(value))
    reserved = {item.rsplit(".", 1)[-1] for item in blocked if "." in item}
    hits = set()
    scope_types = (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)

    def forbidden(name):
        return any(name == item or name.startswith(item + ".") for item in blocked)

    def local_nodes(scope):
        # 嵌套作用域的绑定不可反向覆盖当前作用域。
        result = []
        def descend(node):
            for child in ast.iter_child_nodes(node):
                result.append(child)
                if not isinstance(child, scope_types):
                    descend(child)
        descend(scope)
        return result

    def inspect_scope(scope, inherited):
        aliases = {name: set(paths) for name, paths in inherited.items()}
        nodes = local_nodes(scope)

        def resolve(node):
            if isinstance(node, ast.Name):
                return aliases.get(node.id, {node.id})
            if isinstance(node, ast.Attribute):
                return {f"{base}.{node.attr}" for base in resolve(node.value) if base.count(".") < 16}
            return set()

        def bind(name, paths):
            # 保留同一作用域内可能出现的每个绑定，防止后一次赋值掩盖早先调用。
            aliases.setdefault(name, set()).update(paths)

        for _ in range(len(nodes) + 1):
            previous = {name: set(paths) for name, paths in aliases.items()}
            for node in nodes:
                if isinstance(node, ast.Import):
                    for item in node.names:
                        full = _canonical(item.name)
                        bind(item.asname or item.name.split(".")[0], {full if item.asname else full.split(".")[0]})
                        if forbidden(full):
                            hits.add(full)
                        if not (item.name == "numpy" or item.name.startswith("numpy.")):
                            hits.add(f"不允许导入：{item.name}")
                elif isinstance(node, ast.ImportFrom):
                    module = _canonical(node.module or "")
                    for item in node.names:
                        full = f"{module}.{item.name}"
                        bind(item.asname or item.name, {full})
                        # 星号导入可隐藏重导出的禁用函数，教学编辑器要求显式导入。
                        if forbidden(full) or item.name in reserved or item.name == "*":
                            hits.add(full)
                    if not (module == "np" or module.startswith("np.")):
                        hits.add(f"不允许导入：{node.module}")
                elif isinstance(node, (ast.Assign, ast.AnnAssign, ast.NamedExpr)):
                    paths = resolve(node.value)
                    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                    for target in targets:
                        if isinstance(target, ast.Name):
                            bind(target.id, paths)
            if previous == aliases:
                break
        for node in nodes:
            if isinstance(node, ast.MatMult) and "@" in blocked:
                hits.add("@（矩阵乘法运算符）")
            if isinstance(node, ast.Attribute):
                # 禁用 API 的末级属性名在本练习中保留；即使库经参数、
                # global 或解包传递，也不能借改名直接调用该属性。
                if node.attr in reserved:
                    hits.update(f"{item}（命中属性 {node.attr}）" for item in blocked
                                if item.endswith("." + node.attr))
            if isinstance(node, (ast.Attribute, ast.Name)):
                hits.update(name for name in resolve(node) if forbidden(name))
            if isinstance(node, ast.Attribute) and node.attr.startswith("__"):
                hits.add(f"不允许访问内部属性：{node.attr}")
            if isinstance(node, ast.Name) and node.id in {"eval", "exec", "compile", "__import__", "getattr", "setattr", "globals", "locals", "open"}:
                hits.add(f"不允许动态执行或访问：{node.id}")
            if isinstance(node, scope_types):
                inspect_scope(node, aliases)
    inspect_scope(tree, {"np": {"np"}, "numpy": {"np"}})
    return sorted(hits)


def _same(actual, expected):
    if isinstance(expected, tuple):
        return (isinstance(actual, tuple) and len(actual) == len(expected)
                and all(_same(a, b) for a, b in zip(actual, expected)))
    if expected is None:
        return actual is None
    try:
        a, b = np.asarray(actual), np.asarray(expected)
        return a.shape == b.shape and bool(np.allclose(a, b, rtol=1e-7, atol=1e-8))
    except (TypeError, ValueError):
        return False


def _numpy_import(name, globals=None, locals=None, fromlist=(), level=0):
    if level or not (name == "numpy" or name.startswith("numpy.")):
        raise ImportError("仅允许导入 numpy")
    return __import__(name, globals, locals, fromlist, level)


def _worker(connection, source, exercise, seed):
    try:
        trusted = {"np": np}
        exec(exercise["reference"], trusted)
        exec(exercise["make_cases"], trusted)
        cases = list(trusted["make_cases"](np.random.default_rng(seed)))
        if len(cases) < 20:
            connection.send((False, ["练习配置错误：随机用例不足 20 组。 "]))
            return
        # 先保存预期值，避免学生函数修改输入影响参考输出。
        expected = [trusted[exercise["function_name"]](*args) for args in cases]
        import builtins
        allowed = {name: getattr(builtins, name) for name in (
            "abs", "all", "any", "bool", "complex", "dict", "enumerate", "Exception",
            "float", "int", "isinstance", "len", "list", "max", "min", "range",
            "reversed", "round", "set", "slice", "sorted", "sum", "tuple", "ValueError", "zip"
        )}
        allowed["__import__"] = _numpy_import
        namespace = {"np": np, "__builtins__": allowed}
        exec(compile(source, "<学生代码>", "exec"), namespace)
        function = namespace.get(exercise["function_name"])
        if not callable(function):
            connection.send((False, [f"未找到目标函数 {exercise['function_name']}。 "]))
            return
        feedback = []
        passed = True
        for index, (args, target) in enumerate(zip(cases, expected), 1):
            try:
                actual = function(*(np.array(x, copy=True) if isinstance(x, np.ndarray) else x for x in args))
                ok = _same(actual, target)
                feedback.append(f"用例 {index}：{'通过' if ok else '输出数值、形状或元组结构不符'}。")
                passed = passed and ok
            except Exception as exc:
                passed = False
                feedback.append(f"用例 {index}：执行异常 {type(exc).__name__}：{exc}。")
        connection.send((passed, feedback))
    except BaseException as exc:
        connection.send((False, [f"代码执行失败：{type(exc).__name__}：{exc}。 "]))
    finally:
        connection.close()


def check_code(source, exercise, *, seed=20261024, timeout=10):
    """返回 (通过与否, 中文反馈列表)。AST 是教学检查，非恶意代码安全沙箱。"""
    try:
        hits = banned_hits(source, exercise.get("banned", []))
    except SyntaxError as exc:
        return False, [f"语法错误：第 {exc.lineno} 行，{exc.msg}。"]
    if hits:
        return False, [f"作弊检查失败，命中：{item}。" for item in hits]
    context = mp.get_context("spawn")
    receiver, sender = context.Pipe(duplex=False)
    process = context.Process(target=_worker, args=(sender, source, exercise, seed), daemon=True)
    try:
        process.start()
        sender.close()
        if receiver.poll(timeout):
            try:
                return receiver.recv()
            except EOFError:
                return False, ["执行进程提前退出，请检查代码。"]
        return False, [f"执行超过 {timeout} 秒，请检查无限循环或算法复杂度。"]
    finally:
        if process.pid is not None:
            if process.is_alive():
                process.terminate()
            process.join()
        receiver.close()
        sender.close()
