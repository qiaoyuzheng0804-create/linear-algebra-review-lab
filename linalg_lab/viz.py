"""纯 Matplotlib 图形；调用方负责展示与关闭 Figure。"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib import font_manager

# 按已安装字体选择，不要求外部字体文件。
_AVAILABLE = {font.name for font in font_manager.fontManager.ttflist}
for _name in ("Noto Sans CJK SC", "Noto Sans CJK JP", "WenQuanYi Zen Hei", "SimHei"):
    if _name in _AVAILABLE:
        plt.rcParams["font.sans-serif"] = [_name, "DejaVu Sans"]
        break
plt.rcParams["axes.unicode_minus"] = False


def _axes(ax, points=None):
    extent = max(2.0, float(np.max(np.abs(points))) * 1.2) if points is not None else 3.0
    ax.set(xlim=(-extent, extent), ylim=(-extent, extent), aspect="equal")
    ax.axhline(0, color="gray", lw=.6)
    ax.axvline(0, color="gray", lw=.6)
    ax.grid(alpha=.15)


def _arrow(ax, vector, color, label):
    ax.annotate("", xy=vector, xytext=(0, 0), arrowprops={"arrowstyle": "->", "color": color, "lw": 2})
    ax.text(*vector, label, color=color)


def _matrix(A):
    A = np.asarray(A, dtype=float)
    if A.shape != (2, 2) or not np.isfinite(A).all():
        raise ValueError("请输入有限实数构成的 2×2 矩阵")
    return A


def grid_transform(A):
    A = _matrix(A)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    extent = np.c_[A @ np.array([[x, y] for x in (-2, 2) for y in (-2, 2)]).T, np.eye(2)]
    for ax, matrix, title in zip(axes, (np.eye(2), A), ("原网格", "线性变换后的网格")):
        for value in np.linspace(-2, 2, 9):
            for line in (np.array([[-2, 2], [value, value]]), np.array([[value, value], [-2, 2]])):
                mapped = matrix @ line
                ax.plot(*mapped, color="steelblue", alpha=.5)
        _arrow(ax, matrix[:, 0], "darkorange", "$e_1$" if title == "原网格" else "$Ae_1$")
        _arrow(ax, matrix[:, 1], "purple", "$e_2$" if title == "原网格" else "$Ae_2$")
        _axes(ax, extent)
        ax.set_title(title)
    fig.suptitle("看什么：基向量的像决定整个网格的像，平行关系被保留。")
    fig.tight_layout()
    return fig


def det_area(A):
    A = _matrix(A)
    square = np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
    mapped = square @ A.T
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.add_patch(Polygon(square, alpha=.25, color="gray", label="单位正方形，面积 1"))
    ax.add_patch(Polygon(mapped, alpha=.4, color="steelblue", label=f"像的面积 {abs(np.linalg.det(A)):.3g}"))
    _axes(ax, np.r_[square, mapped])
    ax.legend()
    ax.set_title(f"看什么：面积比为 |det A| = {abs(np.linalg.det(A)):.3g}，符号记录定向。")
    fig.tight_layout()
    return fig


def eigen_demo(A):
    A = _matrix(A)
    circle = np.array([np.cos(np.linspace(0, 2*np.pi, 240)), np.sin(np.linspace(0, 2*np.pi, 240))])
    mapped = A @ circle
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(*circle, label="单位圆")
    ax.plot(*mapped, label="单位圆的像")
    values, vectors = np.linalg.eig(A)
    count = 0
    for j, value in enumerate(values):
        if abs(value.imag) < 1e-9:
            direction = vectors[:, j].real
            line = np.outer(direction, [-2, 2])
            ax.plot(*line, "--", label=f"特征方向，λ={value.real:.2g}")
            count += 1
    if not count:
        ax.text(.02, .03, "此矩阵无实特征方向", transform=ax.transAxes)
    _axes(ax, np.c_[circle, mapped])
    ax.legend()
    ax.set_title("看什么：实特征方向上的向量只伸缩，不离开原直线。")
    fig.tight_layout()
    return fig


def gram_schmidt(vectors, step):
    vectors = np.asarray(vectors, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 2 or not len(vectors):
        raise ValueError("每行给出一个二维向量")
    step = max(1, min(int(step), len(vectors)))
    fig, ax = plt.subplots(figsize=(6, 5))
    orthogonal = []
    for i, vector in enumerate(vectors[:step]):
        residual = vector.copy()
        _arrow(ax, vector, "gray", f"$v_{i+1}$")
        for q in orthogonal:
            residual -= np.dot(residual, q) / np.dot(q, q) * q
        if np.dot(residual, residual) > 1e-20:
            orthogonal.append(residual)
            _arrow(ax, residual, f"C{i % 10}", f"$u_{i+1}$")
        else:
            ax.text(.02, .03 + .06*i, f"第 {i+1} 个残差为零：向量相关", transform=ax.transAxes)
    _axes(ax, vectors)
    ax.set_title(f"看什么：第 {step} 步扣除已有方向上的投影，留下垂直分量。")
    fig.tight_layout()
    return fig


def svd_ellipse(A):
    A = _matrix(A)
    U, singular, Vt = np.linalg.svd(A)
    t = np.linspace(0, 2*np.pi, 240)
    circle = np.array([np.cos(t), np.sin(t)])
    stages = [Vt @ circle, np.diag(singular) @ Vt @ circle, A @ circle]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, points, title in zip(axes, stages, ("① $V^T$ 正交变换", "② $Σ$ 沿坐标轴伸缩", "③ $U$ 正交变换")):
        ax.plot(*circle, "--", color="gray", label="原单位圆")
        ax.plot(*points, color="steelblue", label="当前像")
        _axes(ax, np.c_[circle, *stages])
        ax.set_title(title)
        ax.legend(fontsize=8)
    fig.suptitle("看什么：圆先旋转或反射，再伸缩，最后旋转或反射，得到椭圆。")
    fig.tight_layout()
    return fig


def matrix_power(A, k):
    A = _matrix(A)
    k = int(k)
    if not 0 <= k <= 100:
        raise ValueError("迭代步数应在 0～100 之间")
    points = [np.array([1., .5])]
    for _ in range(k):
        points.append(A @ points[-1])
    points = np.asarray(points)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(*points.T, "o-", markersize=4)
    for i, point in enumerate(points):
        ax.annotate(str(i), point)
    _axes(ax, points)
    ax.set_title("看什么：重复作用 A 后，轨迹可能收缩、增长或沿特征方向排列。")
    fig.tight_layout()
    return fig


def tensor_dims():
    fig, ax = plt.subplots(figsize=(7, 4))
    for i in range(2):
        for j in range(3):
            ax.add_patch(Rectangle((j, 1-i), 1, 1, fill=False, edgecolor="steelblue", lw=2))
            ax.text(j+.5, 1-i+.5, f"$e_{i+1} \\otimes f_{j+1}$", ha="center", va="center")
    ax.set(xlim=(-.2, 3.2), ylim=(-.2, 2.2), aspect="equal")
    ax.axis("off")
    ax.set_title("看什么：2 个基向量与 3 个基向量两两配对，张量积维数为 6。")
    fig.tight_layout()
    return fig
