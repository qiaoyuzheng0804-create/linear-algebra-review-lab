"""章节注册入口：新章导入并加入 CHAPTERS，页面无需另写。"""
from .ch01 import CHAPTER as CH01
from .ch02 import CHAPTER as CH02
from .ch03 import CHAPTER as CH03

CHAPTERS = {chapter["num"]: chapter for chapter in (CH01, CH02, CH03)}
CHAPTER_TITLES = {
    1: "向量空间、线性组合与线性相关",
    2: "矩阵与线性变换",
    3: "线性方程组与高斯消元",
    4: "行列式与体积",
    5: "特征值、特征向量与对角化",
    6: "内积、正交化与最小二乘",
    7: "奇异值分解与低秩近似",
    8: "厄米矩阵、酉矩阵与谱定理",
    9: "矩阵函数与线性动力学",
    10: "张量积与复合量子系统",
}
