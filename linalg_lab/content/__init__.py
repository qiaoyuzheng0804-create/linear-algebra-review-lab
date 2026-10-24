"""章节注册入口：新章导入并加入 CHAPTERS，页面无需另写。"""
from .ch01 import CHAPTER as CH01
from .ch02 import CHAPTER as CH02
from .ch03 import CHAPTER as CH03
from .ch04 import CHAPTER as CH04
from .ch05 import CHAPTER as CH05
from .ch06 import CHAPTER as CH06
from .ch07 import CHAPTER as CH07
from .ch08 import CHAPTER as CH08
from .ch09 import CHAPTER as CH09
from .ch10 import CHAPTER as CH10

CHAPTERS = {chapter["num"]: chapter for chapter in (CH01, CH02, CH03, CH04, CH05, CH06, CH07, CH08, CH09, CH10)}
CHAPTER_TITLES = {
    1: "向量空间、线性组合与线性相关",
    2: "矩阵与线性变换",
    3: "线性方程组与高斯消元",
    4: "行列式",
    5: "特征值与特征向量",
    6: "内积空间、正交化与 Gram–Schmidt",
    7: "对角化与矩阵指数",
    8: "厄米矩阵、酉矩阵与谱定理",
    9: "奇异值分解与低秩近似",
    10: "张量积与复合量子系统",
}
