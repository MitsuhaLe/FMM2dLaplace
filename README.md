# FMM2dLaplace

---

## 求解的问题

光滑边界，2维Dirichlet边界条件的Laplace问题（区域内部）。

## 实现的方法

基于积分方程的偏微分方程快速算法，包含Nyström离散，FMM（快速多极子算法），奇异数值积分（用到Kapur-Rokhlin公式）

---

## 注意事项

fmm算法调用的是[fmmlib2d](https://github.com/zgimbutas/fmmlib2d).
使用main_interior_laplace以及运行文件夹“FMM”中的测试是需要将文件夹“matlab”加入到MATLAB路径
奇异数值积分部分参考了 [martinss/Nystrom/](https://amath.colorado.edu/faculty/martinss/Nystrom/) 中的代码。

## 使用方法
指定边界条件，类似于FMM/bdyfunc1的形式。
具体可以参考FMM/test_main_laplace_interior.
