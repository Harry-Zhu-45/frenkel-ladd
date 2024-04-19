# Compute Frenkel Ladd Free Energy

## 基本信息

该仓库是对以下两篇文章的部分内容复现。

1. Daan Frenkel, Anthony J. C. Ladd 的文章 [New Monte Carlo method to compute the free energy of arbitrary solids. Application to the fcc and hcp phases of hard spheres](https://pubs.aip.org/aip/jcp/article-abstract/81/7/3188/91565/New-Monte-Carlo-method-to-compute-the-free-energy)
2. Carlos Vega, Eva G. Noya 的文章 [Revisiting the Frenkel-Ladd method to compute the free energy of solids: The Einstein molecule approach](https://pubs.aip.org/aip/jcp/article/127/15/154113/914715/Revisiting-the-Frenkel-Ladd-method-to-compute-the)

除了以上两篇文章，还参考了 [Understanding Molecular Simulation - From Algorithms to Applications](https://www.sciencedirect.com/book/9780323902922/understanding-molecular-simulation) 中对固定质心体系自由能的推导。

## 数据可视化 (optional)

简单修改源代码，可以获得蒙特卡洛模拟中的粒子位置数据文件 `.pos`。

`Injavis` (INteractive JAva VISualization 的缩写) 可以显示、分析和操纵粒子模拟数据。

所需的文件可以在 [https://engellab.de/teaching/injavis](https://engellab.de/teaching/injavis) 或 [https://zenodo.org/records/10125525](https://zenodo.org/records/10125525) 获得。

运行该软件包需要 up-to-date 的 `java` 运行时环境 (JRE)：

```bash
java -Xmx4096m -jar injavis.jar
```
