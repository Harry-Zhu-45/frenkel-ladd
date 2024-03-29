# Frenkel Ladd Free Energy

## 基本信息

该仓库是对 Daan Frenkel, Anthony J. C. Ladd 的文章 [New Monte Carlo method to compute the free energy of arbitrary solids. Application to the fcc and hcp phases of hard spheres](https://pubs.aip.org/aip/jcp/article-abstract/81/7/3188/91565/New-Monte-Carlo-method-to-compute-the-free-energy) 的部分内容复现。

除了 F-L 的文章，还参考了 [Understanding Molecular Simulation - From Algorithms to Applications](https://www.sciencedirect.com/book/9780323902922/understanding-molecular-simulation) 中对自由能的推导。

## 数据可视化 (optional)

简单修改源代码，可以获得蒙特卡洛模拟中 每一步 / 每次sample 的粒子位置数据文件 `.pos`。

软件包 `Injavis` (INteractive JAva VISualization 的缩写) 可以显示、分析和操纵粒子模拟数据。

所需的文件可以在 [https://engellab.de/teaching/injavis](https://engellab.de/teaching/injavis) 或 [https://zenodo.org/records/4639570](https://zenodo.org/records/4639570) 获得。

运行该软件包需要 `java` 运行环境：

```bash
java -Xmx4096m -jar injavis.jar
```

## 其他文件说明

可以使用 `minimum_distance.py` 判断两球形粒子之间的最近距离，判断是否存在 overlap。

可以使用 `analyze.py` 处理所得数据。
