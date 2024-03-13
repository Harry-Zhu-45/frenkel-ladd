# README

## 基本信息

该仓库是对 Daan Frenkel; Anthony J. C. Ladd 的文章 [New Monte Carlo method to compute the free energy of arbitrary solids. Application to the fcc and hcp phases of hard spheres](https://pubs.aip.org/aip/jcp/article-abstract/81/7/3188/91565/New-Monte-Carlo-method-to-compute-the-free-energy) 的部分内容复现，算法细节上比较粗糙，有很大的优化空间。

## 编译与使用

`nx, ny, nz` 等变量可以在 `src/myvariables.h` 中设置。

设置好需要模拟的参数后，在根目录下使用 `cmake` 编译

```bash
cmake -B build

cmake --build build
```

编译后的可执行文件位于 `./build/main.out`

可以用 `-l` 指定 `lambda`，用 `-n` 指定 Monte Carlo 模拟步数。例如：

```bash
./build/main.out -l 10 -n 10000
```

## 数据可视化 (optional)

除了 msd 文件，还可以在 `src/main.cpp` 中简单修改，获得 MC 模拟中 sample 的粒子位置数据文件 `.pos`。

一般来说，`.pos` 文件由以下部分组成

盒子可以用长宽高描述：

```text
box Lx Ly Lz
```

或者使用 3x3 的矩阵描述：

```text
boxMatrix Lx 0 0 0 Ly 0 0 0 Lz
```

每一帧内，可以描述粒子的形状、颜色、位置。(如果有朝向的话，可以再末尾加上四元数表示朝向。)

```text
particleType color position_x position_y position_z (quaternion_1 quaternion_2 quaternion_3 quaternion_4)
```

每一帧的结尾为

```text
eof
```

> 重复上述数据部分，可以记录下每一帧的粒子位置。

软件包 Injavis (INteractive JAva VISualization 的缩写) 可以显示、分析和操纵粒子模拟数据。所需的文件可以在 [https://engellab.de/teaching/injavis](https://engellab.de/teaching/injavis) 或 [https://zenodo.org/records/4639570](https://zenodo.org/records/4639570) 获得。

运行该软件包需要 java 运行环境：

```bash
java -Xmx4096m -jar injavis.jar
```

## 其他文件说明

可以使用 `minimum_distance.py` 判断两粒子之间的最近距离。

可以使用 `collect.sh` 连续多次模拟。其中 `values` 为待模拟的 `lambda` 的取值。模拟获得的 msd 结果将默认保存在根目录下的文件 `msd_data.txt` 中。需注意，重复运行脚本将会覆盖之前的结果。

可以使用 `10p-GLI.py` 对所得的 msd 进行 G-L 积分。
