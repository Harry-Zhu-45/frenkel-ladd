# README

## injavis 可视化

打开 injavis 可视化界面

```bash
java -Xmx4096m -jar injavis.jar
```

`.pos` 文件为 injavis 专用数据文件

### pos 文件格式

长方体盒子，一帧开头为

```text
box Lx Ly Lz
```

或者使用 3x3 的矩阵描述

```text
boxMatrix Lx 0 0 0 Ly 0 0 0 Lz
```

每一帧内，可以描述粒子的 形状 颜色 位置

如果有朝向的话，可以再末尾加上四元数表示朝向

```text
particleType color position_x position_y position_z (quaternion_1 quaternion_2 quaternion_3 quaternion_4)
```

一帧结尾为

```text
eof
```

重复上述过程，可以实现动画效果

## 注意

更改 lattice 之后

- 更改 lattice
- 更改 particles
- 更改 number of particles
- 更改 pos 文件显示区域
