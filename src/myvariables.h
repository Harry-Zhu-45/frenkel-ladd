#ifndef MYVARIABLES_H
#define MYVARIABLES_H

const std::string latticeFileName = "fcc-3x4x6.txt"; // lattice 文件名

const double latticeConstant = 2.0; // 晶格常数
const double particleRadius = 0.5;  // 粒子半径
const double temperature = 1.0;     // 温度

extern int numSteps;     // 模拟步数
extern int numParticles; // 粒子数目
extern double lambda;    // 弹簧系数

#endif