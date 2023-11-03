#ifndef MYVARIABLES_H
#define MYVARIABLES_H

#include <string>
const std::string latticeFileName = "fcc-3x4x6.txt"; // lattice 文件名

// 粒子直径作为特征长度 sigma=1，所以半径取 0.5
const double particleRadius = 0.5;

// 密度定义为 rho = N * sigma**3 / V
// 取相对密度为 rho / rho0 = 0.7360
// a = a0 / 0.7360^{1/3} = sqrt(2) * sigma / 0.7360^{1/3} = 1.5663509 * sigma = 1.5663509
const double latticeConstant = 1.5663509;
const double temperature = 1.0; // 温度：k_B T = 1

extern int numSteps;     // 模拟步数
extern int numParticles; // 粒子数目
extern double lambda;    // 弹簧系数

#endif