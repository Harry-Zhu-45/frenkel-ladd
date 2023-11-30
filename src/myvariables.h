#ifndef MYVARIABLES_H
#define MYVARIABLES_H

#include <string>
#include <cmath>

const int nx = 3;
const int ny = 3;
const int nz = 6;
const int numParticles = nx * ny * nz;

// 粒子直径作为特征长度 sigma=1，所以半径取 0.5
const double particleRadius = 0.5;

// 密度定义为 rho = N * sigma**3 / V
// 取相对密度为 rho / rho0 = 0.7360
const double relativeDensity = 0.7360;

// a = a0 / 0.7360^{1/3} = sqrt(2) * sigma / 0.7360^{1/3} = 1.5663509 * sigma = 1.5663509
// const double latticeConstant = 1.5663509;
const double latticeConstant = sqrt(2.0) * 2 * particleRadius / pow(relativeDensity, 1.0 / 3.0);

// 温度：k_B T = 1.0
const double temperature = 1.0;

extern int numSteps;  // 模拟步数
extern double lambda; // 弹簧系数

#endif