#ifndef MYVARIABLES_H
#define MYVARIABLES_H

#include <string>
#include <cmath>

const std::string boxtype = "cubic";
// std::string boxtype = "xyz";

// fcc xyz
const int nx = 4;
const int ny = 4;
const int nz = 12;

// fcc cubic
const int replicate = 2;

const int numParticles = (boxtype == "xyz") ? nx * ny * nz : 4 * replicate * replicate * replicate;

// 粒子直径作为特征长度 sigma=1.0
const double sigma = 1.0;

// 取相对密度为 rho / rho0 = 0.7360
const double relativeDensity = 0.7360;
const double latticeConstant = sqrt(2.0) * sigma / pow(relativeDensity, 1.0 / 3.0);

// 温度：k_B T = 1.0
const double temperature = 1.0;

extern int numSteps;  // 模拟步数
extern double lambda; // 弹簧系数

#endif