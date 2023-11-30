#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <random>

#include "MonteCarloSimulation.h"
#include "myvariables.h"

int numSteps = 10000 * numParticles; // 模拟步数
double lambda = 7.0;                 // 弹簧系数

MonteCarloSimulation::MonteCarloSimulation(int numParticles)
{
    // fcc initialization
    lattice.reserve(numParticles);   // 预分配内存
    particles.reserve(numParticles); // 预分配内存

    double a1 = latticeConstant;
    double a0 = latticeConstant / sqrt(2.0);

    int i = 0; // 一开始没有粒子
    Vector3D centerOfMass = Vector3D(0, 0, 0);

    for (int iz = 0; iz < nz; ++iz)
    {
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                i = i + 1;
                lattice.emplace_back(Vector3D(a0 * ix + (a0 / 2.0) * fmod(iz, 2),
                                              a0 * iy + (a0 / 2.0) * fmod(iz, 2),
                                              (a1 / 2.0) * iz));
                particles.emplace_back(Vector3D(a0 * ix + (a0 / 2.0) * fmod(iz, 2),
                                                a0 * iy + (a0 / 2.0) * fmod(iz, 2),
                                                (a1 / 2.0) * iz));

                centerOfMass = centerOfMass + lattice[i - 1].position; // 计算质心
            }
        }
    }

    // 固定质心于原点
    centerOfMass = centerOfMass * (1 / numParticles);
    for (int i = 0; i < numParticles; ++i)
    {
        lattice[i].position = lattice[i].position - centerOfMass;
        particles[i].position = particles[i].position - centerOfMass;
    }
}

double MonteCarloSimulation::calculateSD()
{
    double totalSD = 0.0;
    for (int i = 0; i < numParticles; ++i)
    {
        // MSD = (r-r0)^2
        totalSD += (particles[i].position - lattice[i].position).magnitude() * (particles[i].position - lattice[i].position).magnitude();
    }
    return totalSD;
}

void MonteCarloSimulation::metropolisStep()
{
    // 创建随机数生成器对象
    std::random_device rd;                                             // 用于种子生成
    std::mt19937 gen(rd());                                            // Mersenne Twister 伪随机数生成器
    std::uniform_real_distribution<double> randomDouble(0.0, 1.0);     // 生成(0, 1)之间的均匀分布随机数
    std::uniform_int_distribution<int> randomInt(0, numParticles - 1); // 生成(0, 1)之间的均匀分布随机数

    int randomIndex = randomInt(gen);                                                                 // 随机选取一个粒子
    Vector3D displacement(randomDouble(gen) - 0.5, randomDouble(gen) - 0.5, randomDouble(gen) - 0.5); // 产生随机位移
    displacement = displacement * 0.05;                                                               // 调整位移幅度
    Vector3D virtualDisplacement = particles[randomIndex].position + displacement;                    // 虚拟位移

    // 判断是否存在重叠
    for (int i = 0; i < numParticles; ++i)
    {
        if (i != randomIndex) // 遍历除选中粒子外的所有粒子
        {
            Vector3D delta = particles[i].position - virtualDisplacement;
            Vector3D displacement[27];

            // 27 种位移，相当于周期性边界条件
            for (int dx = -1; dx <= 1; dx++)
            {
                for (int dy = -1; dy <= 1; dy++)
                {
                    for (int dz = -1; dz <= 1; dz++)
                    {
                        displacement[(dx + 1) * 9 + (dy + 1) * 3 + (dz + 1)] = delta + Vector3D(dx * nx * latticeConstant / sqrt(2.0), dy * ny * latticeConstant / sqrt(2.0), dz * nz * latticeConstant / 2);
                    }
                }
            }

            for (int j = 0; j < 27; j++)
            {
                if (displacement[j].magnitude() < 2 * particleRadius)
                {
                    return; // 结束当前的 metropolisStep
                }
            }
        }
    }

    // 不重叠，计算可能的势能变化
    // deltaEnergy = lambda (2 delta r_i dot delta_i + (N-1)/N delta_i^2)
    double deltaEnergy = lambda * 2 * (particles[randomIndex].position - lattice[randomIndex].position).dot(displacement) + lambda * (numParticles - 1) / numParticles * displacement.magnitude() * displacement.magnitude();

    if (randomDouble(gen) < exp(-deltaEnergy / temperature)) // Metropolis 准则
    {
        // 接受位移
        particles[randomIndex].position = virtualDisplacement;
        for (int i = 0; i < numParticles; ++i)
        {
            particles[i].position = particles[i].position - displacement * (1 / numParticles);
        }
    }
}