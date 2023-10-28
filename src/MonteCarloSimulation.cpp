#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "MonteCarloSimulation.h"
#include "randomDouble.h"
#include "myvariables.h"

int numParticles;                    // 粒子数目
int numSteps = 10000 * numParticles; // 模拟步数
double lambda = 7.0;                 // 弹簧系数

MonteCarloSimulation::MonteCarloSimulation(int numParticles)
{
    std::ifstream inputFile("src/lattice/" + latticeFileName);

    // 检查文件是否打开成功
    if (!inputFile.is_open())
    {
        std::cerr << "Failed to open lattice file." << std::endl;
        return;
    }

    // 读取 lattice 所对应的 numParticles
    inputFile >> numParticles;

    lattice.reserve(numParticles);   // 预分配内存
    particles.reserve(numParticles); // 预分配内存

    for (int i = 0; i < numParticles; ++i)
    {
        double x, y, z;
        inputFile >> x >> y >> z;
        lattice.emplace_back(Vector3D(x, y, z) * latticeConstant);
        particles.emplace_back(Vector3D(x, y, z) * latticeConstant);
    }

    inputFile.close();

    // 整体平移，使质心为 000
    for (int i = 0; i < numParticles; ++i)
    {
        lattice[i].position = lattice[i].position - Vector3D(0.75, 1.25, 1.25) * latticeConstant;
        particles[i].position = particles[i].position - Vector3D(0.75, 1.25, 1.25) * latticeConstant;
    }
}

double MonteCarloSimulation::calculateTotalEnergy()
{
    double totalEnergy = 0.0;
    for (int i = 0; i < numParticles; ++i)
    {
        // 无重叠，计算累加势能 V = lambda * (r-r0)^2
        totalEnergy += lambda * (particles[i].position - lattice[i].position).magnitude() * (particles[i].position - lattice[i].position).magnitude();
    }
    return totalEnergy;
}

void MonteCarloSimulation::metropolisStep()
{
    double energyBefore = calculateTotalEnergy(); // 变化前的势能

    int randomIndex = rand() % numParticles;             // 随机选取一个粒子
    Particle &selectedParticle = particles[randomIndex]; // 引用，方便修改

    // 对选中的粒子进行随机位移
    Vector3D displacement(randomDouble() - 0.5, randomDouble() - 0.5, randomDouble() - 0.5); // 产生随机位移
    displacement = displacement * 0.1;                                                       // 调整位移幅度
    selectedParticle.position = selectedParticle.position + displacement;                    // 位移

    // 判断是否重叠
    // 同时判断选中粒子 x y z 正负方向虚拟位移一个晶胞距离后是否有重叠
    // 一共判断 7 次
    for (int i = 0; i < numParticles && i != randomIndex; ++i) // 遍历除选中粒子外的所有粒子
    {
        Vector3D delta = particles[i].position - particles[randomIndex].position; // 两粒子间的位移
        double distance = delta.magnitude();                                      // 两粒子间的距离
        if (distance < 2 * particleRadius)
        {
            // 重叠：拒绝位移，回退到原位置
            selectedParticle.position = selectedParticle.position - displacement;
            return; // 结束当前的 metropolisStep
        }

        delta = particles[i].position - particles[randomIndex].position + Vector3D(2 * latticeConstant, 0.0, 0.0); // 两粒子间的位移 （x 正方向）
        distance = delta.magnitude();                                                                              // 两粒子间的距离
        if (distance < 2 * particleRadius)
        {
            // 重叠：拒绝位移，回退到原位置
            selectedParticle.position = selectedParticle.position - displacement;
            return; // 结束当前的 metropolisStep
        }

        delta = particles[i].position - particles[randomIndex].position + Vector3D(-2 * latticeConstant, 0.0, 0.0); // 两粒子间的位移 （x 负方向）
        distance = delta.magnitude();                                                                               // 两粒子间的距离
        if (distance < 2 * particleRadius)
        {
            // 重叠：拒绝位移，回退到原位置
            selectedParticle.position = selectedParticle.position - displacement;
            return; // 结束当前的 metropolisStep
        }

        delta = particles[i].position - particles[randomIndex].position + Vector3D(0.0, 3 * latticeConstant, 0.0); // 两粒子间的位移 （y 正方向）
        distance = delta.magnitude();                                                                              // 两粒子间的距离
        if (distance < 2 * particleRadius)
        {
            // 重叠：拒绝位移，回退到原位置
            selectedParticle.position = selectedParticle.position - displacement;
            return; // 结束当前的 metropolisStep
        }

        delta = particles[i].position - particles[randomIndex].position + Vector3D(0.0, -3 * latticeConstant, 0.0); // 两粒子间的位移 （y 负方向）
        distance = delta.magnitude();                                                                               // 两粒子间的距离
        if (distance < 2 * particleRadius)
        {
            // 重叠：拒绝位移，回退到原位置
            selectedParticle.position = selectedParticle.position - displacement;
            return; // 结束当前的 metropolisStep
        }

        delta = particles[i].position - particles[randomIndex].position + Vector3D(0.0, 0.0, 3 * latticeConstant); // 两粒子间的位移 （z 正方向）
        distance = delta.magnitude();                                                                              // 两粒子间的距离
        if (distance < 2 * particleRadius)
        {
            // 重叠：拒绝位移，回退到原位置
            selectedParticle.position = selectedParticle.position - displacement;
            return; // 结束当前的 metropolisStep
        }

        delta = particles[i].position - particles[randomIndex].position + Vector3D(0.0, 0.0, -3 * latticeConstant); // 两粒子间的位移 （z 负方向）
        distance = delta.magnitude();                                                                               // 两粒子间的距离
        if (distance < 2 * particleRadius)
        {
            // 重叠：拒绝位移，回退到原位置
            selectedParticle.position = selectedParticle.position - displacement;
            return; // 结束当前的 metropolisStep
        }
    }

    // (optional) 为了固定质心，全体粒子反向移动 displacement / numParticles
    for (int i = 0; i < numParticles; ++i)
    {
        particles[i].position = particles[i].position - displacement * (1 / numParticles);
    }

    double energyAfter = calculateTotalEnergy(); // 变化后的势能

    if (randomDouble() < exp((energyBefore - energyAfter) / temperature))
    {
        // 接受位移
    }
    else
    {
        // 拒绝位移，回退到原位置
        selectedParticle.position = selectedParticle.position - displacement;
        for (int i = 0; i < numParticles; ++i)
        {
            particles[i].position = particles[i].position + displacement * (1 / numParticles);
        }
    }
}
