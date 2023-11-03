#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <random>

#include "MonteCarloSimulation.h"
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
                        displacement[(dx + 1) * 9 + (dy + 1) * 3 + (dz + 1)] = delta + Vector3D(dx * 2 * latticeConstant, dy * 3 * latticeConstant, dz * 3 * latticeConstant);
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

    // 不重叠，计算势能变化
    double deltaEnergy = lambda * (virtualDisplacement - lattice[randomIndex].position).magnitude() * (virtualDisplacement - lattice[randomIndex].position).magnitude() -
                         lambda * (particles[randomIndex].position - lattice[randomIndex].position).magnitude() * (particles[randomIndex].position - lattice[randomIndex].position).magnitude();

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