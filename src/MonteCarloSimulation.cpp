#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <random>

#include "MonteCarloSimulation.h"
#include "myvariables.h"

int numSteps = 0;    // 模拟步数
double lambda = 0.0; // 弹簧系数

MonteCarloSimulation::MonteCarloSimulation(int numParticles)
{
    lattice.reserve(numParticles);   // 预分配内存
    particles.reserve(numParticles); // 预分配内存

    if (boxtype == "xyz")
    {
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

                    centerOfMass = centerOfMass + lattice[i - 1].position;
                }
            }
        }

        centerOfMass = centerOfMass * (1.0 / numParticles); // 计算质心

        // 固定质心于原点
        for (int i = 0; i < numParticles; ++i)
        {
            lattice[i].position = lattice[i].position - centerOfMass;
            particles[i].position = particles[i].position - centerOfMass;
        }
    }
    else if (boxtype == "cubic")
    {
        double a = latticeConstant;

        for (int iz = 0; iz < replicate; ++iz)
        {
            for (int iy = 0; iy < replicate; ++iy)
            {
                for (int ix = 0; ix < replicate; ++ix)
                {
                    lattice.emplace_back(Vector3D(a * ix, a * iy, a * iz));
                    lattice.emplace_back(Vector3D(a * ix, a * iy + a / 2.0, a * iz + a / 2.0));
                    lattice.emplace_back(Vector3D(a * ix + a / 2.0, a * iy, a * iz + a / 2.0));
                    lattice.emplace_back(Vector3D(a * ix + a / 2.0, a * iy + a / 2.0, a * iz));

                    particles.emplace_back(Vector3D(a * ix, a * iy, a * iz));
                    particles.emplace_back(Vector3D(a * ix, a * iy + a / 2.0, a * iz + a / 2.0));
                    particles.emplace_back(Vector3D(a * ix + a / 2.0, a * iy, a * iz + a / 2.0));
                    particles.emplace_back(Vector3D(a * ix + a / 2.0, a * iy + a / 2.0, a * iz));
                }
            }
        }

        // 计算质心
        Vector3D centerOfMass = Vector3D(0, 0, 0);
        for (int i = 0; i < numParticles; ++i)
        {
            centerOfMass = centerOfMass + lattice[i].position;
        }
        centerOfMass = centerOfMass * (1.0 / numParticles);

        // 固定质心于原点
        for (int i = 0; i < numParticles; ++i)
        {
            lattice[i].position = lattice[i].position - centerOfMass;
            particles[i].position = particles[i].position - centerOfMass;
        }
    }
    else
    {
        std::cerr << "Unknown box type: " << boxtype << std::endl;
        std::exit(1);
    }

    // Initialize counters
    totalMoves = numParticles * numSteps;
    acceptedMoves = 0;
}

double MonteCarloSimulation::getAcceptanceRatio() const
{
    return static_cast<double>(acceptedMoves) / totalMoves;
}

double MonteCarloSimulation::calculateSD()
{
    double totalSD = 0.0;
    for (int i = 0; i < numParticles; ++i)
    {
        totalSD += (particles[i].position - lattice[i].position).dot(particles[i].position - lattice[i].position);
    }
    return totalSD;
}

void MonteCarloSimulation::metropolisStep()
{
    // 创建随机数生成器对象
    std::random_device rd;                         // 用于种子生成
    std::mt19937 gen(rd());                        // Mersenne Twister 伪随机数生成器
    std::uniform_real_distribution<> randomDouble; // 生成属于[0,1)的均匀分布随机数

    // 随机选择粒子，直到每个粒子都被选择一次
    std::vector<bool> particleSelected(numParticles, false); // 初始化为 false，表示所有粒子未被选择

    for (int iteration = 0; iteration < numParticles; ++iteration)
    {
        int randomIndex;
        do
        {
            randomIndex = std::uniform_int_distribution<int>(0, numParticles - 1)(gen);
        } while (particleSelected[randomIndex]);
        particleSelected[randomIndex] = true; // 将该粒子标记为已选择

        // 产生随机位移
        Vector3D displacement(randomDouble(gen) - 0.5, randomDouble(gen) - 0.5, randomDouble(gen) - 0.5);
        double del = 0.12;                 // 随机位移的幅度
        displacement = displacement * del; // 调整随机位移的幅度

        // 先判断 Metroplois 准则，再判断是否接受位移
        // deltaEnergy = lambda (2 delta r_i dot delta_i + (N-1)/N delta_i^2)
        double deltaEnergy = lambda * 2 * (particles[randomIndex].position - lattice[randomIndex].position).dot(displacement) + lambda * (1 - 1.0 / numParticles) * displacement.dot(displacement);

        if (randomDouble(gen) < exp(-deltaEnergy / temperature)) // Metroplois 准则
        {
            Vector3D virtualDisplacement = particles[randomIndex].position + displacement; // 被选中粒子虚拟位移后的位置

            // 判断是否存在重叠
            bool overlap = false; // 初始值为 false

            // 根据粒子位置，将完整盒子分为 2x2x2=8 个小盒子
            int ix = particles[randomIndex].position.x > 0 ? 1 : -1;
            int iy = particles[randomIndex].position.y > 0 ? 1 : -1;
            int iz = particles[randomIndex].position.z > 0 ? 1 : -1;

            for (int i = 0; i < numParticles; ++i)
            {
                if (i != randomIndex) // 遍历除选中粒子外的所有粒子
                {
                    Vector3D delta = particles[i].position - virtualDisplacement; // 粒子间的相对位移

                    // 根据小盒子的划分，只需要判断可能与其相邻的 2x2x2=8 个完整盒子的周期性边界条件
                    Vector3D displacements[8];
                    if (boxtype == "xyz")
                    {
                        for (int dx = 0; dx <= 1; dx++)
                        {
                            for (int dy = 0; dy <= 1; dy++)
                            {
                                for (int dz = 0; dz <= 1; dz++)
                                {
                                    displacements[dx * 4 + dy * 2 + dz] = delta + Vector3D((dx * ix) * nx * latticeConstant / sqrt(2.0), (dy * iy) * ny * latticeConstant / sqrt(2.0), (dz * iz) * nz * latticeConstant / 2.0);
                                }
                            }
                        }
                    }
                    else if (boxtype == "cubic")
                    {
                        for (int dx = 0; dx <= 1; dx++)
                        {
                            for (int dy = 0; dy <= 1; dy++)
                            {
                                for (int dz = 0; dz <= 1; dz++)
                                {
                                    displacements[dx * 4 + dy * 2 + dz] = delta + Vector3D((dx * ix) * replicate * latticeConstant, (dy * iy) * replicate * latticeConstant, (dz * iz) * replicate * latticeConstant);
                                }
                            }
                        }
                    }
                    else
                    {
                        std::cerr << "Unknown box type: " << boxtype << std::endl;
                        std::exit(1);
                    }

                    for (int j = 0; j < 8; j++)
                    {
                        if (displacements[j].dot(displacements[j]) < sigma * sigma)
                        {
                            overlap = true;
                            break;
                        }
                    }
                }
            }

            if (!overlap) // 如果不存在重叠，接受位移
            {
                particles[randomIndex].position = virtualDisplacement;
                for (int i = 0; i < numParticles; ++i)
                {
                    particles[i].position = particles[i].position - displacement * (1.0 / numParticles);
                }
                acceptedMoves++;
            }
        }
    }

    // 重置选择状态
    std::fill(particleSelected.begin(), particleSelected.end(), false);
}
