#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <getopt.h>
#include <fstream>
#include <string>
#include <sstream>
#include <random>

#include "myvariables.h"
#include "Vector3D.h"
#include "Particle.h"
#include "MonteCarloSimulation.h"

int main(int argc, char *argv[])
{
    // 读取 lattice 文件
    std::ifstream inputFile("src/lattice/" + latticeFileName);
    if (!inputFile.is_open())
    {
        // 检查文件是否打开成功，如果没有打开成功，输出错误信息并退出
        std::cerr << "Failed to open lattice.txt" << std::endl;
        return 1;
    }
    inputFile >> numParticles;
    inputFile.close();
    // 关闭 lattice 文件

    // srand(static_cast<unsigned>(time(nullptr))); // 设置随机数种子
    // numSteps = 10000 * numParticles;             // 设置默认模拟步数
    // lambda = 7.0;                                // 设置默认弹簧系数

    // 命令行参数解析
    int opt;
    while ((opt = getopt(argc, argv, "n:l:")) != -1)
    {
        switch (opt)
        {
        case 'n':
            numSteps = std::stoi(optarg) * numParticles; // Parse the number of steps
            break;
        case 'l':
            lambda = std::stod(optarg); // Parse the lambda value
            break;
        default:
            std::cerr << "Usage: " << argv[0] << " [-n numSteps] [-l lambda]" << std::endl;
            return 1;
        }
    }

    MonteCarloSimulation simulation(numParticles); // 初始化模拟

    for (int step = 0; step < numSteps; ++step)
    {
        simulation.metropolisStep(); // 模拟
    }

    std::cout << "numParticles: " << numParticles << std::endl;                      // 输出粒子数
    std::cout << "numSteps: " << numSteps << std::endl;                              // 输出模拟步数
    std::cout << "lambda: " << lambda << std::endl;                                  // 输出弹簧系数
    std::cout << "Total energy: " << simulation.calculateTotalEnergy() << std::endl; // 输出总势能

    // 储存粒子位置为 injavis 可视化软件可以加载的 `.pos` 文件
    // 如果没有就创建一个，如果有就覆盖
    // 文件名与 lattice 文件名相同，后缀名为 .pos
    FILE *fp = fopen((latticeFileName.substr(0, latticeFileName.find_last_of('.')) + ".pos").c_str(), "w");

    // fprintf(fp, "box %f %f %f\n", latticeConstant, latticeConstant, latticeConstant);
    fprintf(fp, "box %f %f %f\n", 2 * latticeConstant, 3 * latticeConstant, 3 * latticeConstant); // box 大小
    fprintf(fp, "def S0 \"sphere %f \"\n", 2 * particleRadius);                                   // 定义粒子形状
    for (int i = 0; i < numParticles; ++i)
    {
        fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z); // 输出粒子位置
    }
    fprintf(fp, "eof\n"); // 结束 .pos 文件
    fclose(fp);

    return 0;
}
