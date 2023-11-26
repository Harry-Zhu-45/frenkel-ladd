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

    double MSD = 0;
    for (int step = 0; step < numSteps; ++step)
    {
        simulation.metropolisStep(); // 模拟

        // 平均最后 10000 步均方位移
        if (step >= numSteps - 10000)
        {
            MSD += simulation.calculateSD();
        }
    }
    MSD = MSD / 10000 / numParticles;

    std::cout << "numParticles: " << numParticles << std::endl; // 粒子数
    std::cout << "numSteps: " << numSteps << std::endl;         // 模拟步数
    std::cout << "lambda: " << lambda << std::endl;             // 弹簧系数
    std::cout << "MSD: " << MSD << std::endl;                   // 均方位移

    // 储存粒子位置为 injavis 可视化软件可以加载的 `.pos` 文件
    // 如果没有就创建一个，如果有就覆盖
    FILE *fp = fopen("fcc-3x4x6.pos", "w");

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
