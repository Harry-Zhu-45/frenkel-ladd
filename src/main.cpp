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
            numSteps = std::stoi(optarg); // Parse the number of steps
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

    // // 把粒子位置储存为 pos 文件
    // std::string posFile = "fcc-" + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + ".pos";
    // FILE *fp = fopen(posFile.c_str(), "w");

    // // 盒子大小
    // fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2);
    // // 粒子形状
    // fprintf(fp, "def S0 \"sphere %f \"\n", sigma);
    // for (int i = 0; i < numParticles; ++i) // 粒子位置
    // {
    //     fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
    // }
    // fprintf(fp, "eof\n");

    double MSD = 0;
    for (int step = 0; step < numSteps; ++step)
    {
        simulation.metropolisStep(); // 模拟

        // fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2);
        // fprintf(fp, "def S0 \"sphere %f \"\n", sigma); // 粒子形状
        // for (int i = 0; i < numParticles; ++i)         // 粒子位置
        // {
        //     fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
        // }
        // fprintf(fp, "eof\n");

        // 平均后半段的均方位移
        if (step >= numSteps / 2)
        {
            MSD += simulation.calculateSD();
        }
    }
    MSD = MSD / (numSteps / 2) / numParticles;

    // fclose(fp); // 关闭文件

    std::cout << "numParticles: " << numParticles << std::endl;                        // 粒子数
    std::cout << "numSteps: " << numSteps << std::endl;                                // 模拟步数
    std::cout << "lambda: " << lambda << std::endl;                                    // 弹簧系数
    std::cout << "MSD: " << MSD << std::endl;                                          // 均方位移
    std::cout << "Acceptance ratio: " << simulation.getAcceptanceRatio() << std::endl; // 接受率

    return 0;
}
