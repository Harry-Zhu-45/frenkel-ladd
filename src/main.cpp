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

    // 把粒子位置储存为 pos 文件
    if (boxtype == "xyz")
    {
        std::string posFile = "fcc-" + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + ".pos";
        FILE *fp = fopen(posFile.c_str(), "w");

        /* 初始粒子位置 */
        fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2);
        fprintf(fp, "def S0 \"sphere %f \"\n", sigma);
        for (int i = 0; i < numParticles; ++i)
        {
            fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
        }
        fprintf(fp, "eof\n");
        /* 初始粒子位置 */

        double MSD = 0;
        int sample_period = 1000;
        for (int step = 0; step < numSteps; ++step)
        {
            simulation.metropolisStep(); // 模拟

            if ((step + 1) % sample_period == 0)
            {
                /* sample MC 模拟过程中的粒子位置 */
                fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2);
                fprintf(fp, "def S0 \"sphere %f \"\n", sigma); // 粒子形状
                for (int i = 0; i < numParticles; ++i)         // 粒子位置
                {
                    fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
                }
                fprintf(fp, "eof\n");
                /* sample MC 模拟过程中的粒子位置 */

                if (step >= numSteps / 2)
                {
                    MSD += simulation.calculateSD();
                }
            }
        }
        MSD = MSD / (numSteps / sample_period / 2) / numParticles;

        /* 末尾粒子位置 */
        fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2);
        fprintf(fp, "def S0 \"sphere %f \"\n", sigma);
        for (int i = 0; i < numParticles; ++i)
        {
            fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
        }
        fprintf(fp, "eof\n");
        /* 末尾粒子位置 */

        fclose(fp); // 关闭文件

        std::cout << "boxtype: " << boxtype << std::endl;                                  // 箱子类型
        std::cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;          // 粒子数
        std::cout << "numSteps: " << numSteps << std::endl;                                // 模拟步数
        std::cout << "lambda: " << lambda << std::endl;                                    // 弹簧系数
        std::cout << "MSD: " << MSD << std::endl;                                          // 均方位移
        std::cout << "Acceptance ratio: " << simulation.getAcceptanceRatio() << std::endl; // 接受率
    }
    else
    {
        std::string posFile = "fcc-cubic-" + std::to_string(numParticles) + ".pos";
        FILE *fp = fopen(posFile.c_str(), "w");

        /* 初始粒子位置 */
        fprintf(fp, "box %f %f %f\n", replicate * latticeConstant, replicate * latticeConstant, replicate * latticeConstant);
        fprintf(fp, "def S0 \"sphere %f \"\n", sigma);
        for (int i = 0; i < numParticles; ++i)
        {
            fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
        }
        fprintf(fp, "eof\n");
        /* 初始粒子位置 */

        double MSD = 0;
        int sample_period = 1000;
        for (int step = 0; step < numSteps; ++step)
        {
            simulation.metropolisStep(); // 模拟

            if ((step + 1) % sample_period == 0)
            {
                /* sample MC 模拟过程中的粒子位置 */
                fprintf(fp, "box %f %f %f\n", replicate * latticeConstant, replicate * latticeConstant, replicate * latticeConstant);
                fprintf(fp, "def S0 \"sphere %f \"\n", sigma); // 粒子形状
                for (int i = 0; i < numParticles; ++i)         // 粒子位置
                {
                    fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
                }
                fprintf(fp, "eof\n");
                /* sample MC 模拟过程中的粒子位置 */

                if (step >= numSteps / 2)
                {
                    MSD += simulation.calculateSD();
                }
            }
        }
        MSD = MSD / (numSteps / sample_period / 2) / numParticles;

        /* 末尾粒子位置 */
        fprintf(fp, "box %f %f %f\n", replicate * latticeConstant, replicate * latticeConstant, replicate * latticeConstant);
        fprintf(fp, "def S0 \"sphere %f \"\n", sigma);
        for (int i = 0; i < numParticles; ++i)
        {
            fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
        }
        fprintf(fp, "eof\n");
        /* 末尾粒子位置 */

        fclose(fp); // 关闭文件

        std::cout << "boxtype: " << boxtype << std::endl;                                  // 箱子类型
        std::cout << "replicate: " << replicate << std::endl;                              // 复制数
        std::cout << "numSteps: " << numSteps << std::endl;                                // 模拟步数
        std::cout << "lambda: " << lambda << std::endl;                                    // 弹簧系数
        std::cout << "MSD: " << MSD << std::endl;                                          // 均方位移
        std::cout << "Acceptance ratio: " << simulation.getAcceptanceRatio() << std::endl; // 接受率
    }

    return 0;
}
