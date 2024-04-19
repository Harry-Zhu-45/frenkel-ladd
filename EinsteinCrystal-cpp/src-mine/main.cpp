#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <getopt.h>
#include <fstream>
#include <string>
#include <sstream>
#include <random>

#include "Vector3D.h"
#include "Particle.h"
#include "MonteCarloSimulation.h"

int main(int argc, char *argv[])
{
    /* 箱子类型 */
    std::string boxtype;
    // fcc xyz
    int nx = 0;
    int ny = 0;
    int nz = 0;
    // fcc cubic
    int replicate = 0;
    double Lx = 0;
    double Ly = 0;
    double Lz = 0;
    /* 箱子类型 */

    int numParticles = (boxtype == "xyz") ? nx * ny * nz : 4 * replicate * replicate * replicate; // 粒子数

    const double sigma = 1.0;

    // 相对密度 rho / rho0
    const double relativeDensity = 0.7360;
    // 晶格常数
    const double latticeConstant = sqrt(2.0) * sigma / pow(relativeDensity, 1.0 / 3.0);

    // 温度：k_B T = 1.0
    const double temperature = 1.0;

    int numSteps = 1000000;   // 模拟步数
    int sample_period = 1000; // 采样周期

    double lambda; // 弹簧系数

    // 命令行参数解析
    int opt;
    while ((opt = getopt(argc, argv, "n:l:b:x:y:z:r:")) != -1)
    {
        switch (opt)
        {
        case 'l':
            lambda = std::stod(optarg); // Parse the lambda value
            break;
        case 'b':
            boxtype = optarg; // Parse the box type
            break;
        case 'x':
            nx = std::stoi(optarg); // Parse the number of particles in x direction
            break;
        case 'y':
            ny = std::stoi(optarg); // Parse the number of particles in y direction
            break;
        case 'z':
            nz = std::stoi(optarg); // Parse the number of particles in z direction
            break;
        case 'r':
            replicate = std::stoi(optarg); // Parse the replicate number
            break;
        default:
            std::cerr << "Wrong argument!" << std::endl;
            std::exit(1);
        }
    }

    // 初始化模拟
    MonteCarloSimulation simulation(boxtype, nx, ny, nz, replicate, latticeConstant,
                                    temperature, numSteps, lambda, sigma);

    // std::string posFile = "trajectory.pos";
    // FILE *fp = fopen(posFile.c_str(), "w");

    // /* 初始粒子位置 */
    // fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2.0);
    // fprintf(fp, "def S0 \"sphere %f \"\n", sigma);
    // for (int i = 0; i < numParticles; ++i)
    // {
    //     fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
    // }
    // fprintf(fp, "eof\n");
    // /* 初始粒子位置 */

    double MSD = 0;
    for (int step = 0; step < numSteps; ++step)
    {
        simulation.metropolisStep(); // 模拟

        // 输出每一步的均方位移
        double msd = simulation.calculateMSD(); // 当前step的均方位移
        std::cout << msd << std::endl;

        if ((step + 1) % sample_period == 0)
        {
            // /* sample MC 模拟过程中的粒子位置 */
            // fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2.0);
            // fprintf(fp, "def S0 \"sphere %f \"\n", sigma); // 粒子形状
            // for (int i = 0; i < numParticles; ++i)         // 粒子位置
            // {
            //     fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
            // }
            // fprintf(fp, "eof\n");
            // /* sample MC 模拟过程中的粒子位置 */

            if (step >= numSteps / 2)
            {
                MSD += simulation.calculateMSD();
            }
        }
    }
    MSD = MSD / (numSteps / sample_period / 2);

    // /* 末尾粒子位置 */
    // fprintf(fp, "box %f %f %f\n", nx * latticeConstant / sqrt(2.0), ny * latticeConstant / sqrt(2.0), nz * latticeConstant / 2.0);
    // fprintf(fp, "def S0 \"sphere %f \"\n", sigma);
    // for (int i = 0; i < numParticles; ++i)
    // {
    //     fprintf(fp, "S0 ffff0000 %f %f %f\n", simulation.particles[i].position.x, simulation.particles[i].position.y, simulation.particles[i].position.z);
    // }
    // fprintf(fp, "eof\n");
    // /* 末尾粒子位置 */

    // fclose(fp); // 关闭文件

    // std::cout << "boxtype: " << boxtype << std::endl;                                  // 箱子类型
    // std::cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;          // 粒子数

    // std::cout << "boxtype: " << boxtype << std::endl;     // 箱子类型
    // std::cout << "replicate: " << replicate << std::endl; // 复制数

    std::cout << "numSteps: " << numSteps << std::endl;                                // 模拟步数
    std::cout << "lambda: " << lambda << std::endl;                                    // 弹簧系数
    std::cout << "MSD: " << MSD << std::endl;                                          // 均方位移
    std::cout << "Acceptance ratio: " << simulation.getAcceptanceRatio() << std::endl; // 接受率

    return 0;
}
