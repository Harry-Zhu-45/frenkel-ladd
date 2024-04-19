#include <getopt.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

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

    // 粒子直径
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
        // case 'n':
        //     numSteps = std::stoi(optarg); // Parse the number of steps
        //     break;
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

    // 粒子数
    const int numParticles = (boxtype == "xyz") ? nx * ny * nz : 4 * replicate * replicate * replicate;

    Eigen::Matrix<double, Eigen::Dynamic, 3> particles = Eigen::MatrixXd::Zero(numParticles, 3); // 初始化粒子位置
    Eigen::Matrix<double, Eigen::Dynamic, 3> lattice = Eigen::MatrixXd::Zero(numParticles, 3);   // 初始化晶格位置

    // 初始化模拟
    if (boxtype == "xyz")
    {
        if (nx * ny * nz == 0)
        {
            std::cerr << "Usage: " << argv[0] << " [-l lambda] [-b boxtype] [-nx nx] [-ny ny] [-nz nz]" << std::endl;
            std::exit(1);
        }

        const double a0 = latticeConstant / sqrt(2.0);
        const double a1 = latticeConstant / 2.0;

        Lx = nx * a0;
        Ly = ny * a0;
        Lz = nz * a1;

        for (int iz = 0; iz < nz; ++iz)
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    lattice.row(iz * nx * ny + iy * nx + ix) << a0 * ix + (a0 / 2.0) * fmod(iz, 2),
                        a0 * iy + (a0 / 2.0) * fmod(iz, 2),
                        a1 * iz;

                    particles.row(iz * nx * ny + iy * nx + ix) << a0 * ix + (a0 / 2.0) * fmod(iz, 2),
                        a0 * iy + (a0 / 2.0) * fmod(iz, 2),
                        a1 * iz;
                }
            }
        }
    }

    else if (boxtype == "cubic")
    {
        if (replicate == 0)
        {
            std::cerr << "Usage: " << argv[0] << " [-l lambda] [-b boxtype] [-r replicate]" << std::endl;
            std::exit(1);
        }

        double a = latticeConstant;

        Lx = replicate * a;
        Ly = replicate * a;
        Lz = replicate * a;

        int i = 0; // 一开始没有粒子
        for (int iz = 0; iz < replicate; ++iz)
        {
            for (int iy = 0; iy < replicate; ++iy)
            {
                for (int ix = 0; ix < replicate; ++ix)
                {
                    lattice.row(4 * (iz * replicate * replicate + iy * replicate + ix)) << a * ix, a * iy, a * iz;
                    lattice.row(4 * (iz * replicate * replicate + iy * replicate + ix) + 1) << a * ix, a * iy + a / 2.0, a * iz + a / 2.0;
                    lattice.row(4 * (iz * replicate * replicate + iy * replicate + ix) + 2) << a * ix + a / 2.0, a * iy, a * iz + a / 2.0;
                    lattice.row(4 * (iz * replicate * replicate + iy * replicate + ix) + 3) << a * ix + a / 2.0, a * iy + a / 2.0, a * iz;

                    particles.row(4 * (iz * replicate * replicate + iy * replicate + ix)) << a * ix, a * iy, a * iz;
                    particles.row(4 * (iz * replicate * replicate + iy * replicate + ix) + 1) << a * ix, a * iy + a / 2.0, a * iz + a / 2.0;
                    particles.row(4 * (iz * replicate * replicate + iy * replicate + ix) + 2) << a * ix + a / 2.0, a * iy, a * iz + a / 2.0;
                    particles.row(4 * (iz * replicate * replicate + iy * replicate + ix) + 3) << a * ix + a / 2.0, a * iy + a / 2.0, a * iz;
                }
            }
        }
    }

    else
    {
        std::cerr << "Unknown box type: " << boxtype << std::endl;
        std::exit(1);
    }

    // 质心初始位于原点
    Eigen::Matrix<double, 1, 3> centerOfMass = Eigen::Matrix<double, 1, 3>::Zero(1, 3);
    for (int i = 0; i < numParticles; ++i)
    {
        centerOfMass += particles.row(i);
    }
    centerOfMass /= numParticles;
    for (int i = 0; i < numParticles; ++i)
    {
        lattice.row(i) -= centerOfMass;
        particles.row(i) -= centerOfMass;
    }

    // // 准备要写入的 pos 文件，但是不一定写入
    // std::ofstream fout;
    // fout.open("trajectory.pos"); // 打开文件

    // /* 初始粒子位置 */
    // fout << "box " << Lx << " " << Ly << " " << Lz << std::endl;
    // fout << "def S0 \"sphere " << sigma << " \"" << std::endl;
    // for (int i = 0; i < numParticles; ++i)
    // {
    //     fout << "S0 ffff0000 " << particles(i, 0) << " " << particles(i, 1) << " " << particles(i, 2) << std::endl;
    // }
    // fout << "eof" << std::endl;
    // /* 初始粒子位置 */

    double MSD = 0; // 均方位移
    for (int step = 0; step < numSteps; ++step)
    {
        Eigen::RowVector3d randomDouble = Eigen::RowVector3d::Random(1, numParticles);
        randomDouble = (randomDouble.array() + 1.0) / 2.0; // 一次性生成 numParticles 个属于 [0,1) 的均匀分布随机数

        for (int selected = 0; selected < numParticles; ++selected)
        {
            // 产生 [-del, del] 之间的随机位移
            Eigen::Matrix<double, 1, 3> displacement = Eigen::Matrix<double, 1, 3>::Random(1, 3);
            double del = 0.06;
            displacement = displacement * del;

            // 先判断 Metroplois 准则，再判断是否接受位移
            // deltaEnergy = lambda (2 delta r_i dot delta_i + (N-1)/N delta_i^2)
            double deltaEnergy = lambda * 2 * (particles.row(selected) - lattice.row(selected)).dot(displacement) + lambda * (1.0 - 1.0 / numParticles) * displacement.dot(displacement);

            if (randomDouble(selected) < exp(-deltaEnergy / temperature)) // Metroplois 准则
            {
                Eigen::Matrix<double, 1, 3> virtualDisplacement = particles.row(selected) + displacement; // 被选中粒子虚拟位移后的位置

                /* 判断是否存在重叠 start */
                bool overlap = false; // 初始值为 false
                for (int i = 0; i < numParticles; ++i)
                {
                    if (i == selected) // 跳过选中粒子
                        continue;

                    Eigen::Matrix<double, 1, 3> delta = particles.row(i) - virtualDisplacement; // 粒子间的相对位移

                    // 划分小盒子，需要判断可能与其相邻的 2x2x2=8 个盒子的周期性边界条件
                    Eigen::Matrix<double, 8, 3> displacements = Eigen::MatrixXd::Zero(8, 3);
                    int ix = particles(selected, 0) > 0 ? 1 : -1;
                    int iy = particles(selected, 1) > 0 ? 1 : -1;
                    int iz = particles(selected, 2) > 0 ? 1 : -1;
                    for (int dx = 0; dx <= 1; dx++)
                    {
                        for (int dy = 0; dy <= 1; dy++)
                        {
                            for (int dz = 0; dz <= 1; dz++)
                            {
                                displacements.row(dx * 4 + dy * 2 + dz) = delta + Eigen::RowVector3d((dx * ix) * Lx, (dy * iy) * Ly, (dz * iz) * Lz);
                            }
                        }
                    }

                    // // 不划分小盒子，需要判断可能与其相邻的 3x3x3=27 个盒子的周期性边界条件
                    // Eigen::MatrixXd displacements = Eigen::MatrixXd::Zero(27, 3);
                    // for (int dx = -1; dx <= 1; dx++)
                    // {
                    //     for (int dy = -1; dy <= 1; dy++)
                    //     {
                    //         for (int dz = -1; dz <= 1; dz++)
                    //         {
                    //             displacements.row((dx + 1) * 9 + (dy + 1) * 3 + dz + 1) = delta + Eigen::RowVector3d(dx * Lx, dy * Ly, dz * Lz);
                    //         }
                    //     }
                    // }

                    if (displacements.rowwise().squaredNorm().minCoeff() < sigma * sigma)
                    {
                        overlap = true;
                        break;
                    }
                }
                /* 判断是否存在重叠 end */

                if (!overlap) // 如果不存在重叠，接受位移
                {
                    particles.row(selected) = virtualDisplacement;
                    for (int i = 0; i < numParticles; ++i)
                    {
                        particles.row(i) -= displacement / numParticles;
                    }
                }
            }
        }

        // 输出每一步的均方位移
        double msd = (particles - lattice).squaredNorm() / numParticles; // 当前step的均方位移
        std::cout << msd << std::endl;

        if ((step + 1) % sample_period == 0) // sample
        {
            // /* sample MC 模拟过程中的粒子位置 */
            // fout << "box " << Lx << " " << Ly << " " << Lz << std::endl;
            // fout << "def S0 \"sphere " << sigma << " \"" << std::endl;
            // for (int i = 0; i < numParticles; ++i)
            // {
            //     fout << "S0 ffff0000 " << particles(i, 0) << " " << particles(i, 1) << " " << particles(i, 2) << std::endl;
            // }
            // fout << "eof" << std::endl;
            // /* sample MC 模拟过程中的粒子位置 */

            if (step >= numSteps / 2)
            {
                MSD += msd;
            }
        }
    }

    MSD = MSD / (numSteps / 2 / sample_period);

    // /* 末尾粒子位置 */
    // fout << "box " << Lx << " " << Ly << " " << Lz << std::endl;
    // fout << "def S0 \"sphere " << sigma << " \"" << std::endl;
    // for (int i = 0; i < numParticles; ++i)
    // {
    //     fout << "S0 ffff0000 " << particles(i, 0) << " " << particles(i, 1) << " " << particles(i, 2) << std::endl;
    // }
    // fout << "eof" << std::endl;
    // /* 末尾粒子位置 */

    // fout.close(); // 关闭文件

    // 输出箱子类型
    if (boxtype == "xyz")
    {
        std::cout << "boxtype: " << boxtype << std::endl;
        std::cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;
    }
    else if (boxtype == "cubic")
    {
        std::cout << "boxtype: " << boxtype << std::endl;
        std::cout << "replicate: " << replicate << std::endl;
    }
    std::cout << "numSteps: " << numSteps << std::endl; // 输出模拟步数
    std::cout << "lambda: " << lambda << std::endl;     // 输出弹簧系数
    std::cout << "MSD: " << MSD << std::endl;           // 输出均方位移

    return 0;
}