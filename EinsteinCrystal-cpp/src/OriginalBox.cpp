#include <getopt.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>

int main(int argc, char *argv[])
{
    srand((unsigned int)time(0)); // 设置随机数种子

    /* 箱子类型 */
    std::string boxtype;
    // fcc xyz
    int nx = 0;
    int ny = 0;
    int nz = 0;
    // fcc cubic
    int replicate = 0;
    /* 箱子类型 */

    // 粒子直径
    const double sigma = 1.0;

    // 相对密度 rho / rho0
    // const double relativeDensity = 1.0;    // 最密堆积
    const double relativeDensity = 0.7360; // 相对密度 1
    // const double relativeDensity = 0.7778; // 相对密度 2

    // 晶格常数
    const double latticeConstant = sqrt(2.0) * sigma / pow(relativeDensity, 1.0 / 3.0);

    // 温度：k_B T = 1.0
    const double temperature = 1.0;

    int numSteps = 100000;  // 模拟步数
    int sample_period = 50; // 采样周期

    double lambda; // 弹簧系数

    // 命令行参数解析
    int opt;
    while ((opt = getopt(argc, argv, "l:b:x:y:z:r:")) != -1)
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

    // 粒子数
    const int numParticles = (boxtype == "cubic") ? 4 * replicate * replicate * replicate : nx * ny * nz;

    Eigen::Matrix<double, Eigen::Dynamic, 3> particles = Eigen::MatrixXd::Zero(numParticles, 3); // 初始化粒子位置
    Eigen::Matrix<double, Eigen::Dynamic, 3> lattice = Eigen::MatrixXd::Zero(numParticles, 3);   // 初始化晶格位置

    Eigen::RowVector3d axis_1 = Eigen::RowVector3d(0, 0, 0); // 初始化晶格轴
    Eigen::RowVector3d axis_2 = Eigen::RowVector3d(0, 0, 0); // 初始化晶格轴
    Eigen::RowVector3d axis_3 = Eigen::RowVector3d(0, 0, 0); // 初始化晶格轴

    Eigen::RowVector3d v1 = Eigen::RowVector3d(0, 0, 0); // 初始化晶格向量
    Eigen::RowVector3d v2 = Eigen::RowVector3d(0, 0, 0); // 初始化晶格向量
    Eigen::RowVector3d v3 = Eigen::RowVector3d(0, 0, 0); // 初始化晶格向量

    // 初始化盒子长度
    double Lx = 0;
    double Ly = 0;
    double Lz = 0;
    double xy = 0;
    double xz = 0;
    double yz = 0;
    double a2x = 0;
    double a3x = 0;

    // 初始化模拟
    if (boxtype == "ABC")
    {
        if (nx * ny * nz == 0)
        {
            std::cerr << "Usage: " << argv[0] << " [-l lambda] [-b boxtype] [-nx nx] [-ny ny] [-nz nz]" << std::endl;
            std::exit(1);
        }

        const double a0 = latticeConstant / sqrt(2.0);

        axis_1 = Eigen::RowVector3d(1, 0, 0) * a0;
        axis_2 = Eigen::RowVector3d(1.0 / 2.0, sqrt(3.0) / 2.0, 0) * a0;
        axis_3 = Eigen::RowVector3d(1.0 / 2.0, sqrt(3.0) / 6.0, sqrt(2.0 / 3.0)) * a0;

        for (int iz = 0; iz < nz; ++iz)
        {
            for (int iy = 0; iy < ny; ++iy)
            {
                for (int ix = 0; ix < nx; ++ix)
                {
                    lattice.row(iz * nx * ny + iy * nx + ix) << ix * axis_1 + iy * axis_2 + iz * axis_3;
                    particles.row(iz * nx * ny + iy * nx + ix) << ix * axis_1 + iy * axis_2 + iz * axis_3;
                }
            }
        }

        // 计算盒子的形状
        v1 = nx * axis_1;
        v2 = ny * axis_2;
        v3 = nz * axis_3;

        Lx = v1.norm();
        a2x = v1.dot(v2) / Lx;
        Ly = sqrt(v2.squaredNorm() - a2x * a2x);
        xy = a2x / Ly;
        Lz = v3.dot(v1.cross(v2)) / v1.cross(v2).norm();
        a3x = v1.dot(v3) / v1.norm();
        xz = a3x / Lz;
        yz = (v2.dot(v3) - a2x * a3x) / (Ly * Lz);
    }

    // else if (boxtype == "AB")
    // {
    //     if (nx * ny * nz == 0)
    //     {
    //         std::cerr << "Usage: " << argv[0] << " [-l lambda] [-b boxtype] [-nx nx] [-ny ny] [-nz nz]" << std::endl;
    //         std::exit(1);
    //     }

    //     const double a0 = latticeConstant / sqrt(2.0);
    //     const double a1 = latticeConstant / 2.0;

    //     Lx = nx * a0;
    //     Ly = ny * a0;
    //     Lz = nz * a1;

    //     for (int iz = 0; iz < nz; ++iz)
    //     {
    //         for (int iy = 0; iy < ny; ++iy)
    //         {
    //             for (int ix = 0; ix < nx; ++ix)
    //             {
    //                 lattice.row(iz * nx * ny + iy * nx + ix) << a0 * ix + (a0 / 2.0) * fmod(iz, 2),
    //                     a0 * iy + (a0 / 2.0) * fmod(iz, 2),
    //                     a1 * iz;

    //                 particles.row(iz * nx * ny + iy * nx + ix) << a0 * ix + (a0 / 2.0) * fmod(iz, 2),
    //                     a0 * iy + (a0 / 2.0) * fmod(iz, 2),
    //                     a1 * iz;
    //             }
    //         }
    //     }
    // }

    else if (boxtype == "cubic")
    {
        if (replicate == 0)
        {
            std::cerr << "Usage: " << argv[0] << " [-l lambda] [-b boxtype] [-r replicate]" << std::endl;
            std::exit(1);
        }

        double a = latticeConstant;

        axis_1 = Eigen::RowVector3d(1, 0, 0) * a;
        axis_2 = Eigen::RowVector3d(0, 1, 0) * a;
        axis_3 = Eigen::RowVector3d(0, 0, 1) * a;

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

        // 计算盒子的形状
        v1 = replicate * axis_1;
        v2 = replicate * axis_2;
        v3 = replicate * axis_3;

        Lx = v1.norm();
        a2x = v1.dot(v2) / Lx;
        Ly = sqrt(v2.squaredNorm() - a2x * a2x);
        xy = a2x / Ly;
        Lz = v3.dot(v1.cross(v2)) / v1.cross(v2).norm();
        a3x = v1.dot(v3) / v1.norm();
        xz = a3x / Lz;
        yz = (v2.dot(v3) - a2x * a3x) / (Ly * Lz);
    }

    else
    {
        std::cerr << "Unknown box type: " << boxtype << std::endl;
        std::exit(1);
    }

    // 整体平移，使质心初始位于 (0, 0, 0)
    Eigen::RowVector3d Rcm = Eigen::RowVector3d::Zero(1, 3);
    for (int i = 0; i < numParticles; ++i)
    {
        Rcm += particles.row(i);
    }
    Rcm /= numParticles;
    for (int i = 0; i < numParticles; ++i)
    {
        lattice.row(i) -= Rcm;
        particles.row(i) -= Rcm;
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
        // 一次性生成 numParticles 个属于 [0,1) 的均匀分布随机数
        Eigen::MatrixXd randomDouble = Eigen::MatrixXd::Random(1, numParticles);
        randomDouble += Eigen::MatrixXd::Ones(1, numParticles);
        randomDouble /= 2.0;

        for (int selected = 0; selected < numParticles; ++selected)
        {
            // 产生 [-del, del] 之间的随机位移
            Eigen::RowVector3d displacement = Eigen::RowVector3d::Random(1, 3);
            double del = 0.05;
            displacement = displacement * del;

            // 先判断 Metroplois 准则，再判断是否接受位移
            // deltaEnergy = lambda (2 delta r_i dot delta_i + (N-1)/N delta_i^2)
            double deltaEnergy = lambda * 2 * (particles.row(selected) - lattice.row(selected)).dot(displacement) + lambda * (1.0 - 1.0 / numParticles) * displacement.squaredNorm();

            if (randomDouble(selected) < exp(-deltaEnergy / temperature)) // Metroplois 准则
            {
                Eigen::RowVector3d virtualDisplacement = particles.row(selected) + displacement; // 被选中粒子虚拟位移后的位置

                /* 判断是否存在重叠 start */
                bool overlap = false; // 初始值为 false
                for (int i = 0; i < numParticles; ++i)
                {
                    if (i == selected) // 跳过选中粒子
                        continue;

                    Eigen::RowVector3d delta = particles.row(i) - virtualDisplacement; // 粒子间的相对位移

                    // // 划分小盒子，需要判断可能与其相邻的 2x2x2=8 个盒子的周期性边界条件
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

                    // 不划分小盒子，需要判断可能与其相邻的 3x3x3=27 个盒子的周期性边界条件
                    // Eigen::MatrixXd displacements = Eigen::MatrixXd::Zero(27, 3);
                    // for (int dx = -1; dx <= 1; dx++)
                    // {
                    //     for (int dy = -1; dy <= 1; dy++)
                    //     {
                    //         for (int dz = -1; dz <= 1; dz++)
                    //         {
                    //             displacements.row((dx + 1) * 9 + (dy + 1) * 3 + dz + 1) = delta + dx * v1 + dy * v2 + dz * v3;
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

        // std::cout << (particles - lattice).squaredNorm() / numParticles << std::endl; // 输出每一步的均方位移

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

            if (step >= 25000)
            {
                MSD += (particles - lattice).squaredNorm() / numParticles;
            }
        }
    }
    MSD /= (numSteps - 25000) / sample_period;

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

    /* 输出 */
    // 输出箱子类型
    if (boxtype == "cubic")
    {
        std::cout << "boxtype: " << boxtype << std::endl;
        std::cout << "replicate: " << replicate << std::endl;
    }
    else
    {
        std::cout << "boxtype: " << boxtype << std::endl;
        std::cout << "nx: " << nx << " ny: " << ny << " nz: " << nz << std::endl;
    }
    std::cout << "numSteps: " << numSteps << std::endl; // 输出模拟步数
    std::cout << "lambda: " << lambda << std::endl;     // 输出弹簧系数
    std::cout << "MSD: " << MSD << std::endl;           // 输出均方位移

    return 0;
}
