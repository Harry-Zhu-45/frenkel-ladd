#ifndef MONTECARLO_SIMULATION_H
#define MONTECARLO_SIMULATION_H

#include <vector>
#include <random>
#include "Particle.h"

class MonteCarloSimulation
{
private:
    std::string boxtype; // Type of box
    int nx, ny, nz;      // Number of unit cells in each direction
    int replicate;       // Number of replications in each direction
    double Lx, Ly, Lz;   // Box dimensions

    int numParticles;       // Number of particles
    double latticeConstant; // Lattice constant
    double temperature;     // Temperature
    int numSteps;           // Number of Monte Carlo steps
    double lambda;          // Spring constant
    double sigma;           // diameter of the particles

    int totalMoves;         // Total number of move attempts
    int acceptedMoves;      // Number of accepted moves
    double acceptanceRatio; // Acceptance ratio
public:
    std::vector<Particle> particles;
    std::vector<Particle> lattice;

    // Constructor
    MonteCarloSimulation(std::string boxtype,
                         int nx, int ny, int nz,
                         int replicate,
                         double latticeConstant,
                         double temperature,
                         int numSteps,
                         double lambda,
                         double sigma)
    {
        this->boxtype = boxtype;
        this->nx = nx;
        this->ny = ny;
        this->nz = nz;
        this->replicate = replicate;

        if (this->boxtype == "xyz")
        {
            this->numParticles = nx * ny * nz;
        }
        else if (this->boxtype == "cubic")
        {
            this->numParticles = 4 * replicate * replicate * replicate;
        }
        else
        {
            std::cerr << "Unknown box type: " << boxtype << std::endl;
            std::exit(1);
        }

        this->latticeConstant = latticeConstant;
        this->temperature = temperature;
        this->numSteps = numSteps;
        this->lambda = lambda;
        this->sigma = sigma;

        this->acceptedMoves = 0;
        this->totalMoves = numParticles * numSteps;

        this->lattice.reserve(numParticles);   // 预分配内存
        this->particles.reserve(numParticles); // 预分配内存

        if (boxtype == "xyz")
        {
            const double a1 = this->latticeConstant;
            const double a0 = this->latticeConstant / sqrt(2.0);

            this->Lx = nx * a0;
            this->Ly = ny * a0;
            this->Lz = nz * a1 / 2.0;

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

            this->Lx = replicate * a;
            this->Ly = replicate * a;
            this->Lz = replicate * a;

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
    }

    double getAcceptanceRatio() const
    {
        return static_cast<double>(acceptedMoves) / totalMoves;
    }

    double calculateMSD()
    {
        double totalSD = 0.0;
        for (int i = 0; i < this->numParticles; ++i)
        {
            totalSD += (particles[i].position - lattice[i].position).dot(particles[i].position - lattice[i].position);
        }
        return totalSD / this->numParticles;
    }

    void metropolisStep()
    {
        // 创建随机数生成器对象
        std::random_device rd;                         // 用于种子生成
        std::mt19937 gen(rd());                        // Mersenne Twister 伪随机数生成器
        std::uniform_real_distribution<> randomDouble; // 生成属于 [0,1) 的均匀分布随机数

        for (int selected = 0; selected < numParticles; ++selected)
        {
            // 产生随机位移
            Vector3D displacement(randomDouble(gen) - 0.5, randomDouble(gen) - 0.5, randomDouble(gen) - 0.5);
            double del = 0.12;                 // 随机位移的幅度
            displacement = displacement * del; // 调整随机位移的幅度

            // 先判断 Metroplois 准则，再判断是否接受位移
            // deltaEnergy = lambda (2 delta r_i dot delta_i + (N-1)/N delta_i^2)
            double deltaEnergy = lambda * 2 * (particles[selected].position - lattice[selected].position).dot(displacement) + lambda * (1.0 - 1.0 / numParticles) * displacement.dot(displacement);

            if (randomDouble(gen) < exp(-deltaEnergy / temperature)) // Metroplois 准则
            {
                Vector3D virtualDisplacement = particles[selected].position + displacement; // 被选中粒子虚拟位移后的位置

                // 判断是否存在重叠
                bool overlap = false; // 初始值为 false

                for (int i = 0; i < numParticles; ++i)
                {
                    if (i != selected) // 遍历除选中粒子外的所有粒子
                    {
                        Vector3D delta = particles[i].position - virtualDisplacement; // 粒子间的相对位移

                        // 根据小盒子的划分，只需要判断可能与其相邻的 2x2x2=8 个完整盒子的周期性边界条件
                        int ix = particles[selected].position.x > 0 ? 1 : -1;
                        int iy = particles[selected].position.y > 0 ? 1 : -1;
                        int iz = particles[selected].position.z > 0 ? 1 : -1;

                        Vector3D displacements[8];
                        for (int dx = 0; dx <= 1; dx++)
                        {
                            for (int dy = 0; dy <= 1; dy++)
                            {
                                for (int dz = 0; dz <= 1; dz++)
                                {
                                    displacements[dx * 4 + dy * 2 + dz] = delta + Vector3D((dx * ix) * this->Lx, (dy * iy) * this->Ly, (dz * iz) * this->Lz);
                                }
                            }
                        }

                        for (int j = 0; j < 8; j++)
                        {
                            if (displacements[j].dot(displacements[j]) < sigma * sigma)
                            {
                                overlap = true;
                                break;
                            }
                        }

                        // // 不划分小盒子，需要判断可能与其相邻的 3x3x3=27 盒子的周期性边界条件
                        // Vector3D displacements[27];
                        // for (int dx = -1; dx <= 1; dx++)
                        // {
                        //     for (int dy = -1; dy <= 1; dy++)
                        //     {
                        //         for (int dz = -1; dz <= 1; dz++)
                        //         {
                        //             displacements[(dx + 1) * 9 + (dy + 1) * 3 + dz + 1] = delta + Vector3D(dx * this->Lx, dy * this->Ly, dz * this->Lz);
                        //         }
                        //     }
                        // }

                        // for (int j = 0; j < 27; j++)
                        // {
                        //     if (displacements[j].dot(displacements[j]) < sigma * sigma)
                        //     {
                        //         overlap = true;
                        //         break;
                        //     }
                        // }
                    }
                }

                if (!overlap) // 如果不存在重叠，接受位移
                {
                    particles[selected].position = virtualDisplacement;
                    for (int i = 0; i < numParticles; ++i)
                    {
                        particles[i].position = particles[i].position - displacement * (1.0 / numParticles);
                    }
                    acceptedMoves++;
                }
            }
        }
    }
};

#endif
