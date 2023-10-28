#ifndef MONTECARLO_SIMULATION_H
#define MONTECARLO_SIMULATION_H

#include <vector>
#include "Particle.h"

class MonteCarloSimulation
{
public:
    std::vector<Particle> particles;
    std::vector<Particle> lattice;
    MonteCarloSimulation(int numParticles);
    double calculateTotalEnergy();
    void metropolisStep();
};

#endif