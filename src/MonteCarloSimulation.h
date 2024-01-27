#ifndef MONTECARLO_SIMULATION_H
#define MONTECARLO_SIMULATION_H

#include <vector>
#include "Particle.h"

class MonteCarloSimulation
{
private:
    int totalMoves;    // Total number of move attempts
    int acceptedMoves; // Number of accepted moves
public:
    std::vector<Particle> particles;
    std::vector<Particle> lattice;
    MonteCarloSimulation(int numParticles);
    double getAcceptanceRatio() const;
    double calculateSD();
    void metropolisStep();
};

#endif