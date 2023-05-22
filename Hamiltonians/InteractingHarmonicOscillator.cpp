#include<memory>
#include <cassert>
#include <iostream>

#include "harmonicoscillator.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

InteractingOscillator::InteractingOscillator(double omega)
{
    assert(omega > 0);
    m_omega  = omega;
}

double InteractingOscillator::computeLocalEnergy(
            class WaveFunction& waveFunction,
            std::vector<std::unique_ptr<class Particle>>& particles
            std::vector<std::vector<double>> rij
        )
{
    /* Here, you need to compute the kinetic and potential energies.
     * Access to the wave function methods can be done using the dot notation
     * for references, e.g., wavefunction.computeDoubleDerivative(particles),
     * to get the Laplacian of the wave function.
     * */

    double r2 = 0;
    unsigned int N = particles.size();
    for (unsigned int i = 0; i < N; i++){
        auto position = particles[i]->getPosition();
        for (unsigned int i = 0; i<particles[0]->getNumberOfDimensions(); i++)
            r2 += position[i]*position[i];
    }
    // m = omega = 1
    double potentialEnergy = 0.5 * r2;
    double kineticEnergy   = waveFunction.computeDoubleDerivative(particles)*-0.5;
    double interactTerm=0; 
    for (unsigned int i= 0;i<N;i++){
        for (unsigned int j>i;j<N;j++){
            interactTerm=interactTerm+rij[i][j]
        }
    }


    return (kineticEnergy + potentialEnergy)/N+interactTerm;
}

std::vector<std::vector<double>>::particlelist(
    std::vector<std::unique_ptr<class Particle>>& particles
)
{
    N=2;
    std::vector<std::vector<double>> rd; 
    for (unsigned int i= 0;i<N;i++){

        std::vector<double> particle1=particles[i]->getPostition();
        for (unsigned int j>i;j<N;j++){
            std::vector<double> particle2=particles[j]->getPostition();
            double r_ij;
            for (unsigned int k=0;k<particles[0]->getNumberOfDimensions(); k++){
                double d_ij=particle1[k]-particle2[k];
                r_ij=r_ij+d_ij*d_ij;


            }
            rd[i][j]=1/r_ij;

        }
    }
    return rd;
}