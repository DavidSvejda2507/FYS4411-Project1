#pragma once
#include<iostream>
#include"wavefunction.h"
#include<memory>
struct rbmParams{
    //a struct to contain the parameters of the boltzmann machine
    std::vector<double> a;
    std::vector<double> b;
    std::vector<std::vector<double>> W;
};
typedef struct rbmParams RBMParams; 
RBMParams initWeights(int nVisible, int Nhidden,  Random * rng);


class RestrictedBoltzmannMachine : public WaveFunction {
    public:
        RestrictedBoltzmannMachine(double sigma, RBMParams * trainableParameters);
        double evaluate(std::vector<std::unique_ptr<class Particle>>& particles);
        void InitialisePositions(std::vector<std::unique_ptr<class Particle>>& particles);
        void adjustPosition(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double> step) ;
        double computeDoubleDerivative(std::vector<std::unique_ptr<class Particle>>& particles);
        std::vector<double> quantumForce(std::vector<std::unique_ptr<class Particle>>& particles, int index);
        std::vector<double> quantumForceMoved(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
        double phiRatio(std::vector<std::unique_ptr<class Particle>>& particles, int index, std::vector<double>& step);
        std::vector<double> getdPhi_dParams(std::vector<std::unique_ptr<class Particle>>& particles);
    private:
    //parameters of restricted boltzmann machine
        double m_sigma ;
        double m_sigma2;
        RBMParams * m_trainableParameters;
        unsigned int m_Nhidden, m_Nvisible;
        double m_gaussianTerm;
        double m_productTerm;
        std::vector<double> m_expBPlusSumXw;//store exp(b_k+m_SumXw)
};