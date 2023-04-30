#pragma once

#include "particle.hpp"
#include "inlet.hpp"

#include "mpi.h"


/**
 * @brief Межпроцессорный обмен.
 * @details Класс-интерфейс для вызова функций MPI обмена.
 * 
 */
class MPIw{
private:
    int rank;
    int size;
    int count;
    
    MPI_Status status;
    
    std::vector<double> ex;
    std::vector<double>::iterator d_it;
    
    std::vector<Particle> mpi_income;
    std::vector<Particle>::iterator iter_0, iter_1;
    
    std::vector<double> penvec;
	std::vector<double> massvec;
    std::vector<std::pair<double&, double&>> sortvec;
    double all_mass = 0;
    double mass = 0;
    double pen = 0;

    bool here, result;

public:
    MPIw(int* argc, char*** argv);

    ~MPIw();

    int getRank() const;
    
    int getSize() const;

    void send(std::vector<Particle>::iterator begin,std::vector<Particle>::iterator end, int address);
    
    void recv(std::vector<Particle>& particles, int address, ParticleSubstance* substance);

    void send(const double& d, int address);

    void recv(double* d, int address);

    void send_vec(const std::vector<double>& vec, int address, int tag_);

    void recv_vec(std::vector<double>& vec, int address, int tag_);

    void exchangeBlocking(std::vector<Particle>& particles, ParticleSubstance* particle_substance);

    double pencalc(const std::vector<Particle>& particles, const Inlet& inlet);

    bool particlesToSolve(std::vector<Particle>& particles);
};