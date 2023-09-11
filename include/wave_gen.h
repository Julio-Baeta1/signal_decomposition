#ifndef WAVE_GEN_H
#define WAVE_GEN_H

/*A very basic wave generator for testing purposes that will generate waves over a set domain based on given period and option.
Will add more functionality in future
options = [sin, cos, square, triangular, sawtooth]*/

#include<Eigen/Dense>
#include<string>
#include<iostream>

namespace SigGen{

enum class SignalType{
    dc=0,sin,cos,square,sawtooth,triangle, size //size is used to get number of elements in SignalType enum class
};

class WaveGen{
    private:
        Eigen::VectorXd domain_sig;
        size_t size;

        Eigen::VectorXd dcWave(double amp);
        Eigen::VectorXd sinWave(double amp, double period);
        Eigen::VectorXd cosWave(double amp, double period);
        Eigen::VectorXd squareWave(double amp, double period);
        Eigen::VectorXd sawtoothWave(double amp, double period);
        Eigen::VectorXd triangleWave(double amp, double period);
        
    public:
        WaveGen(int length=0, double start=0);

        Eigen::VectorXd genWave(double amplitude, double period, std::string opt);
        Eigen::VectorXd genRandomWave(int max_amp, int max_period);

        void setDomain(int new_size, double new_start);
        size_t len() const {return size;}
        
};

}//End SigGen

#endif