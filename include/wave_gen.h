#ifndef WAVE_GEN_H
#define WAVE_GEN_H

/*A very basic wave generator for testing purposes that will generate waves over a set domain based on given period and option.
Will add more functionality in future
options = [sin, cos, square, triangular, sawtooth]*/

#include<Eigen/Dense>
#include<string>
#include<iostream>

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
        WaveGen(int length=0, double start=0, double end=100);
        Eigen::VectorXd genWave(double amplitude, double period, std::string opt);
};


#endif