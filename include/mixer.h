#ifndef MIXER_H
#define MIXER_H

#include <random>
#include"wave_gen.h"

class Mixer{
    private:
        Eigen::MatrixXd raw_sigs,mixed_sigs,mixing_mat; //Will make pointers
        size_t num_sigs, num_samples;
        bool noisy;

    public:
        Mixer(int num_signals=8, int signal_duration=8, Eigen::MatrixXd mixing_matrix=Eigen::MatrixXd::Identity(8,8), bool is_noisy=false);

        Eigen::MatrixXd getRawSignals() const {return raw_sigs;}
        Eigen::MatrixXd getMixedSignals() const {return mixed_sigs;}
        Eigen::MatrixXd getMixingMatrix() const {return mixing_mat;}
        bool isNoisy() const {return noisy;}

        size_t getNumSignals() const {return num_sigs;}
        size_t getNumSamples() const {return num_samples;}

        void genSignals();
        void setMixingMatrix();
        void mixSignals();
};

#endif