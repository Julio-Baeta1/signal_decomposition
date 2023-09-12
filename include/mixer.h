#ifndef MIXER_H
#define MIXER_H

#include <random>
#include <memory>
#include"wave_gen.h"

class Mixer{
    private:
        size_t num_sigs, num_samples;
        std::unique_ptr<Eigen::MatrixXd> mixing_mat;
        std::shared_ptr<Eigen::MatrixXd> raw_sigs,mixed_sigs; 
        

    public:
        Mixer(int num_signals=2, int signal_duration=8);
        Mixer(int num_signals, int signal_duration, Eigen::MatrixXd A);

        Eigen::MatrixXd getRawSignalsValues() const {return *raw_sigs;}
        Eigen::MatrixXd getMixedSignalsValues() const {return *mixed_sigs;}
        Eigen::MatrixXd getMixingMatrixValues() const {return *mixing_mat;}
        std::shared_ptr<Eigen::MatrixXd> getMixedSignalsSharedPtr() {return mixed_sigs;} 
        std::shared_ptr<Eigen::MatrixXd> getRawSignalsSharedPtr() {return raw_sigs;} 

        size_t getNumSignals() const {return num_sigs;}
        size_t getNumSamples() const {return num_samples;}
        void setNumSignals(int new_num_sigs);
        void setNumSamples(int new_num_samps);

        void genSignals();
        void setMixingMatrix(Eigen::MatrixXd A);
        void mixSignals(bool noisy);
};

#endif