#ifndef MIXER_H
#define MIXER_H

#include <random>
#include <memory>
#include"wave_gen.h"
#include <fstream>
#include <sstream>
#include <string>

using Mat = Eigen::MatrixXd;

class Mixer{
    private:
        size_t n_sig, n_samp; //Number of signals and samples respectively
        std::unique_ptr<Mat> mixing_mat;
        std::shared_ptr<Mat> raw_sigs,mixed_sigs;

    public:
        Mixer(int num_signals=2, int signal_duration=8);
        Mixer(int num_signals, int signal_duration, Mat A);

        size_t getNumSignals() const {return n_sig;}
        size_t getNumSamples() const {return n_samp;}
        Mat getRawSignalsValues() const {return *raw_sigs;}
        Mat getMixedSignalsValues() const {return *mixed_sigs;}
        Mat getMixingMatrixValues() const {return *mixing_mat;}

        std::shared_ptr<Mat> getMixedSignalsSharedPtr() {return mixed_sigs;} 
        std::shared_ptr<Mat> getRawSignalsSharedPtr() {return raw_sigs;} 

        void setNumSignalsAndSamples(int new_sigs, int new_samps);
        void setNumSignals(int new_sigs);
        void setNumSamples(int new_samps);

        void setMixingMatrix(Mat A);
        void setMixingFromFile(std::string& file_name);
        void setSourceFromFile(std::string& file_name);

        void genSignals(std::string& gen_file);
        void genSignals(int seed=6, int max_amp=2, int max_period=100);

        void mixSignals(bool noisy=false,int seed=123);
};

#endif