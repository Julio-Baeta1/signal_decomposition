#include"mixer.h"

namespace Mixing{

    //srand( (int)(time(0)) );
    //freq = freqs[i];//rand()% max;
    //mat.row(0) = wave_gen.genWave(1,16,"sin");

Mixer::Mixer(int num_signals, int signal_duration, Eigen::MatrixXd mixing_matrix, bool is_noisy)
{
    num_sigs = num_signals;
    num_samples = signal_duration;
    mixing_mat = mixing_matrix;
    noisy = is_noisy;

    raw_sigs = Eigen::MatrixXd::Zero(num_sigs, num_samples);
    mixed_sigs = Eigen::MatrixXd::Zero(num_sigs, num_samples);
}

void Mixer::genSignals()
{
    int x=0;
}

void Mixer::setMixingMatrix()
{
    int x=0;
}

void mixSignals()
{
    int x=0;
}

}