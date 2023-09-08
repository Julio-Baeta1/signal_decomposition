#include"mixer.h"

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
    //srand( (int)(time(0)) );
    //freq = freqs[i];//rand()% max;
    srand(6);
    SigGen::WaveGen waves(num_samples,0);

    for(int i=0; i<num_sigs; i++)
    {
        raw_sigs.row(i) = waves.genRandomWave(2,100);
    }

}

void Mixer::setMixingMatrix()
{
    int x=0;
}

void mixSignals()
{
    int x=0;
}
