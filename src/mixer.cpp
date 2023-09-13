#include"mixer.h"

Mixer::Mixer(int num_signals, int signal_duration)
{
    num_sigs = num_signals;
    num_samples = signal_duration;
    mixing_mat = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_sigs));
    raw_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_samples));
    mixed_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_samples));
}

Mixer::Mixer(int num_signals, int signal_duration, Eigen::MatrixXd A)
{
    num_sigs = num_signals;
    num_samples = signal_duration;
    
    raw_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_samples));
    mixed_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(A.rows(), num_samples));

    setMixingMatrix(A);
}

void Mixer::setMixingMatrix(Eigen::MatrixXd A)
{
    if(A.cols() != raw_sigs->rows())
        throw std::invalid_argument("Mixing matrix A's number of columns must match the number of signals to be mixed");

    mixing_mat = std::make_unique<Eigen::MatrixXd>(A);

    if(mixed_sigs->rows() != A.rows())
        mixed_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(A.rows(),num_samples));
}

void Mixer::setNumSignals(int new_num_sigs)
{

    if(new_num_sigs < 1)
        throw std::invalid_argument("Invalid number of signals, must be greater than 0");

    num_sigs = new_num_sigs;
    raw_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_samples));
    mixing_mat = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(mixing_mat->rows(), num_sigs));
}

void Mixer::setNumSamples(int new_num_samps)
{
    if(new_num_samps < 1)
        throw std::invalid_argument("Invalid sample size, must be greater than 0");

    num_samples = new_num_samps;
    raw_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_samples));
    mixed_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(mixing_mat->rows(), num_samples));
}

void Mixer::genSignals(int seed, int max_amp, int max_period)
{
    //srand( (int)(time(0)) );
    srand(seed);
    SigGen::WaveGen waves((int)num_samples);

    for(int i=0; i<num_sigs; i++)
    {
        raw_sigs->row(i) = waves.genRandomWave(max_amp,max_period); //default max_amp=2 max_period=100
    }

}


void Mixer::mixSignals(bool noisy,int seed)
{
    if(noisy)
    {
        std::mt19937 r_gen (seed);
        std::normal_distribution<double> r_dis(0.0, 1.0);
        auto norm = [&] () {return r_dis(r_gen);};

        Eigen::MatrixXd v = Eigen::MatrixXd::NullaryExpr(num_sigs,num_samples, norm );
        v += *raw_sigs;

        *mixed_sigs = *mixing_mat *v;
    }
    else
        *mixed_sigs = *mixing_mat * *raw_sigs;
}
