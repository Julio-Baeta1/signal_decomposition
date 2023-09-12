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
        throw std::invalid_argument("Mixing matrix A's columns must match the number of signals to be mixed");

    mixing_mat = std::make_unique<Eigen::MatrixXd>(A);

    if(mixed_sigs->rows() != A.rows())
        mixed_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(A.rows(),num_samples));
}

void Mixer::setNumSignals(int new_num_sigs)
{
    num_sigs = new_num_sigs;
    raw_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_samples));
    if(mixing_mat->cols() != new_num_sigs)
        mixing_mat = std::make_unique<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_sigs));
}

void Mixer::setNumSamples(int new_num_samps)
{
    num_samples = new_num_samps;
    raw_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(num_sigs, num_samples));
    mixed_sigs = std::make_shared<Eigen::MatrixXd>(Eigen::MatrixXd::Zero(mixing_mat->rows(), num_samples));
}

void Mixer::genSignals()
{
    //srand( (int)(time(0)) );
    srand(6);
    SigGen::WaveGen waves(num_samples);

    for(int i=0; i<num_sigs; i++)
    {
        raw_sigs->row(i) = waves.genRandomWave(2,100); //max_amp=2 max_period=100
    }

}


void Mixer::mixSignals(bool noisy)
{
    if(noisy)
    {
        std::mt19937 r_gen (123);
        std::normal_distribution<double> r_dis(0.0, 1.0);
        auto norm = [&] () {return r_dis(r_gen);};

        Eigen::MatrixXd v = Eigen::MatrixXd::NullaryExpr(num_sigs,num_samples, norm );
        v += *raw_sigs;

        *mixed_sigs = *mixing_mat *v;
    }
    else
        *mixed_sigs = *mixing_mat * *raw_sigs;
}
