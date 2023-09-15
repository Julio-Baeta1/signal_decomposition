#include"mixer.h"

Mixer::Mixer(int num_signals, int signal_duration)
{
    mixing_mat = std::make_unique<Mat>(Mat::Zero(num_signals, signal_duration));
    setNumSignalsAndSamples(num_signals,signal_duration);
}

Mixer::Mixer(int num_signals, int signal_duration, Eigen::MatrixXd A)
{
    mixing_mat = std::make_unique<Mat>(Mat::Zero(num_signals, signal_duration));
    setNumSignalsAndSamples(num_signals,signal_duration);
    setMixingMatrix(A);
}

void Mixer::setNumSignalsAndSamples(int new_num_sigs, int new_num_samps)
{
    if(new_num_sigs < 1)
        throw std::invalid_argument("Invalid number of signals, must be greater than 0");
    if(new_num_samps < 1)
        throw std::invalid_argument("Invalid sample size, must be greater than 0");

    n_sig = new_num_sigs;
    n_samp = new_num_samps;
        
    raw_sigs = std::make_shared<Mat>(Mat::Zero(n_sig, n_samp));
    mixing_mat = std::make_unique<Mat>(Mat::Zero(mixing_mat->rows(), n_sig));
    mixed_sigs = std::make_shared<Mat>(Mat::Zero(mixing_mat->rows(), n_samp));
}

void Mixer::setNumSignals(int new_num_sigs)
{

    if(new_num_sigs < 1)
        throw std::invalid_argument("Invalid number of signals, must be greater than 0");

    n_sig = new_num_sigs;
    raw_sigs = std::make_shared<Mat>(Mat::Zero(n_sig, n_samp));
    mixing_mat = std::make_unique<Mat>(Mat::Zero(mixing_mat->rows(), n_sig));
}

void Mixer::setNumSamples(int new_num_samps)
{
    if(new_num_samps < 1)
        throw std::invalid_argument("Invalid sample size, must be greater than 0");

    n_samp = new_num_samps;
    raw_sigs = std::make_shared<Mat>(Mat::Zero(n_sig, n_samp));
    mixed_sigs = std::make_shared<Mat>(Mat::Zero(mixing_mat->rows(), n_samp));
}

void Mixer::setMixingMatrix(Mat A)
{
    if(A.cols() != raw_sigs->rows())
        throw std::invalid_argument("Mixing matrix A's number of columns must match the number of signals to be mixed");

    mixing_mat = std::make_unique<Mat>(A);

    if(mixed_sigs->rows() != A.rows())
        mixed_sigs = std::make_shared<Mat>(Mat::Zero(A.rows(),n_samp));
}



void Mixer::genSignals(std::string& gen_file )
/* Textfile must be of form

number_of_signals number_of_samples
signal_amplitude_1 signal_period_1 signal_type_1
signal_amplitude_2 signal_period_2 signal_type_2

where:
    signal_type is a string
    number_of_signals and number_of_samples are integers > 0
    signal_amplitude and signal_period are double
*/
{
    std::ifstream instr_file (gen_file);
    if(!instr_file)
        throw std::invalid_argument("File with name \"" +gen_file+ "\" not found");

    std::string line,opt;
    int m,n,i{0};
    double amp, tau;

    std::getline (instr_file, line);
    std::stringstream ss(line);
    ss>>m>>n;

    if(ss.fail())
        throw std::invalid_argument("Invalid instruction in text file, line: " + std::to_string(i+2));

    setNumSignalsAndSamples(m,n);
    SigGen::WaveGen waves((int)n_samp);

    while (std::getline (instr_file, line)) {

        std::stringstream ss(line);
        ss >> amp >> tau >> opt;

        if(ss.fail())
            throw std::invalid_argument("Invalid instruction in text file, line: " + std::to_string(i+2));

        raw_sigs->row(i++) = waves.genWave(amp, tau, opt);

        if(i > n_sig)
            throw std::invalid_argument("Mismatch between number of signals (" +std::to_string(n_sig)+ 
                ") and amount of signal instructions in file (" +std::to_string(i)+")");
    }

    instr_file.close();

}

void Mixer::genSignals(int seed, int max_amp, int max_period)
{
    //srand( (int)(time(0)) );
    srand(seed);
    SigGen::WaveGen waves((int)n_samp);

    for(int i=0; i<n_sig; i++)
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

        Mat v = Mat::NullaryExpr(n_sig,n_samp, norm );
        v += *raw_sigs;

        *mixed_sigs = *mixing_mat *v;
    }
    else
        *mixed_sigs = *mixing_mat * *raw_sigs;
}
