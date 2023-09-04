#include "wave_gen.h"

WaveGen::WaveGen(int length, double start, double end)
//Utilises built in linspace to build base signal for all signals
// size is the same for all signals generated from this Wave generator
{
    size = length;
    domain_sig = Eigen::VectorXd::LinSpaced(length,start,end);
}

Eigen::VectorXd WaveGen::dcWave(double amp)
{
        return Eigen::VectorXd::Constant(this->size,amp);
} 

 Eigen::VectorXd WaveGen::sinWave(double amp, double period)
{
        double freq = 1/period;
        Eigen::VectorXd wave = (2.0 * EIGEN_PI* freq *domain_sig);
        wave = wave.array().sin();
        return wave;
} 

 Eigen::VectorXd WaveGen::cosWave(double amp, double period)
{
        double freq = 1/period;
        Eigen::VectorXd wave = (2.0 * EIGEN_PI* freq *this->domain_sig);
        wave = wave.array().cos();
        return wave;
} 

 Eigen::VectorXd WaveGen::squareWave(double amp, double period)
{
        Eigen::VectorXd wave =  Eigen::VectorXd::Zero(size);
        int p = static_cast<int>( period ); //convert period to int to ensure modulus arthmetic and integer comparison is correct
        int half_p = static_cast<int>( period/2 ); //mid-point of square wave

        for(int i=0; i<size; i++)
                wave[i] = (i%p) < half_p ? 1:-1; //if in first half of square wave amp == 1, in second half and amp == -1

        return wave;
} 

 Eigen::VectorXd WaveGen::sawtoothWave(double amp, double period)
{
        Eigen::VectorXd wave =  Eigen::VectorXd::Zero(this->size);
        int p = static_cast<int>( period ); //convert period to int to ensure modulus arthmetic is correct

        for(int i=0; i<size; i++)
                wave[i] = amp*(i%p) / (p-1); //modulus by p and divide by p-1 to get range 0->1 over period which is scaled by amp

        return wave;
} 

 Eigen::VectorXd WaveGen::triangleWave(double amp, double period)
{
        Eigen::VectorXd wave =  Eigen::VectorXd::Zero(this->size);
        int p = static_cast<int>( period ); //convert period to int to ensure modulus arthmetic is correct

        for(int i=0; i<size; i++)
                wave[i] = amp * (abs(2*(i%p)/period - 1)-1); 
                //modulus by p and divide by period/2 to get range 0->~2. Subtract one and take absolute value to get range 0->1 with
                //minumum (0) in the middle and maximum (1) on left.

        return wave;
} 

Eigen::VectorXd WaveGen::genWave(double amplitude, double period, std::string opt)
//Select which wave to generate
{
    if (period == 0 || opt == "dc")
        return dcWave(amplitude);

    if (opt == "sin")
        return sinWave(amplitude,period);
    if (opt == "cos")
        return cosWave(amplitude,period);
    if (opt == "square")
        return squareWave(amplitude,period);
    if (opt == "sawtooth")
        return sawtoothWave(amplitude,period);
    if (opt == "triangle")
        return triangleWave(amplitude,period);

    return Eigen::VectorXd::Zero(this->size);
    
}