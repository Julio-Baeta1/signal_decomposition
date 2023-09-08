#include "wave_gen.h"

namespace SigGen{

WaveGen::WaveGen(int length, double start)
//Utilises built in linspace to build base signal for all signals
// size is the same for all signals generated from this Wave generator
{
    size = length;
    domain_sig = Eigen::VectorXd::LinSpaced(length,start,start+length);
}

void WaveGen::setDomain(int new_size, double new_start)
{
        if(size == 0)
        {
                size = new_size;
                domain_sig = Eigen::VectorXd::LinSpaced(size,new_start,new_start+size);
        }
        else
                throw std::invalid_argument("The domain of WaveGen object is already set");
           
}

Eigen::VectorXd WaveGen::dcWave(double amp)
{
        return Eigen::VectorXd::Constant(this->size,amp);
} 

 Eigen::VectorXd WaveGen::sinWave(double amp, double period)
{
        double freq = 1/period;
        Eigen::VectorXd wave = (2.0 * EIGEN_PI* freq *this->domain_sig);
        wave = wave.array().sin();
        return wave;
} 

 Eigen::VectorXd WaveGen::cosWave(double amp, double period)
{
        double freq = 1/period;
        Eigen::VectorXd wave = (2.0 * EIGEN_PI* freq * this->domain_sig);
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
                wave[i] = amp* (2*abs(2*(i%p)/period - 1)-1); 
                //modulus by p, divide by period/2 and subtract 1 to get range -1->0 with 0 occuring in the middle. Take absolute value
                //so the range is know 0->1. Multiply by 2 and subtract 1 so it is from -1->1, which can be scaled by the amplitude

        return wave;
} 

Eigen::VectorXd WaveGen::genWave(double amplitude, double period, std::string opt)
//Select which wave to generate
{
    if(this->size == 0)
        throw std::length_error("Wavegen object of size zero");

    
    if (period == 0)
        opt = "dc";

    if (opt == "dc")
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
    

    throw std::invalid_argument("invalid option \""+opt+"\" entered");
    
}

Eigen::VectorXd WaveGen::genRandomWave(int max_amp, int max_period)
//Select which wave to generate
{
    if(this->size == 0)
        throw std::length_error("Wavegen object of size zero");

    SignalType rand_wave = static_cast<SignalType>(  rand() % static_cast<int>(SignalType::size) );
    int rand_period = (rand()%max_period)+1; 
    int rand_amp = (rand()%max_amp+1)+1;
    
    switch (rand_wave){ 
    
        case SignalType::dc:
                return dcWave(rand_amp);
        case SignalType::sin:
                return sinWave(rand_amp,rand_period);
        case SignalType::cos:
                return cosWave(rand_amp,rand_period);
        case SignalType::square:
                return squareWave(rand_amp,rand_period);
        case SignalType::sawtooth:
                return sawtoothWave(rand_amp,rand_period);
        case SignalType::triangle:
                return triangleWave(rand_amp,rand_period);

        default:
                throw std::invalid_argument("Random signal generator selected a none SignalType option");
    }
    
}

}//End SigGen namespace