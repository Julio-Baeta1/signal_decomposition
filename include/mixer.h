#ifndef MIXER_H
#define MIXER_H

#include"wave_gen.h"

namespace Mixing{

enum class SignalType{
    dc=1,sin,cos,square,sawtooth,triangle
};

constexpr const char* SignalTypeToString(SignalType st) noexcept
{
    switch (st)
    {
        case SignalType::dc: return "dc";
        case SignalType::sin: return "sin";
        case SignalType::cos: return "cos";
        case SignalType::square: return "square";
        case SignalType::sawtooth: return "sawtooth";
        case SignalType::triangle: return "triangle";
        //No default to throw error
    }
}

class Mixer{
    private:
        Eigen::MatrixXd raw_sigs,mixed_sigs,mixing_mat; //Will make pointers
        size_t num_sigs, num_samples;
        bool noisy;

    public:
        Mixer(int num_signals=8, int signal_duration=8, Eigen::MatrixXd mixing_matrix, bool is_noisy=false);

        Eigen::MatrixXd getRawSignals() const {return raw_sigs;}
        Eigen::MatrixXd getMixedSignals() const {return mixed_sigs;}
        Eigen::MatrixXd getMixingMatrix() const {return mixing_mat;}
        bool isNoisy() const {return noisy;}

        void genSignals();
        void setMixingMatrix();
        void mixSignals();
};

}
#endif