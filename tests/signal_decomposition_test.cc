#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../include/wave_gen.h"
#include "../include/mixer.h"

using namespace SigGen;

constexpr double MIN_PRECISION = 1.0E-6 ; //Used as tolerance for deciding if two matrices are equivalent
constexpr double MAX_PRECISION = 1.0 ; //If matrices are not equivalent, we want the tolerance to be high for greater certainty

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//wave_gen

TEST(WaveGenTest, DefaultConstructor) {
  WaveGen gen_default;
  size_t len = gen_default.len();
  EXPECT_EQ(len,0);
}

TEST(WaveGenTest, ParameterConstructor) {
  WaveGen gen_parm(16,0);
  size_t len = gen_parm.len();
  EXPECT_EQ(len,16);
}

TEST(WaveGenTest, SetDomainEmptyWaveGen) {
  WaveGen gen_default;
  gen_default.setDomain(16,0);
  size_t len = gen_default.len();
  EXPECT_EQ(len,16);
}

TEST(WaveGenTest, SetDomainFullWaveGen) {
  WaveGen gen_parm(16,0);
  EXPECT_THROW(gen_parm.setDomain(16,0),std::invalid_argument);
}

TEST(WaveGenTest, GenerateDCWave) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(1,0,"dc");
  Eigen::VectorXd exp_wave = Eigen::VectorXd::Ones(16);
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateSinWave) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(1,16,"sin");
  Eigen::VectorXd exp_wave = (2.0 * EIGEN_PI* 1/16 *Eigen::VectorXd::LinSpaced(16,0,16)).array().sin();
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateCosWave) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(1,16,"cos");
  Eigen::VectorXd exp_wave = (2.0 * EIGEN_PI* 1/16 *Eigen::VectorXd::LinSpaced(16,0,16)).array().cos();
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateSquareWave) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(1,4,"square");
  Eigen::VectorXd exp_wave{{1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1}};
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateSawtoothWave) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(1,5,"sawtooth");
  Eigen::VectorXd exp_wave{{0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1,0}};
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateTriangleWave) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(3,6,"triangle");
  Eigen::VectorXd exp_wave{{3, 1, -1, -3, -1, 1, 3, 1, -1, -3, -1, 1, 3, 1, -1, -3}};
  EXPECT_TRUE(test_wave.isApprox(exp_wave)); //Possibly Fix to use Expect_Double_True
}

TEST(WaveGenTest, GenerateZeroPeriodWave) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(1,0,"cos");
  Eigen::VectorXd exp_wave = Eigen::VectorXd::Ones(16);
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateWaveUsingInvalidOption) {
  WaveGen gen_parm(16,0);
  EXPECT_THROW(
    Eigen::VectorXd test_wave = gen_parm.genWave(1,16,"s");
    ,std::invalid_argument);
}

TEST(WaveGenTest, GenerateWaveUsingInvalidOptionAndZeroPeriod) {
  WaveGen gen_parm(16,0);
  Eigen::VectorXd test_wave = gen_parm.genWave(1,0,"s");
  Eigen::VectorXd exp_wave = Eigen::VectorXd::Ones(16);
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateWaveUsingZeroSizeWaveGen) {
  WaveGen gen_parm;
  EXPECT_THROW(
    Eigen::VectorXd test_wave = gen_parm.genWave(1,0,"cos");
    ,std::length_error);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//mixer

TEST(MixerTest, DefaultConstructor) {
  Mixer mixer_default;
  Eigen::MatrixXd zeros_2_by_8 = Eigen::MatrixXd::Zero(2, 8);

  size_t d_num_sigs = mixer_default.getNumSignals();
  size_t d_num_sams = mixer_default.getNumSamples();
  Eigen::MatrixXd d_mixing = mixer_default.getMixingMatrixValues();
  Eigen::MatrixXd d_raw = mixer_default.getRawSignalsValues();
  Eigen::MatrixXd d_mixed = mixer_default.getMixedSignalsValues();

  EXPECT_EQ(d_num_sigs,(size_t)2);
  EXPECT_EQ(d_num_sams,(size_t)8);
  EXPECT_EQ(d_mixing,Eigen::MatrixXd::Zero(2, 2));
  EXPECT_EQ(d_raw,zeros_2_by_8);
  EXPECT_EQ(d_mixed,zeros_2_by_8);
  
}

TEST(MixerTest, ParamaterConstructor) {
  Eigen::MatrixXd A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  Eigen::MatrixXd zeros_2_by_16 = Eigen::MatrixXd::Zero(2, 16);

  size_t d_num_sigs = mixer_parm.getNumSignals();
  size_t d_num_sams = mixer_parm.getNumSamples();
  Eigen::MatrixXd d_mixing = mixer_parm.getMixingMatrixValues();
  Eigen::MatrixXd d_raw = mixer_parm.getRawSignalsValues();
  Eigen::MatrixXd d_mixed = mixer_parm.getMixedSignalsValues();

  EXPECT_EQ(d_num_sigs,(size_t)2);
  EXPECT_EQ(d_num_sams,(size_t)16);
  EXPECT_EQ(d_mixing,A);
  EXPECT_EQ(d_raw,zeros_2_by_16);
  EXPECT_EQ(d_mixed,zeros_2_by_16);
}

TEST(MixerTest, SetMixingMatrix) {
  Mixer mixer_default;
  Eigen::MatrixXd A{{1,2},{3,0.05}};
  mixer_default.setMixingMatrix(A);
  Eigen::MatrixXd d_mixing = mixer_default.getMixingMatrixValues();

  EXPECT_EQ(d_mixing,A);
}

TEST(MixerTest, SetMixingMatrixIncorrectDims) {
  Mixer mixer_default;
  Eigen::MatrixXd A{{1,2,3},{4,5,6},{7,8,9}};
  
  EXPECT_THROW(
    mixer_default.setMixingMatrix(A);
    ,std::invalid_argument);
}

TEST(MixerTest, SetNumSignals) {
  Eigen::MatrixXd A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  mixer_parm.setNumSignals(3);

  size_t d_num_sigs = mixer_parm.getNumSignals();
  size_t d_num_sams = mixer_parm.getNumSamples();
  Eigen::MatrixXd d_mixing = mixer_parm.getMixingMatrixValues();
  Eigen::MatrixXd d_raw = mixer_parm.getRawSignalsValues();
  Eigen::MatrixXd d_mixed = mixer_parm.getMixedSignalsValues();

  EXPECT_EQ(d_num_sigs,(size_t)3);
  EXPECT_EQ(d_num_sams,(size_t)16);
  EXPECT_EQ(d_mixing,Eigen::MatrixXd::Zero(2, 3));
  EXPECT_EQ(d_raw,Eigen::MatrixXd::Zero(3, 16));
  EXPECT_EQ(d_mixed,Eigen::MatrixXd::Zero(2, 16));
}

TEST(MixerTest, SetNumSamples) {
  Eigen::MatrixXd A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  mixer_parm.setNumSamples(5);

  size_t d_num_sigs = mixer_parm.getNumSignals();
  size_t d_num_sams = mixer_parm.getNumSamples();
  Eigen::MatrixXd d_mixing = mixer_parm.getMixingMatrixValues();
  Eigen::MatrixXd d_raw = mixer_parm.getRawSignalsValues();
  Eigen::MatrixXd d_mixed = mixer_parm.getMixedSignalsValues();

  EXPECT_EQ(d_num_sigs,(size_t)2);
  EXPECT_EQ(d_num_sams,(size_t)5);
  EXPECT_EQ(d_mixing,A);
  EXPECT_EQ(d_raw,Eigen::MatrixXd::Zero(2, 5));
  EXPECT_EQ(d_mixed,Eigen::MatrixXd::Zero(2, 5));
}

TEST(MixerTest, SetNumSignalsInvalidSize) {
  Eigen::MatrixXd A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  
  EXPECT_THROW(
    mixer_parm.setNumSignals(0);
    ,std::invalid_argument);

}

TEST(MixerTest, SetNumSamplesInvalidSize) {
  Eigen::MatrixXd A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);

  EXPECT_THROW(
    mixer_parm.setNumSamples(-5);
    ,std::invalid_argument);

}

TEST(MixerTest, GenerateRawSignals) {
  Mixer mixer_parm(4,4);

  mixer_parm.genSignals(6);
  Eigen::MatrixXd* gen_sigs = mixer_parm.getRawSignalsSharedPtr().get();
  //Eigen::MatrixXd gen_sigs = mixer_parm.getRawSignalsValues();
  Eigen::MatrixXd exp_sigs{{0, 0.0410959, 0.0821918,  0.123288},
                           {0,  0.930874,  0.680173, -0.433884},
                           {0, 0.0909091,  0.181818,  0.272727},
                           {1,  0.992193,  0.968893,  0.930465}};

  EXPECT_TRUE(exp_sigs.isApprox(*gen_sigs,MIN_PRECISION)); //second arg (p) percision tolerance ∥v−w∥ ⩽ p min(∥v∥,∥w∥)
}

TEST(MixerTest, MixSignalsNoNoise) {
  Eigen::MatrixXd A{{0.1,0.3,0.15,0.95},{0.4,0.5,0.6,0.5},{0.07,1.2,0.01,0.9},{0.01,0.85,0.15,0.23}};
  Mixer mixer_parm(4,4,A);

  mixer_parm.genSignals(6);
  mixer_parm.mixSignals(false,123);
  Eigen::MatrixXd* mixed = mixer_parm.getMixedSignalsSharedPtr().get();

  Eigen::MatrixXd exp_mixed{{0, 0.0410959, 0.0821918,  0.123288},
                           {0,  0.930874,  0.680173, -0.433884},
                           {0, 0.0909091,  0.181818,  0.272727},
                           {1,  0.992193,  0.968893,  0.930465}};
  exp_mixed = A*exp_mixed; 

  EXPECT_TRUE(exp_mixed.isApprox(*mixed,MIN_PRECISION)); //second arg (p) percision tolerance ∥v−w∥ ⩽ p min(∥v∥,∥w∥)
}

TEST(MixerTest, MixSignalsWithNoise) {
  Eigen::MatrixXd A{{0.1,0.3,0.15,0.95},{0.4,0.5,0.6,0.5},{0.07,1.2,0.01,0.9},{0.01,0.85,0.15,0.23}};
  Mixer mixer_parm(4,4,A);

  mixer_parm.genSignals(6);
  mixer_parm.mixSignals(true,123);
  Eigen::MatrixXd* mixed = mixer_parm.getMixedSignalsSharedPtr().get();

  Eigen::MatrixXd exp_mixed{{0, 0.0410959, 0.0821918,  0.123288},
                           {0,  0.930874,  0.680173, -0.433884},
                           {0, 0.0909091,  0.181818,  0.272727},
                           {1,  0.992193,  0.968893,  0.930465}};
  exp_mixed = A*exp_mixed; 

  EXPECT_FALSE(exp_mixed.isApprox(*mixed,MAX_PRECISION)); //second arg (p) percision tolerance ∥v−w∥ ⩽ p min(∥v∥,∥w∥)
}