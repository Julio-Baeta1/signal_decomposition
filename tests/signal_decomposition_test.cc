#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <filesystem>
#include "../include/wave_gen.h"
#include "../include/mixer.h"

using namespace SigGen;
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;

constexpr double MIN_PRECISION = 1.0E-6 ; //Used as tolerance for deciding if two matrices are equivalent
constexpr double MAX_PRECISION = 5.0E-1 ; //If matrices are not equivalent, we want the tolerance to be high for greater certainty

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
  Vec test_wave = gen_parm.genWave(1,0,"dc");
  Vec exp_wave = Vec::Ones(16);
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateSinWave) {
  WaveGen gen_parm(16,0);
  Vec test_wave = gen_parm.genWave(1,16,"sin");
  Vec exp_wave = (2.0 * EIGEN_PI* 1/16 *Vec::LinSpaced(16,0,16)).array().sin();
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateCosWave) {
  WaveGen gen_parm(16,0);
  Vec test_wave = gen_parm.genWave(1,16,"cos");
  Vec exp_wave = (2.0 * EIGEN_PI* 1/16 *Vec::LinSpaced(16,0,16)).array().cos();
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateSquareWave) {
  WaveGen gen_parm(16,0);
  Vec test_wave = gen_parm.genWave(1,4,"square");
  Vec exp_wave{{1,1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,1,-1,-1}};
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateSawtoothWave) {
  WaveGen gen_parm(16,0);
  Vec test_wave = gen_parm.genWave(1,5,"sawtooth");
  Vec exp_wave{{0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1,0}};
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateTriangleWave) {
  WaveGen gen_parm(16,0);
  Vec test_wave = gen_parm.genWave(3,6,"triangle");
  Vec exp_wave{{3, 1, -1, -3, -1, 1, 3, 1, -1, -3, -1, 1, 3, 1, -1, -3}};
  EXPECT_TRUE(test_wave.isApprox(exp_wave)); //Possibly Fix to use Expect_Double_True
}

TEST(WaveGenTest, GenerateZeroPeriodWave) {
  WaveGen gen_parm(16,0);
  Vec test_wave = gen_parm.genWave(1,0,"cos");
  Vec exp_wave = Vec::Ones(16);
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateWaveUsingInvalidOption) {
  WaveGen gen_parm(16,0);
  EXPECT_THROW(
    Vec test_wave = gen_parm.genWave(1,16,"s");
    ,std::invalid_argument);
}

TEST(WaveGenTest, GenerateWaveUsingInvalidOptionAndZeroPeriod) {
  WaveGen gen_parm(16,0);
  Vec test_wave = gen_parm.genWave(1,0,"s");
  Vec exp_wave = Vec::Ones(16);
  EXPECT_EQ(test_wave,exp_wave);
}

TEST(WaveGenTest, GenerateWaveUsingZeroSizeWaveGen) {
  WaveGen gen_parm;
  EXPECT_THROW(
    Vec test_wave = gen_parm.genWave(1,0,"cos");
    ,std::length_error);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//mixer

TEST(MixerTest, DefaultConstructor) {
  Mixer mixer_default;
  Mat zeros_2_by_8 = Mat::Zero(2, 8);

  EXPECT_EQ(mixer_default.getNumSignals(),(size_t)2);
  EXPECT_EQ(mixer_default.getNumSamples(),(size_t)8);
  EXPECT_EQ(mixer_default.getMixingMatrixValues(),Mat::Zero(2, 2));
  EXPECT_EQ(mixer_default.getRawSignalsValues(),zeros_2_by_8);
  EXPECT_EQ(mixer_default.getMixedSignalsValues(),zeros_2_by_8); 
}

TEST(MixerTest, ParamaterConstructor) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  Mat zeros_2_by_16 = Mat::Zero(2, 16);

  EXPECT_EQ(mixer_parm.getNumSignals(),(size_t)2);
  EXPECT_EQ(mixer_parm.getNumSamples(),(size_t)16);
  EXPECT_EQ(mixer_parm.getMixingMatrixValues(),A);
  EXPECT_EQ(mixer_parm.getRawSignalsValues(),zeros_2_by_16);
  EXPECT_EQ(mixer_parm.getMixedSignalsValues(),zeros_2_by_16);
}

TEST(MixerTest, SetNumSignalsAndSamples) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  mixer_parm.setNumSignalsAndSamples(3,3);

  EXPECT_EQ(mixer_parm.getNumSignals(),(size_t)3);
  EXPECT_EQ(mixer_parm.getNumSamples(),(size_t)3);
  EXPECT_EQ(mixer_parm.getMixingMatrixValues(),Mat::Zero(2, 3));
  EXPECT_EQ(mixer_parm.getRawSignalsValues(),Mat::Zero(3, 3));
  EXPECT_EQ(mixer_parm.getMixedSignalsValues(),Mat::Zero(2, 3));
}

TEST(MixerTest, SetNumSignalsInvalidNumSignalsSize) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  
  EXPECT_THROW(
    mixer_parm.setNumSignalsAndSamples(-1,2);
    ,std::invalid_argument);
}

TEST(MixerTest, SetNumSignalsInvalidNumSamplesSize) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  
  EXPECT_THROW(
    mixer_parm.setNumSignalsAndSamples(2,-5);
    ,std::invalid_argument);
}

TEST(MixerTest, SetNumSignals) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  mixer_parm.setNumSignals(3);

  EXPECT_EQ(mixer_parm.getNumSignals(),(size_t)3);
  EXPECT_EQ(mixer_parm.getNumSamples(),(size_t)16);
  EXPECT_EQ(mixer_parm.getMixingMatrixValues(),Mat::Zero(2, 3));
  EXPECT_EQ(mixer_parm.getRawSignalsValues(),Mat::Zero(3, 16));
  EXPECT_EQ(mixer_parm.getMixedSignalsValues(),Mat::Zero(2, 16));
}

TEST(MixerTest, SetNumSignalsInvalidSize) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  
  EXPECT_THROW(
    mixer_parm.setNumSignals(0);
    ,std::invalid_argument);
}

TEST(MixerTest, SetNumSamples) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);
  mixer_parm.setNumSamples(5);

  EXPECT_EQ(mixer_parm.getNumSignals(),(size_t)2);
  EXPECT_EQ(mixer_parm.getNumSamples(),(size_t)5);
  EXPECT_EQ(mixer_parm.getMixingMatrixValues(),A);
  EXPECT_EQ(mixer_parm.getRawSignalsValues(),Mat::Zero(2, 5));
  EXPECT_EQ(mixer_parm.getMixedSignalsValues(),Mat::Zero(2, 5));
}

TEST(MixerTest, SetNumSamplesInvalidSize) {
  Mat A{{1,2},{3,0.05}};
  Mixer mixer_parm(2,16,A);

  EXPECT_THROW(
    mixer_parm.setNumSamples(-5);
    ,std::invalid_argument);

}

TEST(MixerTest, SetMixingMatrix) {
  Mixer mixer_default;
  Mat A{{1,2},{3,0.05}};
  mixer_default.setMixingMatrix(A);

  EXPECT_EQ(mixer_default.getMixingMatrixValues(),A);
}

TEST(MixerTest, SetMixingMatrixIncorrectDims) {
  Mixer mixer_default;
  Mat A{{1,2,3},{4,5,6},{7,8,9}};
  
  EXPECT_THROW(
    mixer_default.setMixingMatrix(A);
    ,std::invalid_argument);
}

TEST(MixerTest, GenerateRawSignalsFromFileInstructions) {
  Mixer mixer_parm(4,4);
  std::string test_file = "../test_mixer_sig_gen_from_file.txt";
  mixer_parm.genSignals(test_file);

  Mat* gen_sigs = mixer_parm.getRawSignalsSharedPtr().get();
  Mat exp_sigs{{0, 0.166769, 0.328867, 0.481754},
                {1, 1, 1, -1},
                {2, 1.338262, -0.209056, -1.618034},
                {0, 0.0238095, 0.047619, 0.0714286}};

  EXPECT_TRUE(exp_sigs.isApprox(*gen_sigs,MIN_PRECISION)); //second arg (p) percision tolerance ∥v−w∥ ⩽ p min(∥v∥,∥w∥)
}

TEST(MixerTest, GenerateRawSignalsFromFileInstructionsFileNotFound) {
  Mixer mixer_parm(4,4);
  std::string test_file = "../test_file.txt";
  
  EXPECT_THROW(
    mixer_parm.genSignals(test_file);
    ,std::invalid_argument);
}

TEST(MixerTest, GenerateRawSignalsFromFileInstructionsMissingArgHeader) {
  Mixer mixer_parm(4,4);
  std::string test_file = "../test_mixer_sig_gen_from_file_missing_arg_header.txt";
  
  EXPECT_THROW(
    mixer_parm.genSignals(test_file);
    ,std::invalid_argument);
}

TEST(MixerTest, GenerateRawSignalsFromFileInstructionsMissingArgInstruction) {
  Mixer mixer_parm(4,4);
  std::string test_file = "../test_mixer_sig_gen_from_file_missing_arg_instruction.txt";
  
  EXPECT_THROW(
    mixer_parm.genSignals(test_file);
    ,std::invalid_argument);
}

TEST(MixerTest, GenerateRawSignalsFromFileInstructionsFileTooManySignals) {
  Mixer mixer_parm(4,4);
  std::string test_file = "../test_mixer_sig_gen_from_file_too_many_sigs";
  
  EXPECT_THROW(
    mixer_parm.genSignals(test_file);
    ,std::invalid_argument);
}

TEST(MixerTest, GenerateRawSignals) {
  //Tests if signal is generated not quality of signals
  Mixer mixer_parm(4,4);
  Mat zero_sigs = Mat::Zero(4,4);

  mixer_parm.genSignals(6);
  Mat* gen_sigs = mixer_parm.getRawSignalsSharedPtr().get();

  EXPECT_FALSE(zero_sigs.isApprox(*gen_sigs,MAX_PRECISION)); //second arg (p) percision tolerance ∥v−w∥ ⩽ p min(∥v∥,∥w∥)
}

TEST(MixerTest, MixSignalsNoNoise) {
  Mat A{{0.1,0.3,0.15,0.95},{0.4,0.5,0.6,0.5},{0.07,1.2,0.01,0.9},{0.01,0.85,0.15,0.23}};
  Mixer mixer_parm(4,4);

  std::string test_file = "../test_mixer_sig_gen_from_file.txt";
  mixer_parm.genSignals(test_file);
  mixer_parm.setMixingMatrix(A);

  mixer_parm.mixSignals(false,123);
  Mat* mixed = mixer_parm.getMixedSignalsSharedPtr().get();

  Mat exp_mixed{{0, 0.166769, 0.328867, 0.481754},
                {1, 1, 1, -1},
                {2, 1.338262, -0.209056, -1.618034},
                {0, 0.0238095, 0.047619, 0.0714286}};
  exp_mixed.rowwise().normalize();

  exp_mixed = A*exp_mixed; 

  EXPECT_TRUE(exp_mixed.isApprox(*mixed,MIN_PRECISION)); //second arg (p) percision tolerance ∥v−w∥ ⩽ p min(∥v∥,∥w∥)
}

TEST(MixerTest, MixSignalsWithNoise) {
  Mat A{{0.1,0.3,0.15,0.95},{0.4,0.5,0.6,0.5},{0.07,1.2,0.01,0.9},{0.01,0.85,0.15,0.23}};
  Mixer mixer_parm(4,4);

  std::string test_file = "../test_mixer_sig_gen_from_file.txt";
  mixer_parm.genSignals(test_file);
  mixer_parm.setMixingMatrix(A);

  mixer_parm.mixSignals(true,123);
  Mat* mixed = mixer_parm.getMixedSignalsSharedPtr().get();

  Mat exp_mixed{{0, 0.166769, 0.328867, 0.481754},
                {1, 1, 1, -1},
                {2, 1.338262, -0.209056, -1.618034},
                {0, 0.0238095, 0.047619, 0.0714286}};
  exp_mixed.rowwise().normalize();

  exp_mixed = A*exp_mixed; 

  EXPECT_FALSE(exp_mixed.isApprox(*mixed,MAX_PRECISION)); //second arg (p) percision tolerance ∥v−w∥ ⩽ p min(∥v∥,∥w∥)
}