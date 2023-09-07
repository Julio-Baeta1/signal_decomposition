#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../include/wave_gen.h"


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




