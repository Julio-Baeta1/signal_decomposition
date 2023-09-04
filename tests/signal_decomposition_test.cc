#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../include/wave_gen.h"

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

TEST(EigenTest, BasicEqual) {
  Eigen::VectorXd x1 = Eigen::VectorXd::Zero(3);
  x1 = x1+Eigen::VectorXd::Constant(3,1);
  Eigen::VectorXd x2 = Eigen::VectorXd::Ones(3);
  // Expect equality.
  EXPECT_EQ(x1, x2);
}

TEST(Wave_gen_Test, DCWave) {
  Eigen::VectorXd x1 = Eigen::VectorXd::Constant(16,5) ;
  WaveGen wave_gen(16,0,15);
  Eigen::VectorXd x2 = wave_gen.genWave(5,0,"dc");
  // Expect equality.
  EXPECT_EQ(x1, x2);
}


