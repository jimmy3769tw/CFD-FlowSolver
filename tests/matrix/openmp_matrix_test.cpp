#include <gtest/gtest.h>

#include "../../cfd_mx/inc/matrix/csr_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/ell_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab_restart.hpp"
#include "../validation_tool/l2norm_validation.hpp"
#include "io_2d_plot3d.hpp"
#include "poisson.hpp"
HeatConduction heat_conduction(10);

vector<int>
RandomVector(int size, int range_begin = 0) {
  vector<int> result(size);
  for (auto& x : result) {
    x = rand() % size + range_begin;
  }
  return result;
}

vector<double> RandomVectorDouble(int size, int range_begin = 0) {
  vector<double> result(size);
  for (auto& x : result) {
    x = rand() % size + range_begin + rand()/size;
  }
  return result;
}



TEST(csr_mat, resize_row_col) {
  mat::CsrMat<double> matA(5, 4);
  EXPECT_EQ(matA.row(), 5);
  EXPECT_EQ(matA.col(), 4);
  matA.resize(4, 3);
  EXPECT_EQ(matA.row(), 4);
  EXPECT_EQ(matA.col(), 3);
}

TEST(csr_mat, test2) {
  mat::CsrMat<double> matA;
  matA.Set(heat_conduction.GetCsr());
  EXPECT_EQ(heat_conduction.GetCsr(), matA.GetCsr());
}

TEST(csr_ell_mat, solver_bicgstabRe2) {
  mat::CsrMat<double> matA;
  matA.Set(heat_conduction.GetCsr());
  solver::BicgstabRestart<mat::CsrMat<double>> solverA(matA);

  mat::EllMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  solver::BicgstabRestart<mat::EllMat<double>> solverB(matB);

  auto resultA = heat_conduction.x;
  auto resultB = heat_conduction.x;
  solverA.setTolerance(0.0001);
  solverB.setTolerance(0.0001);
  EXPECT_EQ(solverA(heat_conduction.rhs, resultA),
            solverB(heat_conduction.rhs, resultB));
  EXPECT_EQ(resultA, resultB);
}

TEST(csr_ell_mat, mult) {
  mat::EllMat<double> matA;
  matA.Set(heat_conduction.GetCsr());

  mat::CsrMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  auto rhs = RandomVectorDouble(heat_conduction.rhs.size());

  EXPECT_EQ(matA.row(), matB.row());
  EXPECT_EQ(matA.col(), matB.col());
  EXPECT_EQ(matA * rhs, matB * rhs);
}

TEST(bicgstab, solver_2) {
  mat::CsrMat<double> matA;
  matA.Set(heat_conduction.GetCsr());
  solver::Bicgstab<mat::CsrMat<double>> solverA(matA);

  mat::EllMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  solver::Bicgstab<mat::EllMat<double>> solverB(matB);

  auto resultA = heat_conduction.x;
  auto resultB = heat_conduction.x;
  solverA.setTolerance(0.0001);
  solverB.setTolerance(0.0001);
  solverA(heat_conduction.rhs, resultA);
  solverB(heat_conduction.rhs, resultB);
  EXPECT_TRUE(L2Norm(resultA, resultB) < 0.0001);
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
