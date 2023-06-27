#include <gtest/gtest.h>

#include "../../cfd_mx/inc/matrix/csr_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/csr_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/ell_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab_restart.hpp"
#include "../validation_tool/l2norm_validation.hpp"
#include "../validation_tool/random_vector.hpp"
#include "io_2d_plot3d.hpp"
#include "poisson.hpp"


namespace{

class OpenMPMatrixTest : public ::testing::Test {
};

TEST(OpenMPMatrixTest, Bicgstab) {
  HeatConduction heat_conduction(10);

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

TEST(OpenMPMatrixTest, solver_bicgstabRe2) {
  HeatConduction heat_conduction(10);

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

  solverA(heat_conduction.rhs, resultA);
  solverB(heat_conduction.rhs, resultB);

  EXPECT_EQ(solverA(heat_conduction.rhs, resultA),
            solverB(heat_conduction.rhs, resultB));
  EXPECT_EQ(resultA, resultB);
}

TEST(OpenMPMatrixTest, get_csr) {
  HeatConduction heat_conduction(10);

  mat::EllMat<double> matA;
  matA.Set(heat_conduction.GetCsr());

  mat::CsrMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  auto rhs = GetRandomVectorDouble(heat_conduction.rhs.size());

  EXPECT_EQ(matA.row(), matB.row());
  EXPECT_EQ(matA.col(), matB.col());
  EXPECT_EQ(matA * rhs, matB * rhs);
}


} // namespace
