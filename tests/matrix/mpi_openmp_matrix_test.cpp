#include <gtest/gtest.h>

#include "../../cfd_mx/inc/matrix/csr_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/ell_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab_restart.hpp"
#include "../../cfd_mx/inc/matrix/solver/mpi/bicgstab_mpi.hpp"
#include "../../cfd_mx/inc/matrix/solver/mpi/bicgstab_restart_mpi.hpp"
#include "../../cfd_mx/inc/mpi_tool/mpi_complex.hpp"
#include "../../cfd_mx/inc/mpi_tool/mpi_getter.hpp"
#include "../../cfd_mx/inc/mpi_tool/mpi_tool.hpp"
#include "../validation_tool/l2norm_validation.hpp"
#include "io_2d_plot3d.hpp"
#include "poisson.hpp"

mpi::MpiComplex mpi_complex_;

HeatConduction heat_conduction(100);

TEST(csr_mpi, BicgstabRestart) {

  mat::CsrMat<double> matA;
  matA.Set(heat_conduction.GetCsr());
  matA.InitMpi(mpi_complex_.GetCommWorld(), matA.col());
  solver::BicgstabRestartMpi<mat::CsrMat<double>> solverA(matA);

  mat::CsrMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  solver::BicgstabRestart<mat::CsrMat<double>> solverB(matB);
  solverA.SetTolerance(1.0e-9);
  solverB.SetTolerance(1.0e-9);

  for (int i = 0; i < 100; i++) {
    auto resultA = heat_conduction.x;
    auto resultB = heat_conduction.x;

    solverA(heat_conduction.rhs, resultA);
    solverB(heat_conduction.rhs, resultB);
    EXPECT_TRUE(L2Norm(resultA, resultB) < 1.0e-8);
  }
}

TEST(csr_mpi, Bicgstab) {
  mat::CsrMat<double> matA;
  matA.Set(heat_conduction.GetCsr());
  matA.InitMpi(mpi_complex_.GetCommWorld(), matA.col());
  solver::BicgstabMpi<mat::CsrMat<double>> solverA(matA);

  mat::CsrMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  solver::Bicgstab<mat::CsrMat<double>> solverB(matB);
  solverA.SetTolerance(1.0e-9);
  solverB.SetTolerance(1.0e-9);

  for (int i = 0; i < 100; i++) {
    auto resultA = heat_conduction.x;
    auto resultB = heat_conduction.x;

    solverA(heat_conduction.rhs, resultA);
    solverB(heat_conduction.rhs, resultB);
    EXPECT_TRUE(L2Norm(resultA, resultB) < 1.0e-8);
  }
}

TEST(EllMat, BicgstabRestartMpi) {
  mat::EllMat<double> matA;
  matA.Set(heat_conduction.GetCsr());
  matA.InitMpi(mpi_complex_.GetCommWorld(), matA.col());
  solver::BicgstabRestartMpi<mat::EllMat<double>> solverA(matA);

  mat::CsrMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  solver::BicgstabRestart<mat::CsrMat<double>> solverB(matB);
  solverA.SetTolerance(1.0e-9);
  solverB.SetTolerance(1.0e-9);

  for (int i = 0; i < 100; i++) {
    auto resultA = heat_conduction.x;
    auto resultB = heat_conduction.x;

    solverA(heat_conduction.rhs, resultA);
    solverB(heat_conduction.rhs, resultB);
    EXPECT_TRUE(L2Norm(resultA, resultB) < 1.0e-8);
  }
}

TEST(EllMat, Bicgstab) {
  mat::EllMat<double> matA;
  matA.Set(heat_conduction.GetCsr());
  matA.InitMpi(mpi_complex_.GetCommWorld(), matA.col());
  solver::BicgstabMpi<mat::EllMat<double>> solverA(matA);

  mat::CsrMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  solver::Bicgstab<mat::CsrMat<double>> solverB(matB);
  solverA.SetTolerance(1.0e-9);
  solverB.SetTolerance(1.0e-9);

  for (int i = 0; i < 100; i++) {
    auto resultA = heat_conduction.x;
    auto resultB = heat_conduction.x;

    solverA(heat_conduction.rhs, resultA);
    solverB(heat_conduction.rhs, resultB);
    EXPECT_TRUE(L2Norm(resultA, resultB) < 1.0e-8);
  }
}

int main(int argc, char **argv) {
  mpi_complex_ = mpi::MpiComplex(argc, argv);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
