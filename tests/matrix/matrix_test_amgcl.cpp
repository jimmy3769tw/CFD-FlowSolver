#include <gtest/gtest.h>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>

#include "../../cfd_mx/inc/matrix/csr_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/ell_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab_restart.hpp"
#include "../../cfd_mx/inc/mpi_tool/mpi_complex.hpp"
#include "../../cfd_mx/inc/mpi_tool/mpi_getter.hpp"
#include "../../cfd_mx/inc/mpi_tool/mpi_tool.hpp"
#include "../validation_tool/l2norm_validation.hpp"
#include "poisson.hpp"

HeatConduction heat_conduction(10);

//   the solver backend:
typedef amgcl::backend::builtin<double> SBackend;
//   the preconditioner backend:
typedef amgcl::backend::builtin<double> PBackend;

//=========== Compose the solver type ===========//

typedef amgcl::make_solver<
    amgcl::amg<PBackend, amgcl::coarsening::smoothed_aggregation,
               amgcl::relaxation::gauss_seidel>,
    amgcl::solver::bicgstab<SBackend> >
    Solver;


TEST(csr_mat, get_csr) {
  auto data = heat_conduction.GetCsr();
  int row = heat_conduction.rhs.size();
  std::vector<double> amgcl_x = heat_conduction.x;
  std::vector<double> resultB = heat_conduction.x;
  auto ptr = std::get<0>(data);      // size_t
  auto indices = std::get<1>(data);  // szie_t
  auto values = std::get<2>(data);   // double
  
  auto amgcl_mat = std::tie(row, ptr, indices, values);
  Solver solve(amgcl_mat);
  auto [iters, error] = solve(amgcl_mat, heat_conduction.rhs, amgcl_x);

  mat::CsrMat<double> matB;
  matB.Set(heat_conduction.GetCsr());
  solver::BicgstabRestart<mat::CsrMat<double>> solverB(matB);
  solverB.SetTolerance(1.0e-9);

  solverB(heat_conduction.rhs, resultB);

  EXPECT_TRUE(L2Norm(amgcl_x, resultB) < 1.0e-6);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
