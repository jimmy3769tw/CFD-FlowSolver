#pragma once
#include <algorithm>

#include "backend/domain.hpp"
#include "backend/physical_variables.hpp"
#include "backend/print_openmp.hpp"
#include "backend/simulation.hpp"
#include "dfib/update_u_and_force.hpp"
#include "grid/structured_grid.hpp"
#include "io/csv/csv_structured_grid.hpp"
#include "io/io_tools.hpp"
#include "io/plot3d/write_qfile.hpp"
#include "io/plot3d/write_xfile.hpp"
#include "source/convection_and_difussion.hpp"
#include "boundary_condition/boundary_condition.hpp"
#include "boundary_condition/boundary_condition_copy.hpp"
#include "matrix/solver/bicgstab_restart.hpp"
#include "matrix/solver/bicgstab.hpp"
#include "pressure/pressure_mat.hpp"
#include "omp.h"
#include "analyses/central_profile.hpp"
#include "source/smagorinsky_model.hpp"

using solverType = solver::BicgstabRestart<MatType>;

namespace projection_method {
class CpuOpenMp{
 public:
  CpuOpenMp(Simulation &simu, StructuredGrid &grid)
      : simu_(simu),
        grid_(&grid),
        vel_(StaggeredVelocity(grid)),
        intermediate_vel_(StaggeredVelocity(grid)),
        dfib_(ImmersedBoundary(grid)),
        pressure_(Pressure(grid)),
        pressure_mat_(PressureMat(grid)),
        global_domain_(grid) {
    vel_.FillVel(simu.ini_condition.u, simu.ini_condition.v, simu.ini_condition.w);
    CreatOutputFile();
    PrintIsOpenmpExist();
    std::vector<int> grid_size{grid.nx, grid.ny, grid.nz};
    global_domain_.Init({1, 1, 1}, {0, 0, 0});
    int ompThreads = omp_get_max_threads();
    omp_set_num_threads(ompThreads);
    vel_ = intermediate_vel_ = vel_;
    csv::WriteCsvFile("Information/gA.csv", *grid_);
    plot3d::write_xfile(grid);
    UpdateAllVelocityOnBoundary(global_domain_, vel_, pressure_, grid);
    CopyVelocityOnBoundary(global_domain_, vel_, intermediate_vel_, grid);
  }

  void solve(){
    pressure_solver_ = solverType(pressure_mat_.mat_a);
    pressure_solver_.SetTolerance(simu_.poisson_criteria);
    CreatEta();
    WriteQfile(dfib_, simu_, pressure_, vel_, *grid_);
    for (; !simu_.tva.IsLoopFinish(); simu_.tva.AddLoop()) {
      CalProjectionMethod();
    }
  }

private :
  Simulation simu_;
  StructuredGrid* grid_;
  StaggeredVelocity vel_;
  StaggeredVelocity intermediate_vel_;
  Pressure pressure_;
  LocalDomain global_domain_;
  ImmersedBoundary dfib_;
  PressureMat pressure_mat_;
  solverType pressure_solver_;

  void CalProjectionMethod() {

    // CalSmagorinskyModel(simu_, vel_, global_domain_, *grid_);

    CalConvectionAndDiffusion(simu_, vel_, intermediate_vel_, global_domain_, *grid_);

    pressure_mat_.CalMatB(intermediate_vel_, simu_.tva.GetDt(), global_domain_);
    auto [iters, error] =
        pressure_solver_(pressure_mat_.mat_b, pressure_mat_.x_result);
    std::cout << "[iters, error]\t" << iters << '\t' << error << std::endl;
    pressure_mat_.ConvertResultToPressure(pressure_, global_domain_);

    UpdateForceAndVelocity(dfib_, simu_.tva.GetDt(), pressure_, intermediate_vel_, vel_,
                           global_domain_, *grid_);


    UpdateAllVelocityOnBoundary(global_domain_, vel_, pressure_, *grid_);
    // BC_staggered_copy(global_domain_, vel_, intermediate_vel_, grid_);

    if (simu_.tva.IsWritingTime()) {
      WriteQfile(dfib_, simu_, pressure_, vel_, *grid_);
    }
    getCentProfile(simu_, vel_, global_domain_, *grid_);
  }

  void CreatEta() {
    if (grid_->cal_eta_imp) {
      grid_->cal_eta_imp->solver(dfib_, global_domain_);
    } 
  }
};

}  // namespace projection_method
