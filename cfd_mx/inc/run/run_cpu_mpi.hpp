#pragma once
#include <algorithm>

#include "../backend/domain.hpp"
#include "../backend/physical_variables.hpp"
#include "../backend/print_openmp.hpp"
#include "../backend/simulation.hpp"
#include "../boundary_condition/boundary_condition.hpp"
#include "../boundary_condition/boundary_condition_copy.hpp"
#include "../dfib/update_u_and_force.hpp"
#include "../grid/structured_grid.hpp"
#include "../io/csv/csv_structured_grid.hpp"
#include "../io/io_tools.hpp"
#include "../io/plot3d/write_qfile.hpp"
#include "../io/plot3d/write_xfile.hpp"
#include "../matrix/solver/mpi/bicgstab_restart_mpi.hpp"
#include "../mpi_tool/mpi_tool.hpp"
#include "../pressure/pressure_mat.hpp"
#include "../source/convection_and_difussion.hpp"
#include "../matrix/solver/mpi/bicgstab_mpi.hpp"
#include "mpi.h"
#include "omp.h"

// MpiTool(bool &reorder, int argc, char **argv, StructuredGrid &grid,
//         LocalDomain &domain)
//     : grid_(&grid), domain_(domain) 

    
namespace projection_method {
  class CpuOpenMpMpi {
    private:
    using SolverType = solver::BicgstabMpi<MatType>;

   public:
    ~CpuOpenMpMpi() { MPI_Finalize(); }
    CpuOpenMpMpi(Simulation &simu, StructuredGrid &grid)
        : simu_(simu),
          grid_(&grid),
          vel_(StaggeredVelocity(grid)),
          intermediate_vel_(StaggeredVelocity(grid)),
          mpi_tool_(mpi::MpiTool(simu.GetArgc(), simu.GetArgv(), grid)),
          dfib_(ImmersedBoundary(grid)),
          pressure_mat_(PressureMat(grid)),
          local_domain_(LocalDomain(grid)),
          global_domain_(LocalDomain(grid)),
          pressure_(Pressure(grid))

    {
      vel_.FillVel(simu.ini_condition.u, simu.ini_condition.v, simu.ini_condition.w);

      if (mpi_tool_.IsMaster()) {
        CreatOutputFile();
        PrintIsOpenmpExist();
      }
      std::vector<int> grid_size{grid.nx, grid.ny, grid.nz};
      global_domain_.Init({1, 1, 1}, {0, 0, 0});
      local_domain_.Init({mpi::GetSize(), 1, 1}, {mpi::GetRank(), 0, 0});
      mpi_tool_.SetDomain(local_domain_);
      simu.pid = mpi_tool_.GetRank();
      if (mpi_tool_.IsMaster()) {
        csv::WriteCsvFile("Information/grid.csv", *grid_);
        plot3d::write_xfile(grid);
      }
      UpdateAllVelocityOnBoundary(global_domain_, vel_, pressure_, grid);
      CopyVelocityOnBoundary(global_domain_, vel_, intermediate_vel_, grid);
      auto [start, end] = local_domain_.GetOneDimStartEnd();
      pressure_mat_.mat_a.InitMpi(mpi_tool_.GetCommWorld(), start, end);
    }

    void solve() {
      pressure_solver_ = SolverType(pressure_mat_.mat_a);
      pressure_solver_.SetTolerance(simu_.poisson_criteria);
      CreatEta();
      if (mpi_tool_.IsMaster()) {
      WriteQfile(dfib_, simu_, pressure_, vel_, *grid_);
      }
      for (; !simu_.tva.IsLoopFinish(); simu_.tva.AddLoop()) {
        CalProjectionMethod();
      }
      mpi_tool_.Barrier();
    }

   private:
    mpi::MpiTool mpi_tool_;
    bool reorder_ = true;
    Simulation simu_;
    StructuredGrid *grid_;
    StaggeredVelocity vel_;
    StaggeredVelocity intermediate_vel_;
    Pressure pressure_;
    LocalDomain global_domain_;
    LocalDomain local_domain_;
    ImmersedBoundary dfib_;
    PressureMat pressure_mat_;
    SolverType pressure_solver_;

    void CalProjectionMethod() {
    
      #if defined(TERBULENCE_SMAGORINSKY)
            CalSmagorinskyModel(ShareM, simu_, vel_, t1, local_domain_, *grid_);
      #endif

      CalConvectionAndDiffusion(simu_, vel_, intermediate_vel_, local_domain_, *grid_);
      mpi_tool_.SendRecv(grid_->no_ghost_cell, vel_.u)
          .SendRecv(grid_->no_ghost_cell, vel_.v)
          .SendRecv(grid_->no_ghost_cell, vel_.w)
          .Barrier();

      pressure_mat_.CalMatB(intermediate_vel_, simu_.tva.GetDt(), local_domain_);
      auto [iters, error] =
          pressure_solver_(pressure_mat_.mat_b, pressure_mat_.x_result);
     
      if(mpi_tool_.IsMaster()) {
        std::cout << simu_.tva.GetLoop() << "\t [iters, error]\t" << iters << '\t' << error << std::endl;
      }

      pressure_mat_.ConvertResultToPressure(pressure_, local_domain_);
      mpi_tool_.Barrier().SendRecvRight(1, pressure_.p).Barrier();

      UpdateForceAndVelocity(dfib_, simu_.tva.GetDt(), pressure_, intermediate_vel_,
                             vel_, local_domain_, *grid_);

      UpdateAllVelocityOnBoundary(local_domain_, vel_, pressure_, *grid_);

      // ! Do not 
      // BC_staggered_copy(global_domain_, vel_, intermediate_vel_, grid_);
      // ! Do not 

      if (simu_.tva.IsWritingTime()) {
      mpi_tool_.Barrier()
          .CollectToMaster(vel_.u)
          .CollectToMaster(vel_.v)
          .CollectToMaster(vel_.w)
          .CollectToMaster(pressure_.p)
          .Barrier();

      if (mpi_tool_.IsMaster()) {
        WriteQfile(dfib_, simu_, pressure_, vel_, *grid_);
        }
      }
    }

    void CreatEta() {

    }
  };

}  // namespace projection_method
