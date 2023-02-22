#pragma once

#include "../grid/structured_grid.hpp"
#include "../matrix/csr_sparse_mat.hpp"
#include "../matrix/ell_sparse_mat.hpp"
#include "../pressure/pressure_mat.hpp"
#include "pressure_mat_csr.hpp"
#include <vector>

using MatType = typename mat::CsrMat<double>;

class PressureMat{
  public:
  std::vector<double> mat_b;
  std::vector<double> x_result;
  MatType mat_a;

  PressureMat(StructuredGrid & grid) : grid_(&grid) {
    mat_b.resize(grid_->cal_no_grid);
    x_result.resize(grid_->cal_no_grid, 0.0);
    mat_a.Set(GetCsrMat(grid));
  }

  void ConvertResultToPressure(Pressure& pressure, const CalDomain& domain) {
#pragma omp parallel for
    for (size_t i = domain.x_start; i < domain.x_end; ++i)
      for (size_t j = domain.y_start; j < domain.y_end; ++j)
        for (size_t k = domain.z_start; k < domain.z_end; ++k) {
          pressure(i, j, k) = x_result[grid_->icelCal(i, j, k)];
        }
  }

  void CalMatB(StaggeredVelocity& curr_vel, double dt,const CalDomain& domain) {
    double ddt = 1.0 / dt;

#pragma omp parallel for
    for (auto i = domain.x_start; i < domain.x_end; ++i) {
      for (auto j = domain.y_start; j < domain.y_end; ++j) {
        for (auto k = domain.z_start; k < domain.z_end; ++k) {
          auto row = grid_->icelCal(i, j, k);
          mat_b[row] =
              -ddt *
              ((curr_vel.U(i, j, k) - curr_vel.U(i - 1, j, k)) / grid_->dx[i] +
               (curr_vel.V(i, j, k) - curr_vel.V(i, j - 1, k)) / grid_->dy[j] +
               (curr_vel.W(i, j, k) - curr_vel.W(i, j, k - 1)) / grid_->dz[k]);
        }
      }
    }
  }

 private:
  StructuredGrid* grid_;
};

