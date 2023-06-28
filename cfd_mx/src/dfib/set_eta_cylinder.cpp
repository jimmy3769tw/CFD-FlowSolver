#include "dfib/set_eta_cylinder.hpp"

void CylinderZ::solver(ImmersedBoundary& dfib, LocalDomain& domain) {
    std::cout << "radius: " << radius << std::endl;
    std::cout << "[no_sub_nx, no_sub_ny] " << no_sub_nx << ", " << no_sub_ny << std::endl;
    std::cout << "[no_sub_nx, no_sub_ny] " << no_sub_nx << ", " << no_sub_ny << std::endl;
    

    std::vector<double> sdx(no_sub_nx + 1);
    std::vector<double> sdy(no_sub_ny + 1);
    for (auto i = domain.x_start; i < domain.x_end; ++i) {

        const double dxg = grid_->dx[i] / double(no_sub_nx);

        for (auto j = domain.y_start; j < domain.y_end; ++j) {
          const double grid_2d_diagonal =
              sqrt(std::pow(grid_->dx[i], 2) + std::pow(grid_->dy[j], 2)) * 0.5;
          const double distance = sqrt(std::pow(grid_->Xc(i) - center_x, 2) +
                                       std::pow(grid_->Yc(j) - center_y, 2));
          if (std::abs(distance - radius) < grid_2d_diagonal) {
            size_t eta_cnt = 0;
            const double dyg = grid_->dy[j] / double(no_sub_ny);
            for (size_t ii = 0; ii < no_sub_nx + 1; ii++)
              for (size_t jj = 0; jj < no_sub_nx + 1; jj++) {
                sdy[jj] = grid_->x_pos[j] + jj * dyg;
                sdx[ii] = grid_->x_pos[i] + ii * dxg;
              }

            for (size_t ii = 0; ii < no_sub_nx; ii++)
              for (size_t jj = 0; jj < no_sub_nx; jj++) {
                auto distance_sub =
                    sqrt(std::pow((sdx[ii] + sdx[ii + 1]) * 0.5 - center_x, 2) +
                         std::pow((sdy[jj] + sdy[jj + 1]) * 0.5 - center_y, 2));

                if (distance_sub <= radius) {
                  ++eta_cnt;
                }
              }
            dfib.Eta(i, j, domain.z_start) =
                double(eta_cnt) / no_sub_nx / no_sub_ny;
          }
          else if (distance > radius) {
            dfib.Eta(i, j, domain.z_start) = 0.;
          }
          else {
            dfib.Eta(i, j, domain.z_start) = 1.;
          }
        }
    }

    for (auto i = domain.x_start; i < domain.x_end; ++i) {
        for (auto j = domain.y_start; j < domain.y_end; ++j) {
          for (auto k = domain.z_start + 1; k < domain.z_end; ++k) {
            dfib.Eta(i, j, k) = dfib.Eta(i, j, domain.z_start);
          }
        }
    }
}

