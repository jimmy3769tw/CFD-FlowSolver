#pragma once

#include "../backend/physical_variables.hpp"
#include "../backend/domain.hpp"
#include "../grid/structured_grid.hpp"

void UpdateForceAndVelocity(ImmersedBoundary& dfib, const double dt,
                            const Pressure& pressure, const StaggeredVelocity& old_vel,
                            StaggeredVelocity& curr_vel, const CalDomain& location, StructuredGrid& grid) {

    auto nx = grid.nx;
    auto ny = grid.ny;
    auto nz = grid.nz;
    auto gc = grid.no_ghost_cell;

    double uS = 0.0, vS = 0.0, wS = 0.0;

    int one = 1;
    int d = 0;

    // * ---------------- In x direction ----------------
    if (location.x_end == nx - gc) {
      d = one;
    } else {
      d = 0;
    }

#pragma omp parallel for firstprivate(dt, location)
    for (auto i = location.x_start; i < location.x_end - d; ++i)
      for (auto j = location.y_start; j < location.y_end; ++j)
        for (auto k = location.z_start; k < location.z_end; ++k) {
          auto temp_u = old_vel.U(i, j, k)- dt * (pressure(i+1, j, k) - pressure(i, j, k)) / grid.staggered_dx[i];
          auto eta = dfib.Eta(i, j, k) + ( dfib.Eta(i+1, j, k) - dfib.Eta(i, j, k)) *
                                        grid.GetInterpolationX(i);
          curr_vel.U(i, j, k) = dfib.Eta(i, j, k) * uS + (1.0 - eta) * temp_u;
          dfib.Fx(i, j, k) = (curr_vel.U(i, j, k) - temp_u) / dt;
        }

    if (location.y_end == ny - gc) {
      d = one;
    } else {
      d = 0;
    }

#pragma omp parallel for firstprivate(dt, location)
    for (auto i = location.x_start; i < location.x_end; ++i)
      for (auto j = location.y_start; j < location.y_end - d; ++j)
        for (auto k = location.z_start; k < location.z_end; ++k) {
          auto temp_v = old_vel.V(i, j, k) - dt * (pressure(i, j+1, k) - pressure(i, j, k)) / grid.staggered_dy[j];
          auto eta = dfib.Eta(i, j, k) +
                     (dfib.Eta(i, j+1, k) - dfib.Eta(i, j, k)) * grid.GetInterpolationY(j);
          curr_vel.V(i, j, k) = dfib.Eta(i, j, k) * vS + (1.0 - eta) * temp_v;
          dfib.Fy(i, j, k)= (curr_vel.V(i, j, k) - temp_v) / dt;
        }

    // * ---------------- In z direction ----------------
    if (location.z_end == nz - gc) {
      d = one;
    } else {
      d = 0;
    }

#pragma omp parallel for firstprivate(dt, location)
    for (auto i = location.x_start; i < location.x_end; ++i)
      for (auto j = location.y_start; j < location.y_end; ++j)
        for (auto k = location.z_start; k < location.z_end - d; ++k) {
          auto temp_w = old_vel.W(i, j, k)- dt * (pressure(i, j, k+1) - pressure(i, j, k)) / grid.staggered_dz[k];
          auto eta = dfib.Eta(i, j, k) +
                     (dfib.Eta(i, j, k+1)- dfib.Eta(i, j, k)) * grid.GetInterpolationZ(k);
          curr_vel.W(i, j, k) = dfib.Eta(i, j, k) * wS + (1.0 - eta) * temp_w;
          dfib.Fz(i, j, k) = (curr_vel.W(i, j, k) - temp_w) / dt;
        }
}

void RemoveVelocity(const ImmersedBoundary& dfib, StaggeredVelocity& vel,
                    const CalDomain& location) {
#pragma omp parallel firstprivate(location)
    {
#pragma omp for simd nowait
      for (auto i = location.x_start; i < location.x_end; ++i)
        for (auto j = location.y_start; j < location.y_end; ++j)
          for (auto k = location.z_start; k < location.z_end; ++k) {
            vel.U(i, j, k) *= 0.5 * (dfib.Eta(i, j, k) + dfib.Eta(i + 1, j, k));
            vel.V(i, j, k) *= 0.5 * (dfib.Eta(i, j, k) + dfib.Eta(i, j + 1, k));
            vel.W(i, j, k) *= 0.5 * (dfib.Eta(i, j, k) + dfib.Eta(i, j, k + 1));
          }
    }
}
