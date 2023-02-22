#pragma once
#include <cmath>

#include "controlPanel.hpp"

void DfibInit(ImmersedBoundary& dfib, CalDomain& domain, grid& gridA) {
  auto [nx, ny, nz, gC] = gridA.nxyzgC;
  for (auto i = domain.x_start; i < domain.x_end; ++i)
    for (auto j = domain.y_start; j < domain.y_end; ++j)
      for (auto k = domain.z_start; k < domain.z_end; ++k) {
        dfib.eta[gridA.icel(i, j, k)] = 0.0;

        if (i == 2) {
          dfib.eta[gridA.icel(i - 1, j, k)] = 0.0;
          dfib.eta[gridA.icel(i - 2, j, k)] = 0.0;
        }

        if (i == nx - 3) {
          dfib.eta[gridA.icel(i + 1, j, k)] = 0.0;
          dfib.eta[gridA.icel(i + 2, j, k)] = 0.0;
        }

        if (j == 2) {
          dfib.eta[gridA.icel(i, j - 1, k)] = 0.0;
          dfib.eta[gridA.icel(i, j - 2, k)] = 0.0;
        }

        if (j == nx - 3) {
          dfib.eta[gridA.icel(i, j + 1, k)] = 0.0;
          dfib.eta[gridA.icel(i, j + 2, k)] = 0.0;
        }

        if (k == 2) {
          dfib.eta[gridA.icel(i, j, k - 1)] = 0.0;
          dfib.eta[gridA.icel(i, j, k - 2)] = 0.0;
        }

        if (k == nz - 3) {
          dfib.eta[gridA.icel(i, j, k + 1)] = 0.0;
          dfib.eta[gridA.icel(i, j, k + 2)] = 0.0;
        }
      }
}

// * for X direction
void CalEtaCylinderXdirection(ImmersedBoundary& dfib, CalDomain& domain,
                              StructGrid& grid) {
  auto [nx, ny, nz] = grid.nxyz;

  double center_z = grid.lz / 2.0, center_y = grid.ly / 2.0;

  double radius = grid.lz / 2.0;  // * setting

  int no_sub = 100;

  vector<double> sdy(no_sub), sdz(no_sub);

  int one = 1;
  int zero = 0;

  for (size_t j = 0; j < ny - 4; ++j)
    for (size_t k = 0; k < nz - 4; ++k) {
      const int icel = grid.icel(zero, j, k);
      const double Diagonal = sqrt(  // gridA.Dx[0]*gridA.Dx[0] +
                                  grid.Dy[j] * grid.Dy[j] + grid.Dz[k] * grid.Dz[k]) *
                              0.5;

      const double distance =
          sqrt(pow(((grid.Z[k] + grid.Z[k + 1]) * 0.5) - center_z, 2) +
               pow(((grid.Y[j] + grid.Y[j + 1]) * 0.5) - center_y, 2));

      if (abs(distance - radius) < Diagonal) {
        size_t xi = 0;
        const double dyg = grid.Dy[j] / double(no_sub);
        const double dzg = grid.Dz[k] / double(no_sub);

        for (size_t jj; jj < no_sub; jj++)
          for (size_t kk; kk < no_sub; kk++) {
            sdy[jj] = grid.Z[j] + (jj)*dyg;
            sdz[jj] = grid.Z[k] + (kk)*dzg;
          }

        for (size_t jj; jj < no_sub - 1; jj++)
          for (size_t kk; kk < no_sub - 1; kk++) {
            auto Sdist = sqrt(sdy[jj] + sdy[jj + 1]);
            if (Sdist <= radius) {
              ++xi;
            }
          }
        dfib.eta[icel] = double(xi) * pow(no_sub, -2);

      } else if (distance > radius) {
        dfib.eta[icel] = 0.;
      } else {
        dfib.eta[icel] = 1.;
      }
    }

  for (size_t i = 1; i < nx - 1; ++i)
    for (size_t j = 0; j < ny - 1; ++j)
      for (size_t k = 0; k < nz - 1; ++k) {
        const int icel = grid.icel(i, j, k);
        dfib.eta[icel] = grid.icel(one, j, k);
      }
}

void DFIB_CylinderZ(ImmersedBoundary& dfib, CalDomain& domain, grid& grid) {
  auto [nx, ny, nz, gC] = grid.nxyzgC;

  // * ------------------------- setting
  double center_x = 6.5, center_y = grid.ly / 2.0, radius = 0.5;
  dfib.cylinderDimension = radius * 2;
  dfib.cylinderCenter.resize(2);
  dfib.cylinderCenter[0] = center_x;
  dfib.cylinderCenter[1] = center_y;
  int nSubGrids = 100;
  // * ------------------------- setting

  vector<double> sdy(nSubGrids + 1);
  vector<double> sdx(nSubGrids + 1);
  double Sdist;

  for (size_t i = gC; i < nx - gC; ++i)
    for (size_t j = gC; j < ny - gC; ++j) {
      const int icel = i * nz * ny + j * nz;
      const double Diagonal = sqrt(grid.Dx[i] * grid.Dx[i] + grid.Dy[j] * grid.Dy[j]) *
                              0.5;  // FIXME :: this is three way??
      const double distance =
          sqrt(pow(((grid.X[i] + grid.X[i + 1]) * 0.5) - center_x, 2) +
               pow(((grid.Y[j] + grid.Y[j + 1]) * 0.5) - center_y, 2));

      if (std::abs(distance - radius) < Diagonal) {
        size_t xi = 0;
        const double dyg = grid.Dy[j] / double(nSubGrids);
        const double dxg = grid.Dx[i] / double(nSubGrids);
        // * | | | |
        for (size_t ii = 0; ii < nSubGrids + 1; ii++)
          for (size_t jj = 0; jj < nSubGrids + 1; jj++) {
            sdy[jj] = grid.Y[j] + (jj)*dyg;
            sdx[ii] = grid.X[i] + (ii)*dxg;
          }

        for (size_t ii = 0; ii < nSubGrids; ii++)
          for (size_t jj = 0; jj < nSubGrids; jj++) {
            auto Sdist = sqrt(pow((sdx[ii] + sdx[ii + 1]) * 0.5 - center_x, 2) +
                              pow((sdy[jj] + sdy[jj + 1]) * 0.5 - center_y, 2));

            if (Sdist <= radius) {
              ++xi;
            }
          }
        dfib.eta[icel] = double(xi) * pow(nSubGrids, -2);

      } else if (distance > radius) {
        dfib.eta[icel] = 0.;
      } else {
        dfib.eta[icel] = 1.;
      }
    }

  for (size_t i = 0; i < nx - 1; ++i)
    for (size_t j = 0; j < ny - 1; ++j)
      for (size_t k = 1; k < nz - 1; ++k)
        dfib.eta[grid.icel(i, j, k)] = dfib.eta[grid.icel(i, j, 0)];
}

void PotentialFlow(Simulation& simu, ImmersedBoundary& dfib, velocity& vel,
                   pressure& pre, CalDomain& domain, grid& grid) {
  auto [nx, ny, nz, gC] = grid.nxyzgC;

  double U_infty = 1;
  double pi = acos(-1);

  auto D_c = dfib.cylinderDimension;
  auto x_D = dfib.cylinderCenter[0];
  auto y_D = dfib.cylinderCenter[1];
  for (auto i = grid.gC; i < grid.nx - grid.gC; ++i)
    for (auto j = grid.gC; j < grid.ny - grid.gC; ++j) {
      auto tX = (grid.X[i] + grid.X[i + 1]) * 0.5 - x_D;
      auto tY = (grid.Y[j] + grid.Y[j + 1]) * 0.5 - y_D;
      auto tempX = pow((tX), 2);
      auto tempY = pow((tY), 2);
      // vel.u[grid.icel(i,j,gC)] += -D_c * (tempX + tempY) / (2 * pi) * (tempX -
      // tempY) ; vel.v[grid.icel(i,j,gC)] += -D_c * (tempX + tempY) / (2.* pi) *
      // (2*(tX*tY)); pre.p[grid.icel(i,j,1)] = -0.5* simu.nu *(
      //     pow(vel.u[grid.icel(i,j,1)], 2.0) + pow(vel.v[grid.icel(i,j,1)], 2.0)
      //     );
    }

  for (auto i = gC + 0; i < nx - 2; ++i)
    for (auto j = gC + 0; j < ny - 2; ++j)
      for (auto k = gC + 1; k < nz - 2; ++k) {
        auto ick1 = grid.icel(i, j, gC);
        auto ic = grid.icel(i, j, k);
        vel.u[ic] = vel.u[ick1];
        pre.p[ic] = pre.p[ick1];
        vel.v[ic] = vel.v[ick1];
      }
}