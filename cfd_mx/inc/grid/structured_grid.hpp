#pragma once
#include <dirent.h>
#include <sys/types.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>
#include <memory>

#include "boundary_condition/boundary_condition_selector.hpp"

class CalEtaInterface;

// namespace::structureGrid
class StructuredGrid {
 private:
  std::vector<double> delta_;
  std::vector<std::string> direction_order_{"[X]", "[Y]", "[Z]"};
  bool isSetLen_ = false;
  std::string gridTypeName_;
  virtual void InitPos() = 0;

 protected:
  double lx_, ly_, lz_;

 public:
  explicit StructuredGrid(int cal_nx, int cal_ny, int cal_nz)
      : cal_nx(cal_nx), cal_ny(cal_ny), cal_nz(cal_nz) {
    resize();
  }

  BoundaryConditionSelector bc_selector;

  std::shared_ptr<CalEtaInterface> cal_eta_imp;

  StructuredGrid& setLen(int lx, int ly, int lz) {
    lx_ = lx;
    ly_ = ly;
    lz_ = lz;
    isSetLen_ = true;
    return *this;
  }
  

  std::tuple<double, double, double> GetLen(){
    return std::make_tuple(lx_, ly_, lz_);
  }


  void Init() {
    CheckParameter();
    InitPos();
    InitGhostCellPos();
    InitGridDistance();
    InitCentPos();
  }
  
  int cal_nx, cal_ny, cal_nz;
  const int no_ghost_cell = 2;
  const int cal_no_grid = cal_nx * cal_nz * cal_ny;
  const int cal_ny_times_nz = cal_ny * cal_nz;
  const int nx = cal_nx + 2 * no_ghost_cell;
  const int ny = cal_ny + 2 * no_ghost_cell;
  const int nz = cal_nz + 2 * no_ghost_cell;
  const int ny_times_nz = ny*nz;
  const int nx_times_ny = ny*nx;
  const int no_grid = nx * ny * nz;
  
  std::vector<double> dx, dy, dz;
  std::vector<double> staggered_dx, staggered_dy, staggered_dz;
  std::vector<double> x_pos, y_pos, z_pos;
  std::vector<double> x_cent_pos, y_cent_pos, z_cent_pos;

  void resize() {
    delta_.resize(cal_no_grid);

    x_cent_pos.resize(cal_nx);
    y_cent_pos.resize(cal_ny);
    z_cent_pos.resize(cal_nz);

    dx.resize(nx);
    dy.resize(ny);
    dz.resize(nz);

    staggered_dx.resize(nx);
    staggered_dy.resize(ny);
    staggered_dz.resize(nz);

    x_pos.resize(nx + 1);
    y_pos.resize(ny + 1);
    z_pos.resize(nz + 1);
  }

  void CheckParameter() {
    if (!isSetLen_) {
      throw std::invalid_argument(
          "Parameter settings have not been completed.\n");
    }
  }


  void InitCentPos() {
    for (int i = 0; i < nx - 2 * no_ghost_cell; ++i) {
      x_cent_pos.at(i) =
          (x_pos.at(i + no_ghost_cell) + x_pos.at(i + no_ghost_cell + 1)) / 2.0;
    }
    for (int j = 0; j < ny - 2 * no_ghost_cell; ++j) {
      y_cent_pos.at(j) =
          (y_pos.at(j + no_ghost_cell) + y_pos.at(j + no_ghost_cell + 1)) / 2.0;
    }
    for (int k = 0; k < nz - 2 * no_ghost_cell; ++k) {
      z_cent_pos.at(k) =
          (z_pos.at(k + no_ghost_cell) + z_pos.at(k + no_ghost_cell + 1)) / 2.0;
    }
  }

  double& Xc(int i) {
    return x_cent_pos.at(i-no_ghost_cell);
  }
  
  double& Yc(int i) {
    return y_cent_pos.at(i-no_ghost_cell);
  }

  double& Zc(int i) {
    return z_cent_pos.at(i-no_ghost_cell);
  }

  template <typename T>
  inline size_t icel(T i, T j, T k) {
    return i * ny_times_nz + j * nz + k;
  }

  template <typename T>
  std::tuple<int, int, int, int, int, int, int, int, int, int, int, int>
  Get12Nb(T& icel) {
    return std::make_tuple(icel - ny_times_nz, icel + ny_times_nz, icel - nz,
                      icel + nz, icel - 1, icel + 1, icel - 2 * ny_times_nz,
                      icel + 2 * ny_times_nz, icel - 2 * nz, icel + 2 * nz,
                      icel - 2, icel + 2);
  }

  int nzDir = (cal_nz + 2);
  int nznyDir = (cal_nz + 2) * (cal_ny + 2);

  template <typename T>
  inline int icelDir(T i, T j, T k) const {
    return (i - no_ghost_cell) * nznyDir + (j - no_ghost_cell) * nzDir +
           (k - no_ghost_cell);
  }


  template <typename T>
  inline size_t icelCal(T i, T j, T k) const {
    return (i - no_ghost_cell) * cal_ny_times_nz + (j - no_ghost_cell) * cal_nz + (k - no_ghost_cell);
  }

  std::string show() {
    std::stringstream ss;
    ss << ", [nx, ny, nz]=";
    ss << cal_nx << ", ";
    ss << cal_ny << ", ";
    ss << cal_nz << ", ";
    return ss.str();
  }

  inline auto GetInterpolationX(int i) { return (staggered_dx[i] * 0.5) / dx[i]; }
  inline auto GetInterpolationY(int j) { return (staggered_dy[j] * 0.5) / dy[j]; }
  inline auto GetInterpolationZ(int k) { return (staggered_dz[k] * 0.5) / dz[k]; }

  void InitGhostCellPos() {
    x_pos[1] = 2 * x_pos[2] - x_pos[3];
    y_pos[1] = 2 * y_pos[2] - y_pos[3];
    z_pos[1] = 2 * z_pos[2] - z_pos[3];

    x_pos[0] = 2 * x_pos[1] - x_pos[2];
    y_pos[0] = 2 * y_pos[1] - y_pos[2];
    z_pos[0] = 2 * z_pos[1] - z_pos[2];

    x_pos[nx - 1] = 2 * x_pos[nx - 2] - x_pos[nx - 3];
    y_pos[ny - 1] = 2 * y_pos[ny - 2] - y_pos[ny - 3];
    z_pos[nz - 1] = 2 * z_pos[nz - 2] - z_pos[nz - 3];

    x_pos[nx] = 2 * x_pos[nx - 1] - x_pos[nx - 2];
    y_pos[ny] = 2 * y_pos[ny - 1] - y_pos[ny - 2];
    z_pos[nz] = 2 * z_pos[nz - 1] - z_pos[nz - 2];
  }

  bool InitGridDistance() {

    for (auto i = no_ghost_cell; i < nx - no_ghost_cell - 1; ++i) {
      staggered_dx[i] = (x_pos[i + 2] - x_pos[i]) / 2.0;
      dx[i] = (x_pos[i + 1] - x_pos[i]);
    }

    for (auto j = no_ghost_cell; j < ny - no_ghost_cell - 1; ++j) {
      staggered_dy[j] = (y_pos[j + 2] - y_pos[j]) / 2.0;
      dy[j] = (y_pos[j + 1] - y_pos[j]);
    }

    for (auto k = no_ghost_cell; k < nz - no_ghost_cell - 1; ++k) {
      staggered_dz[k] = (z_pos[k + 2] - z_pos[k]) / 2.0;
      dz[k] = (z_pos[k + 1] - z_pos[k]);
    }

    // Ghost grid
    dx[0] = dx[1] = dx[2];
    dy[0] = dy[1] = dy[2];
    dz[0] = dz[1] = dz[2];

    staggered_dx[0] = staggered_dx[1] = staggered_dx[2];
    staggered_dy[0] = staggered_dy[1] = staggered_dy[2];
    staggered_dz[0] = staggered_dz[1] = staggered_dz[2];

    dx[nx - 1] = dx[nx - 2] = dx[nx - 3] = x_pos[nx - 2] - x_pos[nx - 3];
    dy[ny - 1] = dy[ny - 2] = dy[ny - 3] = y_pos[ny - 2] - y_pos[ny - 3];
    dz[nz - 1] = dz[nz - 2] = dz[nz - 3] = z_pos[nz - 2] - z_pos[nz - 3];

    staggered_dx[nx - 1] = staggered_dx[nx - 2] = staggered_dx[nx - 3] = staggered_dx[nx - 4];
    staggered_dy[ny - 1] = staggered_dy[ny - 2] = staggered_dy[ny - 3] = staggered_dy[ny - 4];
    staggered_dz[nz - 1] = staggered_dz[nz - 2] = staggered_dz[nz - 3] = staggered_dz[nz - 4];

    return true;
  }
};

// void checkgC() {
//   // check == x_pos[nx - no_ghost_cell] == lx
//   // -----------------------
//   if (abs(x_pos[nx - no_ghost_cell] - lx) > 1.0e-5) {
//     cout << nx - no_ghost_cell << "\t" << x_pos[nx - no_ghost_cell] << "\t"
//          << lx << std::endl;
//     throw std::invalid_argument("x_pos[nx - no_ghost_cell]");
//   }

//   if (abs(y_pos[ny - no_ghost_cell] - ly) > 1.0e-5) {
//     cout << ny - no_ghost_cell << "\t" << y_pos[ny - no_ghost_cell] << "\t"
//          << ly << std::endl;
//     throw std::invalid_argument("y_pos[ny - no_ghost_cell]");
//   }

//   if (abs(z_pos[nz - no_ghost_cell] - lz) > 1.0e-5) {
//     cout << nz - no_ghost_cell << "\t" << z_pos[nz - no_ghost_cell] << "\t"
//          << lz << std::endl;
//     throw std::invalid_argument("z_pos[nz - no_ghost_cell]");
//   }
//   // -----------------------

//   // -----------------------
//   for (auto i = 0; i < nx; ++i) {
//     if (x_pos[i] >= x_pos[i + 1]) throw std::invalid_argument("x_pos[i] i:");
//   }

//   for (auto j = 0; j < ny; ++j) {
//     if (y_pos[j] >= y_pos[j + 1]) throw std::invalid_argument("y_pos[j]");
//   }

//   for (auto k = 0; k < nz; ++k) {
//     if (z_pos[k] >= z_pos[k + 1]) throw std::invalid_argument("z_pos[k]");
//   }
// }

