#pragma once
#include "../grid/structured_grid.hpp"

auto GetCsrMat(StructuredGrid& grid) {
  auto cal_nx = grid.cal_nx;
  auto cal_ny = grid.cal_ny;
  auto cal_nz = grid.cal_nz;
  auto no_ghost_cell = grid.no_ghost_cell;

  auto xp = [&](auto i) {
    return 1.0 / grid.dx.at(i + no_ghost_cell) /
           grid.staggered_dx.at(i + no_ghost_cell);
  };

  auto yp = [&](auto j) {
    return 1.0 / grid.dy.at(j + no_ghost_cell) /
           grid.staggered_dy.at(j + no_ghost_cell);
  };
  auto zp = [&](auto k) {
    return 1.0 / grid.dz.at(k + no_ghost_cell) /
           grid.staggered_dz.at(k + no_ghost_cell);
  };

  auto xm = [&](auto i) {
    return 1.0 / grid.dx.at(i + no_ghost_cell) /
           grid.staggered_dx.at(i + no_ghost_cell - 1);
  };

  auto ym = [&](auto j) {
    return 1.0 / grid.dy.at(j + no_ghost_cell) /
           grid.staggered_dy.at(j + no_ghost_cell - 1);
  };

  auto zm = [&](auto k) {
    return 1.0 / grid.dz.at(k + no_ghost_cell) /
           grid.staggered_dz.at(k + no_ghost_cell - 1);
  };

  std::vector<int> ptr;
  std::vector<int> idx;
  std::vector<double> val;
  ptr.reserve(grid.cal_no_grid + 1);
  ptr.push_back(0);
  idx.reserve(grid.cal_no_grid * 7);
  val.reserve(grid.cal_no_grid * 7);

  int x_start = 0, y_start = 0, z_start = 0;
  auto x_end = cal_nx - 1;
  auto y_end = cal_ny - 1;
  auto z_end = cal_nz - 1;

  for (int i = 0, i_count = 0; i < cal_nx; ++i) {
    for (int j = 0; j < cal_ny; ++j) {
      for (int k = 0; k < cal_nz; ++k, ++i_count) {
        // Boundary point. Use Num condition.
        double aa = 0.0, bb = 0.0, cc = 0.0;
        if (i == x_start) {
          aa = -xm(i);
        }

        if (i == x_end) {
          aa = -xp(i);
        }

        if (j == y_start) {
          bb = -ym(j);
        }

        if (j == y_end) {
          bb = -yp(j);
        }

        if (k == z_start) {
          cc = -zm(k);
        }

        if (k == z_end) {
          cc = -zp(k);
        }
        auto ijk = xm(i) + ym(j) + zm(k) + xp(i) + yp(j) + zp(k) + aa + bb + cc;

        if (i != x_start) {
          idx.push_back(i_count - grid.cal_ny_times_nz);
          val.push_back(-xm(i));
        }

        if (j != y_start) {
          idx.push_back(i_count - grid.cal_nz);
          val.push_back(-ym(j));
        }

        if (k != z_start) {
          idx.push_back(i_count - 1);
          val.push_back(-zm(k));
        }

        idx.push_back(i_count);
        val.push_back(ijk);

        if (k != z_end) {
          idx.push_back(i_count + 1);
          val.push_back(-zp(k));
        }

        if (j != y_end) {
          idx.push_back(i_count + grid.cal_nz);
          val.push_back(-yp(j));
        }

        if (i != x_end) {
          idx.push_back(i_count + grid.cal_ny_times_nz);
          val.push_back(-xp(i));
        }

        ptr.push_back(idx.size());
      }
    }
  }

  return std::make_tuple(ptr, idx, val);
}
