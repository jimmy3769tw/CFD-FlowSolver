#pragma once
#include <utility>
#include <vector>

#include "../backend/physical_variables.hpp"
#include "../grid/structured_grid.hpp"


void UpdateAllVelocityOnBoundary(const CalDomain& location,
                                 StaggeredVelocity& vel, 
                                 Pressure& pre,
                                 StructuredGrid& grid) {
//   auto [num, dir] = getNumDIr();
  auto num = grid.bc_selector.num;
  auto dir = grid.bc_selector.dir;



#pragma omp parallel firstprivate(location) shared(grid)
  {
    auto nx = grid.nx;
    auto ny = grid.ny;
    auto nz = grid.nz;
    /*

    * |----|-----|----|-----|----|--    --|----|---|---|---|------|-----|------|
    * 0----1-----2----3-----|----|-- --|----|---|---|---n-4---n-3---n-2----n-1
    ! g1---g0---c0----c1----|----|--    --|----|---|---|---c1----c0----g0-----g1

    */
    static const std::vector<int> b0_g = {1, 0}, b0_c = {2, 3, 4},
                                  b2_g = {1, 0}, b2_c = {2, 3, 4},
                                  b4_g = {1, 0}, b4_c = {2, 3, 4},
                                  b1_g = {nx - 2, nx - 1},
                                  b1_c = {nx - 3, nx - 4, nx - 5},
                                  b3_g = {ny - 2, ny - 1},
                                  b3_c = {ny - 3, ny - 4, ny - 5},
                                  b5_g = {nz - 2, nz - 1},
                                  b5_c = {nz - 3, nz - 4, nz - 5};

    int d = 1;

#pragma omp for nowait
    for (auto j = location.y_start - d; j < location.y_end + d; ++j)
      for (auto k = location.z_start - d; k < location.z_end + d;
           ++k) {  // ====================== BC_0 ======================
        const int b = 0;
        vel.u[grid.icel(b0_g[1], j, k)] = vel.u[grid.icel(b0_g[0], j, k)] =
            dir[b][0] + num[b][0] * vel.u[grid.icel(b0_c[0], j, k)];
        vel.v[grid.icel(b0_g[1], j, k)] = vel.v[grid.icel(b0_g[0], j, k)] =
            dir[b][1] * 2. +
            (2.0 * num[b][1] - 1.0) * vel.v[grid.icel(b0_c[0], j, k)];
        vel.w[grid.icel(b0_g[1], j, k)] = vel.w[grid.icel(b0_g[0], j, k)] =
            dir[b][2] * 2. +
            (2.0 * num[b][2] - 1.0) * vel.w[grid.icel(b0_c[0], j, k)];
      }

#pragma omp for nowait
    for (auto j = location.y_start - d; j < location.y_end + d; ++j)
      for (auto k = location.z_start - d; k < location.z_end + d;
           ++k) {  // ====================== BC_1 ======================
        const int b = 1;
        vel.u[grid.icel(b1_g[1], j, k)] = vel.u[grid.icel(b1_g[0], j, k)] =
            vel.u[grid.icel(b1_c[0], j, k)] =
                dir[b][0] + num[b][0] * vel.u[grid.icel(b1_c[1], j, k)];
        vel.v[grid.icel(b1_g[1], j, k)] = vel.v[grid.icel(b1_g[0], j, k)] =
            dir[b][1] * 2. +
            (2.0 * num[b][1] - 1.0) * vel.v[grid.icel(b1_c[0], j, k)];
        vel.w[grid.icel(b1_g[1], j, k)] = vel.w[grid.icel(b1_g[0], j, k)] =
            dir[b][2] * 2. +
            (2.0 * num[b][2] - 1.0) * vel.w[grid.icel(b1_c[0], j, k)];
      }

#pragma omp for nowait
    for (auto i = location.x_start - d; i < location.x_end + d; ++i)
      for (auto k = location.z_start - d; k < location.z_end + d;
           ++k) {  // ====================== BC_2 ======================
        const int b = 2;
        vel.u[grid.icel(i, b2_g[1], k)] = vel.u[grid.icel(i, b2_g[0], k)] =
            dir[b][0] * 2. +
            (2.0 * num[b][0] - 1.0) * vel.u[grid.icel(i, b2_c[0], k)];
        vel.v[grid.icel(i, b2_g[1], k)] = vel.v[grid.icel(i, b2_g[0], k)] =
            dir[b][1] + num[b][1] * vel.v[grid.icel(i, b2_c[0], k)];
        vel.w[grid.icel(i, b2_g[1], k)] = vel.w[grid.icel(i, b2_g[0], k)] =
            dir[b][2] * 2. +
            (2.0 * num[b][2] - 1.0) * vel.w[grid.icel(i, b2_c[0], k)];
      }

#pragma omp for nowait
    for (auto i = location.x_start - d; i < location.x_end + d; ++i)
      for (auto k = location.z_start - d; k < location.z_end + d;
           ++k) {  // ====================== BC_3 ======================
        const int b = 3;
        vel.v[grid.icel(i, b3_g[1], k)] = vel.v[grid.icel(i, b3_g[0], k)] =
            vel.v[grid.icel(i, b3_c[0], k)] =
                dir[b][1] + num[b][1] * vel.v[grid.icel(i, b3_c[1], k)];

        vel.u[grid.icel(i, b3_g[1], k)] = vel.u[grid.icel(i, b3_g[0], k)] =
            dir[b][0] * 2. +
            (2.0 * num[b][0] - 1.0) * vel.u[grid.icel(i, b3_c[0], k)];
        vel.w[grid.icel(i, b3_g[1], k)] = vel.w[grid.icel(i, b3_g[0], k)] =
            dir[b][2] * 2. +
            (2.0 * num[b][2] - 1.0) * vel.w[grid.icel(i, b3_c[0], k)];
      }

#pragma omp for nowait
    for (auto i = location.x_start - d; i < location.x_end + d; ++i)
      for (auto j = location.y_start - d; j < location.y_end + d;
           ++j) {  // ====================== BC_4 ======================
        const int b = 4;
        vel.u[grid.icel(i, j, b4_g[1])] = vel.u[grid.icel(i, j, b4_g[0])] =
            dir[b][0] * 2. +
            (2.0 * num[b][0] - 1.0) * vel.u[grid.icel(i, j, b4_c[0])];
        vel.v[grid.icel(i, j, b4_g[1])] = vel.v[grid.icel(i, j, b4_g[0])] =
            dir[b][1] * 2. +
            (2.0 * num[b][1] - 1.0) * vel.v[grid.icel(i, j, b4_c[0])];
        vel.w[grid.icel(i, j, b4_g[1])] = vel.w[grid.icel(i, j, b4_g[0])] =
            dir[b][2] + num[b][2] * vel.w[grid.icel(i, j, b4_c[0])];
      }

#pragma omp for nowait
    for (auto i = location.x_start - d; i < location.x_end + d; ++i)
      for (auto j = location.y_start - d; j < location.y_end + d;
           ++j) {  // ====================== BC_5 ======================
        const int b = 5;
        vel.w[grid.icel(i, j, b5_g[1])] = vel.w[grid.icel(i, j, b5_g[0])] =
            vel.w[grid.icel(i, j, b5_c[0])] =
                dir[b][2] + num[b][2] * vel.w[grid.icel(i, j, b5_c[1])];

        vel.u[grid.icel(i, j, b5_g[1])] = vel.u[grid.icel(i, j, b5_g[0])] =
            dir[b][0] * 2. +
            (2.0 * num[b][0] - 1.0) * vel.u[grid.icel(i, j, b5_c[0])];
        vel.v[grid.icel(i, j, b5_g[1])] = vel.v[grid.icel(i, j, b5_g[0])] =
            dir[b][1] * 2. +
            (2.0 * num[b][1] - 1.0) * vel.v[grid.icel(i, j, b5_c[0])];
      }

// * =================== pressure =================
#pragma omp for nowait
    for (auto j = location.y_start - d; j < location.y_end + d; ++j)
      for (auto k = location.z_start - d; k < location.z_end + d; ++k) {
        pre.p[grid.icel(b0_g[1], j, k)] = pre.p[grid.icel(b0_g[0], j, k)] =
            pre.p[grid.icel(b0_c[0], j, k)];
        pre.p[grid.icel(b1_g[1], j, k)] = pre.p[grid.icel(b1_g[0], j, k)] =
            pre.p[grid.icel(b1_c[0], j, k)];
      }

#pragma omp for nowait
    for (auto i = location.x_start - d; i < location.x_end + d; ++i)
      for (auto k = location.z_start - d; k < location.z_end + d; ++k) {
        pre.p[grid.icel(i, b2_g[1], k)] = pre.p[grid.icel(i, b2_g[0], k)] =
            pre.p[grid.icel(i, b2_c[0], k)];
        pre.p[grid.icel(i, b3_g[1], k)] = pre.p[grid.icel(i, b3_g[0], k)] =
            pre.p[grid.icel(i, b3_c[0], k)];
      }

#pragma omp for
    for (auto i = location.x_start - d; i < location.x_end + d; ++i)
      for (auto j = location.y_start - d; j < location.y_end + d; ++j) {
        pre.p[grid.icel(i, j, b4_g[1])] = pre.p[grid.icel(i, j, b4_g[0])] =
            pre.p[grid.icel(i, j, b4_c[0])];
        pre.p[grid.icel(i, j, b5_g[1])] = pre.p[grid.icel(i, j, b5_g[0])] =
            pre.p[grid.icel(i, j, b5_c[0])];
      }

  }  // Omp fork join
}