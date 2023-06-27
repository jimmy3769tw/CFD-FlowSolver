#pragma once
#include "boundary_condition.hpp"

bool CopyVelocityOnBoundary(const CalDomain& location,
                            const StaggeredVelocity& old_vel,
                            StaggeredVelocity& curr_vel, StructuredGrid& grid) {
  int nx = grid.nx;
  int ny = grid.ny;
  int nz = grid.nz;

  const std::vector<int> b0_g = {1, 0}, b0_c = {2, 3, 4}, b2_g = {1, 0},
                         b2_c = {2, 3, 4}, b4_g = {1, 0}, b4_c = {2, 3, 4},
                         b1_g = {nx - 2, nx - 1},
                         b1_c = {nx - 3, nx - 4, nx - 5},
                         b3_g = {ny - 2, ny - 1},
                         b3_c = {ny - 3, ny - 4, ny - 5},
                         b5_g = {nz - 2, nz - 1},
                         b5_c = {nz - 3, nz - 4, nz - 5};
  int d = 1;

  for (int j = location.y_start - d; j < location.y_end + d; ++j)
    for (int k = location.z_start - d; k < location.z_end + d;
         ++k) {  // ====================== BC_0 and BC_1 ======================
      curr_vel.u[grid.icel(b0_g[0], j, k)] = old_vel.u[grid.icel(b0_g[0], j, k)];
      curr_vel.u[grid.icel(b0_g[1], j, k)] = old_vel.u[grid.icel(b0_g[1], j, k)];

      curr_vel.v[grid.icel(b0_g[0], j, k)] = old_vel.v[grid.icel(b0_g[0], j, k)];
      curr_vel.v[grid.icel(b0_g[1], j, k)] = old_vel.v[grid.icel(b0_g[1], j, k)];

      curr_vel.w[grid.icel(b0_g[0], j, k)] = old_vel.w[grid.icel(b0_g[0], j, k)];
      curr_vel.w[grid.icel(b0_g[1], j, k)] = old_vel.w[grid.icel(b0_g[1], j, k)];

      curr_vel.u[grid.icel(b1_g[1], j, k)] =
          old_vel.u[grid.icel(b1_g[1], j, k)];  // omission?
      curr_vel.u[grid.icel(b1_c[0], j, k)] = old_vel.u[grid.icel(b1_c[0], j, k)];
      curr_vel.u[grid.icel(b1_g[0], j, k)] = old_vel.u[grid.icel(b1_g[0], j, k)];

      curr_vel.v[grid.icel(b1_g[0], j, k)] = old_vel.v[grid.icel(b1_g[0], j, k)];
      curr_vel.v[grid.icel(b1_g[1], j, k)] = old_vel.v[grid.icel(b1_g[1], j, k)];

      curr_vel.w[grid.icel(b1_g[0], j, k)] = old_vel.w[grid.icel(b1_g[0], j, k)];
      curr_vel.w[grid.icel(b1_g[1], j, k)] = old_vel.w[grid.icel(b1_g[1], j, k)];
    }

  for (int i = location.x_start - d; i < location.x_end + d; ++i)
    for (int k = location.z_start - d; k < location.z_end + d;
         ++k) {  // ====================== BC_2 and BC_3 ======================
      curr_vel.u[grid.icel(i, b2_g[0], k)] = old_vel.u[grid.icel(i, b2_g[0], k)];
      curr_vel.u[grid.icel(i, b2_g[1], k)] = old_vel.u[grid.icel(i, b2_g[1], k)];

      curr_vel.v[grid.icel(i, b2_g[0], k)] = old_vel.v[grid.icel(i, b2_g[0], k)];
      curr_vel.v[grid.icel(i, b2_g[1], k)] = old_vel.v[grid.icel(i, b2_g[1], k)];

      curr_vel.w[grid.icel(i, b2_g[0], k)] = old_vel.w[grid.icel(i, b2_g[0], k)];
      curr_vel.w[grid.icel(i, b2_g[1], k)] = old_vel.w[grid.icel(i, b2_g[1], k)];

      curr_vel.u[grid.icel(i, b3_g[0], k)] = old_vel.u[grid.icel(i, b3_g[0], k)];
      curr_vel.u[grid.icel(i, b3_g[1], k)] = old_vel.u[grid.icel(i, b3_g[1], k)];

      curr_vel.v[grid.icel(i, b3_c[0], k)] = old_vel.v[grid.icel(i, b3_c[0], k)];
      curr_vel.v[grid.icel(i, b3_g[0], k)] = old_vel.v[grid.icel(i, b3_g[0], k)];
      curr_vel.v[grid.icel(i, b3_g[1], k)] =
          old_vel.v[grid.icel(i, b3_g[1], k)];  // omission?

      curr_vel.w[grid.icel(i, b3_g[0], k)] = old_vel.w[grid.icel(i, b3_g[0], k)];
      curr_vel.w[grid.icel(i, b3_g[1], k)] = old_vel.w[grid.icel(i, b3_g[1], k)];
    }

  for (int i = location.x_start - d; i < location.x_end + d; ++i)
    for (int j = location.y_start - d; j < location.y_end + d;
         ++j) {  // ====================== BC_4 and BC_5 ======================
      curr_vel.u[grid.icel(i, j, b4_g[0])] = old_vel.u[grid.icel(i, j, b4_g[0])];
      curr_vel.u[grid.icel(i, j, b4_g[1])] = old_vel.u[grid.icel(i, j, b4_g[1])];

      curr_vel.v[grid.icel(i, j, b4_g[0])] = old_vel.v[grid.icel(i, j, b4_g[0])];
      curr_vel.v[grid.icel(i, j, b4_g[1])] = old_vel.v[grid.icel(i, j, b4_g[1])];

      curr_vel.w[grid.icel(i, j, b4_g[0])] = old_vel.w[grid.icel(i, j, b4_g[0])];
      curr_vel.w[grid.icel(i, j, b4_g[1])] = old_vel.w[grid.icel(i, j, b4_g[1])];

      curr_vel.u[grid.icel(i, j, b5_g[0])] = old_vel.u[grid.icel(i, j, b5_g[0])];
      curr_vel.u[grid.icel(i, j, b5_g[1])] = old_vel.u[grid.icel(i, j, b5_g[1])];

      curr_vel.v[grid.icel(i, j, b5_g[0])] = old_vel.v[grid.icel(i, j, b5_g[0])];
      curr_vel.v[grid.icel(i, j, b5_g[1])] = old_vel.v[grid.icel(i, j, b5_g[1])];

      curr_vel.w[grid.icel(i, j, b5_c[0])] = old_vel.w[grid.icel(i, j, b5_c[0])];
      curr_vel.w[grid.icel(i, j, b5_g[0])] = old_vel.w[grid.icel(i, j, b5_g[0])];
      curr_vel.w[grid.icel(i, j, b5_g[1])] =
          old_vel.w[grid.icel(i, j, b5_g[1])];  // omission?
    }

  return true;
}

bool BC_updateSlid(const CalDomain& location, StaggeredVelocity& vel,
                   StructuredGrid& grid) {
  auto ii = [&](auto& i, auto& j, auto& k) { return grid.icel(i, j, k); };

  // auto [num, dir] = getNumDIr();

  auto num = grid.bc_selector.num;
  auto dir = grid.bc_selector.dir;



  std::vector<int> b0_g = {1, 0}, b0_c = {2, 3, 4},

                   b2_g = {1, 0}, b2_c = {2, 3, 4},

                   b4_g = {1, 0}, b4_c = {2, 3, 4};
  int d = 1;  // o

  for (int j = location.y_start - d; j < location.y_end + d; ++j)
    for(int k = location.z_start-d ; k < location.z_end+d; ++k)
    {// ====================== BC_0 ======================
        const int b = 0;
        vel.u[grid.icel(b0_g[1],j,k)] = vel.u[grid.icel(b0_g[0],j,k)] =  dir[b][0]    +      num[b][0]      * vel.u[grid.icel(b0_c[0],j,k)];
        vel.v[grid.icel(b0_g[1],j,k)] = vel.v[grid.icel(b0_g[0],j,k)] =  dir[b][1]*2. + (2.0*num[b][1]-1.0) * vel.v[grid.icel(b0_c[0],j,k)];
        vel.w[grid.icel(b0_g[1],j,k)] = vel.w[grid.icel(b0_g[0],j,k)] =  dir[b][2]*2. + (2.0*num[b][2]-1.0) * vel.w[grid.icel(b0_c[0],j,k)];
    }

    for(int i = location.x_start-d ; i < location.x_end+d; ++i )
    for(int k = location.z_start-d ; k < location.z_end+d; ++k )
    {// ====================== BC_2 ======================
        const int b = 2;
        vel.u[grid.icel(i,b2_g[1],k)] = vel.u[grid.icel(i,b2_g[0],k)] =  dir[b][0]*2. + (2.0*num[b][0]-1.0) * vel.u[grid.icel(i,b2_c[0],k)];
        vel.v[grid.icel(i,b2_g[1],k)] = vel.v[grid.icel(i,b2_g[0],k)] =  dir[b][1]    +      num[b][1]      * vel.v[grid.icel(i,b2_c[0],k)];
        vel.w[grid.icel(i,b2_g[1],k)] = vel.w[grid.icel(i,b2_g[0],k)] =  dir[b][2]*2. + (2.0*num[b][2]-1.0) * vel.w[grid.icel(i,b2_c[0],k)];
    }

    for(int i = location.x_start-d ; i < location.x_end+d; ++i )
    for(int j = location.y_start-d ; j < location.y_end+d; ++j )
    {// ====================== BC_4 ======================
        const int b = 4;
        vel.u[grid.icel(i,j,b4_g[1])] = vel.u[grid.icel(i,j,b4_g[0])] =  dir[b][0]*2. + (2.0*num[b][0]-1.0) * vel.u[grid.icel(i,j,b4_c[0])];
        vel.v[grid.icel(i,j,b4_g[1])] = vel.v[grid.icel(i,j,b4_g[0])] =  dir[b][1]*2. + (2.0*num[b][1]-1.0) * vel.v[grid.icel(i,j,b4_c[0])];
        vel.w[grid.icel(i,j,b4_g[1])] = vel.w[grid.icel(i,j,b4_g[0])] =  dir[b][2]    +      num[b][2]      * vel.w[grid.icel(i,j,b4_c[0])];
    }
    return true;
}
