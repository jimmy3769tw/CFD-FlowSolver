#pragma once
#include <string>
#include <sys/types.h>
#include <dirent.h>
#include "../../grid/structured_grid.hpp"
#include "../../backend/physical_variables.hpp"
#include "../../backend/simulation.hpp"

void WriteQfile(ImmersedBoundary &dfib, Simulation &simu, Pressure &pressure,
                StaggeredVelocity &curr_vel, StructuredGrid &grid) {
  auto gc = grid.no_ghost_cell;
  auto nx = grid.nx;
  auto ny = grid.ny;
  auto nz = grid.nz;
  auto cal_nx = grid.cal_nx;
  auto cal_ny = grid.cal_ny;
  auto cal_nz = grid.cal_nz;

  // double temp = 1.0;
  float mach = 0;
  float alpha = 0;
  float reyn = simu.GetReynoldsNumber();
  float time = simu.tva.GetTime();
  int io = simu.tva.GetFileNo();

  std::ofstream file;

  int no_block = 1;
  int temp_nx = grid.cal_nx;
  int temp_ny = grid.cal_ny;
  int temp_nz = grid.cal_nz;

  string filename = "mx_out/P3D";

  if (io < 10) {
    filename += "0";
  }

  filename += std::to_string(io);
  filename += ".q";
  file.open(filename, ofstream::binary);
  cout << "(Output Plaot3D):" << io << endl;

  file.write((char *)(&no_block), sizeof(int));
  file.write((char *)(&temp_nx), sizeof(int));
  file.write((char *)(&temp_ny), sizeof(int));
  file.write((char *)(&temp_nz), sizeof(int));

  file.write((char *)(&mach), sizeof(float));
  file.write((char *)(&alpha), sizeof(float));
  file.write((char *)(&reyn), sizeof(float));
  file.write((char *)(&time), sizeof(float));

  auto ii = [&](auto i, auto j, auto k) { return grid.icel(i, j, k); };

  for (size_t k = gc; k < cal_nz + gc; ++k)
    for (size_t j = gc; j < cal_ny + gc; ++j)
      for (size_t i = gc; i < cal_nx + gc; ++i) {
        float out = pressure.p[ii(i, j, k)];
        file.write((char *)(&out), sizeof(float));
      }

  for (size_t k = gc; k < cal_nz + gc; ++k)
    for (size_t j = gc; j < cal_ny + gc; ++j)
      for (size_t i = gc; i < cal_nx + gc; ++i) {
        float out = 0.5 * (curr_vel.u[ii(i - 1, j, k)] + curr_vel.u[ii(i, j, k)]);
        file.write((char *)(&out), sizeof(float));
      }

  for (size_t k = gc; k < cal_nz + gc; ++k)
    for (size_t j = gc; j < cal_ny + gc; ++j)
      for (size_t i = gc; i < cal_nx + gc; ++i) {
        float out = 0.5 * (curr_vel.v[ii(i, j - 1, k)] + curr_vel.v[ii(i, j, k)]);
        file.write((char *)(&out), sizeof(float));
      }
  for (size_t k = gc; k < cal_nz + gc; ++k)
    for (size_t j = gc; j < cal_ny + gc; ++j)
      for (size_t i = gc; i < cal_nx + gc; ++i) {
        float out = 0.5 * (curr_vel.w[ii(i, j, k - 1)] + curr_vel.w[ii(i, j, k)]);
        file.write((char *)(&out), sizeof(float));
      }

  for (size_t k = gc; k < cal_nz + gc; ++k)
    for (size_t j = gc; j < cal_ny + gc; ++j)
      for (size_t i = gc; i < cal_nx + gc; ++i) {
        float out = dfib.eta[ii(i, j, k)];
        file.write((char *)(&out), sizeof(float));
      }
}
