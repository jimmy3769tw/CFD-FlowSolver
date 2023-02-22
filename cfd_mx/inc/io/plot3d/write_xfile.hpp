#pragma once
#include "../../grid/structured_grid.hpp"
namespace plot3d {

auto write_xfile(StructuredGrid& grid) {
   vector<vector<vector<float> > > x_out(
      grid.cal_nx, vector<vector<float> >(grid.cal_ny, vector<float>(grid.cal_nz)));
  vector<vector<vector<float> > > y_out(
      grid.cal_nx, vector<vector<float> >(grid.cal_ny, vector<float>(grid.cal_nz)));
  vector<vector<vector<float> > > z_out(
      grid.cal_nx, vector<vector<float> >(grid.cal_ny, vector<float>(grid.cal_nz)));

  int Nblock = 1, tempNX = grid.cal_nx, tempNY = grid.cal_ny, tempNZ = grid.cal_nz;

  for (size_t i = 0; i <  grid.cal_nx; ++i)
    for (size_t j = 0; j <  grid.cal_ny; ++j)
      for (size_t k = 0; k <  grid.cal_nz; ++k) {
        x_out[i][j][k] = grid.x_cent_pos.at(i);
        y_out[i][j][k] = grid.y_cent_pos.at(j);
        z_out[i][j][k] = grid.z_cent_pos.at(k);
      }

  FILE* fptr;

  fptr = fopen("mx_out/P3D.x", "wb");

  fwrite(&Nblock, sizeof(int), 1, fptr);
  fwrite(&tempNX, sizeof(int), 1, fptr);
  fwrite(&tempNY, sizeof(int), 1, fptr);
  fwrite(&tempNZ, sizeof(int), 1, fptr);

  for (size_t k = 0; k < grid.cal_nz; ++k)
    for (size_t j = 0; j < grid.cal_ny; ++j)
      for (size_t i = 0; i < grid.cal_nx; ++i)
        fwrite(&x_out[i][j][k], sizeof(float), 1, fptr);

  for (size_t k = 0; k < grid.cal_nz; ++k)
    for (size_t j = 0; j < grid.cal_ny; ++j)
      for (size_t i = 0; i < grid.cal_nx; ++i)
        fwrite(&y_out[i][j][k], sizeof(float), 1, fptr);

  for (size_t k = 0; k < grid.cal_nz; ++k)
    for (size_t j = 0; j < grid.cal_ny; ++j)
      for (size_t i = 0; i < grid.cal_nx; ++i)
        fwrite(&z_out[i][j][k], sizeof(float), 1, fptr);

  fclose(fptr);
  return 0;
}
}

