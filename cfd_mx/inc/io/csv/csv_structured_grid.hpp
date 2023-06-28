#pragma once
#include <fstream>
#include <string>
#include <vector>

#include "../../grid/structured_grid.hpp"
namespace csv{
bool WriteCsvFile(std::string filename, StructuredGrid& grid) {
  std::ofstream file;
  auto nx = grid.nx;
  auto ny = grid.ny;
  auto nz = grid.nz;
  auto no_ghost_cell = grid.no_ghost_cell;
  auto dy = grid.dy;
  auto dx = grid.dx;
  auto dz = grid.dz;
  auto y_cent_pos = grid.y_cent_pos;
  auto x_cent_pos = grid.x_cent_pos;
  auto z_cent_pos = grid.z_cent_pos;
  auto staggered_dx = grid.staggered_dx;
  auto staggered_dy = grid.staggered_dy;
  auto staggered_dz = grid.staggered_dz;
  auto x_pos = grid.x_pos;
  auto y_pos = grid.y_pos;
  auto z_pos = grid.z_pos;

  file.open(filename, std::ios::out);  // | ios::app
  // * --------------------- x---------------------
  for (size_t i = 0; i < nx + 1; ++i) {
    file << "x_pos[" << i << "] ,";
    file << x_pos[i] << ", ";
    if (i < nx - no_ghost_cell + 2) {
      if (i < 2) {
        file << "x_cent_pos[no_ghost_cell], not define, ";
      } else if (i < nx - no_ghost_cell) {
        file << "xc[" << i - 2 << "] ," << x_cent_pos[i - 2] << ", ";
      } else {
        file << "x_cent_pos[no_ghost_cell], not define, ";
      }
    } else {
      file << " , , ";
    }
    if (i < nx) {
      file << "dx[" << i << "] ," << dx[i] << ", "
           << "staggered_dx[" << i << "] ," << staggered_dx[i];

    } else {
      file << " , , ";
    }

    file << std::endl;
  }

  file << std::endl;
  // * --------------------- y---------------------
  for (size_t j = 0; j < ny + 1; ++j) {
    file << "y_pos[" << j << "] ,";
    file << y_pos[j] << ", ";

    if (j < ny - no_ghost_cell + 2) {
      if (j < 2) {
        file << "y_cent_pos[no_ghost_cell], not define, ";
      } else if (j < ny - no_ghost_cell) {
        file << "y_cent_pos[" << j - 2 << "] ," << y_cent_pos[j - 2] << ", ";
      } else {
        file << "y_cent_pos[no_ghost_cell], not define, ";
      }

    } else {
      file << " , , ";
    }

    if (j < ny) {
      file << "dy[" << j << "] ," << dy[j] << ", "
           << "staggered_dy[" << j << "] ," << staggered_dy[j];
    } else {
      file << " , , ";
    }

    file << std::endl;
  }

  file << std::endl;
  // *---------------------------- z --------------------
  for (size_t k = 0; k < nz + 1; ++k) {
    file << "z_pos[" << k << "] ," << z_pos[k] << ", ";

    if (k < nz - no_ghost_cell + 2) {
      if (k > 1 || k < nz - no_ghost_cell) {
        file << "zc[" << k - 2 << "] ," << z_cent_pos[k - 2] << ", ";
      } else {
        file << "dzs[no_ghost_cell], not define,";
      }

    } else {
      file << " , , ";
    }

    if (k < nz) {
      file << "dz[" << k << "] ," << dz[k] << ", "
           << "dzs[" << k << "] ," << staggered_dz[k];
    } else {
      file << " , , ";
    }

    file << std::endl;
  }
  return true;
}
}

