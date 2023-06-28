#pragma once
#include <dirent.h>
#include <sys/types.h>

#include "backend/physical_variables.hpp"
#include "backend/simulation.hpp"
#include "grid/structured_grid.hpp"
#include <string>
#include <fstream>

auto ReadQfile(StructuredGrid &grid, std::string filename) {
  //*  prepare the data to read
  float mach, alpha, reyn, time;
  int Nblock, nxCal, nyCal, nzCal;
  auto cal_nz = grid.cal_nz;
  auto cal_ny = grid.cal_ny;
  auto cal_nx = grid.cal_nx;
  auto gc = grid.no_ghost_cell;

  const uint nodata = 5;
  std::vector<std::vector<float> > data(
      nodata, std::vector<float>(grid.cal_no_grid, 0.0));
  // * read data;

  std::ifstream file;
  file.open(filename, std::ifstream::binary);

  auto ii = [&](auto i, auto j, auto k) {
    return i * nzCal * nyCal + j * nzCal + k;
  };

  if (!file.is_open()) {
    file.close();
    throw std::invalid_argument("can't Read The File !\n");
  }


  file.read((char *)(&Nblock), sizeof(int));
  file.read((char *)(&nxCal), sizeof(int));
  file.read((char *)(&nyCal), sizeof(int));
  file.read((char *)(&nzCal), sizeof(int));

  file.read((char *)(&mach), sizeof(float));
  file.read((char *)(&alpha), sizeof(float));
  file.read((char *)(&reyn), sizeof(float));
  file.read((char *)(&time), sizeof(float));

  for (auto &x : data) {
    for (size_t k = 0; k < nzCal; ++k)
      for (size_t j = 0; j < nyCal; ++j)
        for (size_t i = 0; i < nxCal; ++i)
          file.read((char *)(&x[ii(i, j, k)]), sizeof(float));
  }

  file.close();
  // float mach, alpha, reyn, time;

  return std::make_tuple(Nblock, mach, alpha, reyn, time, data);
}
