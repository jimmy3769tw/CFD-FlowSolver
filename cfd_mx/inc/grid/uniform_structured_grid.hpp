
#pragma once
#include "structured_grid.hpp"

class UniformStructuredGrid : public StructuredGrid {
  public:
   UniformStructuredGrid(int cal_nx, int cal_ny, int cal_nz)
       : StructuredGrid(cal_nx, cal_ny, cal_nz) {
   }

  private:
   virtual void InitPos() override {
     auto dx = lx_ / (cal_nx);
     auto dy = ly_ / (cal_ny);
     auto dz = lz_ / (cal_nz);
     x_pos[0] = -2.0 * dx;
     for (size_t i = 1; i < nx + 1; ++i) {
       x_pos.at(i) = x_pos.at(i - 1) + dx;
     }

     y_pos[0] = -2.0 * dy;
     for (size_t j = 1; j < ny + 1; ++j) {
       y_pos.at(j) = y_pos.at(j - 1) + dy;
     }

     z_pos[0] = -2.0 * dz;
     for (size_t k = 1; k < nz + 1; ++k) {
       z_pos.at(k) = z_pos.at(k - 1) + dz;
     }
   }
};