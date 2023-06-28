#pragma once
#include "set_eta.hpp"
#include <vector>

class CylinderZ : public CalEtaInterface {
 public:

  CylinderZ(StructuredGrid& grid) : grid_(&grid), CalEtaInterface() {
    std::tie(lx, ly, lz) = grid_->GetLen();
  }

  void solver(ImmersedBoundary& dfib, LocalDomain& domain) override;

  CylinderZ& SetCylinderCenter(double center_x, double center_y) {
    this->center_x = center_x;
    this->center_y = center_y;
    return *this;
  }

  CylinderZ& SetSubGrid(int nx, int ny) {
    no_sub_nx = nx;
    no_sub_ny = ny;
    return *this;
  }

  CylinderZ& SetRadius(double rad) { 
    radius = rad;
    return *this;
  }

 private:
  StructuredGrid* grid_;
  int no_sub_nx = 100;
  int no_sub_ny = 100;
  int ly;
  int lx;
  int lz;

  double center_x{0};
  double center_y{0};
  double radius{0};
  int gc = grid_->no_ghost_cell;
  double nx = grid_->nx;
  double ny = grid_->ny;
  double nz = grid_->nz;
};
