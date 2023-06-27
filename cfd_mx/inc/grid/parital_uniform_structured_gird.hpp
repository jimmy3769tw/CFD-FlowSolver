
#pragma once
#include "structured_grid.hpp"

// UniformStructuredGrid
class PartialUniformStructuredGrid : public StructuredGrid {
 public:
  PartialUniformStructuredGrid(int cal_nx, int cal_ny, int cal_nz)
      : StructuredGrid(cal_nx, cal_ny, cal_nz) {}

 private:
  virtual void InitPos() override {
    const double lx_small = 0.3 * lx_, ly_small = 0.3 * ly_, lz_small = 0.5 * lz_;
    const int nx_small = nx - (0.5 * nx / 2 * 2),
              ny_small = ny - (0.5 * ny / 2 * 2), nz_small = 0;
    const int offset_nx_small = -4, offset_ny_small = 0, offset_nz_small = 0;

    PartialUniformConstruct(x_pos, lx_small, nx_small, offset_nx_small, lx_, nx);
    PartialUniformConstruct(y_pos, ly_small, ny_small, offset_ny_small, ly_, ny);
    PartialUniformConstruct(z_pos, lz_small, nz_small, offset_nz_small, lz_, nz);
  }

  void PartialUniformConstruct(vector<double> &pos, const double small_len,
                               int small_grid, int offset_small_grid,
                               double total_len, int total_grid) {
    auto cal_grid = total_grid - 2 * this->no_ghost_cell;

    pos.at(this->no_ghost_cell) = 0.0;
    if (small_grid == 0) {
      for (size_t i = this->no_ghost_cell; i < total_grid - this->no_ghost_cell;
           i++) {
        pos.at(i + 1) = pos.at(i) + total_len / cal_grid;
      }
      return;
    }

    if (small_grid != 0) {
      if ((total_grid - small_grid) % 2 != 0) {
        small_grid++;
        // throw std::invalid_argument("(total_grid-n_sml)%2 != 0");
      }

      const double diff_small = small_len / small_grid;
      const double diff = (total_len - small_len) / (cal_grid - small_grid);

      int domain[3];
      domain[0] = (cal_grid - small_grid) / 2;
      domain[1] = domain[0] + small_grid;
      domain[2] = domain[1] + domain[0];

      domain[0] += offset_small_grid;
      domain[1] += offset_small_grid;

      if (domain[0] <= 0) {
        throw std::invalid_argument("cal_grid <= small_grid");
      } 

      for (size_t i = 0; i < domain[0]; i++)
        pos.at(i + this->no_ghost_cell + 1) =
            pos.at(i + this->no_ghost_cell) + diff;

      for (size_t i = domain[0]; i < domain[1]; ++i)
        pos.at(i + this->no_ghost_cell + 1) =
            pos.at(i + this->no_ghost_cell) + diff_small;

      for (size_t i = domain[1]; i < domain[2]; i++)
        pos.at(i + this->no_ghost_cell + 1) =
            pos.at(i + this->no_ghost_cell) + diff;

    }
  }
};
