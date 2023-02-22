#pragma once
#include "../grid/structured_grid.hpp"

class GlobalDomain {
 public:
  GlobalDomain(StructuredGrid &grid) : grid_(&grid) {}

  std::vector<int> start;
  std::vector<int> end;
  std::vector<int> len;

  void resize(int size) {
    len.resize(size);
    end.resize(size);
    start.resize(size);
  }

  void DivideLen(int dims, int n) {
    for (int i = 0; i < dims; ++i) {
      len[i] = (n - n % dims) / (dims);
    }

    for (int i = 0; i < n % dims; ++i) {
      len[i] += 1;
    }
  }

  void DivideDomain(int dims, int n) {
    if (dims == 0) {
      throw std::invalid_argument("DivideDomain !!!!\n");
    }
    resize(dims);
    DivideLen(dims, n);
    start[0] = grid_->no_ghost_cell;
    for (int i = 1; i < dims; i++) {
      start[i] = start[i - 1] + len[i - 1];
    }

    for (int i = 0; i < dims; i++) {
      end[i] = start[i] + len[i];
    }
  }

  private:
  StructuredGrid* grid_;
};