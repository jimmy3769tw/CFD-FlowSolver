#pragma once
#include <iostream>
#include <tuple>
#include <vector>
#include "../grid/structured_grid.hpp"
#include "global_domain.hpp"

class CalDomain{
 public:
  CalDomain(StructuredGrid &grid)
      : grid_(&grid),
        global_x_(GlobalDomain(grid)),
        global_y_(GlobalDomain(grid)),
        global_z_(GlobalDomain(grid)) {}
  int x_start, y_start, z_start, x_end, y_end, z_end;

  void Init(std::vector<int> dims, std::vector<int> rank) {
    global_x_.DivideDomain(dims[0], grid_->cal_nx);
    global_y_.DivideDomain(dims[1], grid_->cal_ny);
    global_z_.DivideDomain(dims[2], grid_->cal_nz);
    InitLocal(rank[0], rank[1], rank[2]);
  }

  std::pair<std::vector<int>, std::vector<int> > GetOneDimStartEnd() {
    int slid = grid_->cal_ny_times_nz;
    auto beg = global_x_.start;
    auto end = global_x_.end;

    for (auto &x : beg) {
      x -= grid_->no_ghost_cell;
    }
    for (auto &x : end) {
      x -= grid_->no_ghost_cell;
    }
    for (auto &x : beg) {
      x *= slid;
    }
    for (auto &x : end) {
      x *= slid;
    }

    return std::make_pair(beg, end);
  }

  GlobalDomain GetGlobal_x() { return global_x_; }

 private:
  StructuredGrid *grid_;
  GlobalDomain global_x_;
  GlobalDomain global_y_;
  GlobalDomain global_z_;
  void InitLocal(int x_rank, int y_rank, int z_rank) {
    x_start = global_x_.start.at(x_rank);
    y_start = global_y_.start.at(y_rank);
    z_start = global_z_.start.at(z_rank);
    x_end = global_x_.end.at(x_rank);
    y_end = global_y_.end.at(y_rank);
    z_end = global_z_.end.at(z_rank);
  }
};

// void CalDomain::show_begin(){
// 	cout << "show begin";
// 	for (auto x :i_begin_table){
// 		std::cout << x << ", ";
// 	}

// 	for (auto x :j_begin_table){
// 		std::cout << x << ", ";
// 	}

// 	for (auto x :k_begin_table){
// 		std::cout << x << ", ";
// 	}
// }

