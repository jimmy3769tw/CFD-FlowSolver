#pragma once

#include <iostream>
#include <tuple>
#include <vector>

class StructuredGrid;

class GlobalDomain {
 public:
 
  GlobalDomain(StructuredGrid& grid) : grid_(&grid) {}
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

 void DivideDomain(int dims, int n);
  private:
  StructuredGrid* grid_;
};

class LocalDomain{
 public:
  LocalDomain(StructuredGrid& grid);

  int x_start, y_start, z_start, x_end, y_end, z_end;

  void Init(std::vector<int> dims, std::vector<int> rank);

  std::pair<std::vector<int>, std::vector<int> > 
  GetOneDimStartEnd();

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

// void LocalDomain::show_begin(){
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

