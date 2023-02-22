#pragma once
#include <cstdlib>  // for div
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include "../mpi_tool/mpi_getter.hpp"
#include "mpi.h"
using namespace std;
namespace mat {

class MpiVectorTool {
 public: 

  inline int beg() { return beg_; }
  inline int end() { return end_; }
  void Barrier() { MPI_Barrier(comm_world_); }
  void Init(MPI_Comm &comm_world, std::vector<int> start,
            std::vector<int> end) {
    comm_world_ = comm_world;
    rank_ = mpi::GetRank(comm_world_);
    size_ = mpi::GetSize(comm_world_);
    InitLocal(start, end);
    InitGlobal(start, end);
  }

  void Allocate(std::vector<double> &x) {
    MPI_Barrier(comm_world_);
    for (size_t i = 0; i < size_; ++i) {
      MPI_Bcast((void *)&x[table_beg_[i]], table_len_[i], MPI_DOUBLE, i, comm_world_);
      MPI_Barrier(comm_world_);
    }
    MPI_Barrier(comm_world_);
  }

  double AllreduceSum(double &x) {
    MPI_Barrier(comm_world_);
    double temp_g = x;
    MPI_Allreduce(&x, &temp_g, 1, MPI_DOUBLE, MPI_SUM, comm_world_);
    return temp_g;
    MPI_Barrier(comm_world_);
  }

  int AllreduceSum(int &x) {
    MPI_Barrier(comm_world_);
    int temp_g = x;
    MPI_Allreduce(&x, &temp_g, 1, MPI_INT, MPI_SUM, comm_world_);
    return temp_g;
    MPI_Barrier(comm_world_);
  }

  inline double InnerProduct(const vector<double> &a, const vector<double> &b) {
    MPI_Barrier(comm_world_);
    double temp_l{0.0};
    auto A = a.data() + beg_;
    auto B = b.data() + beg_;

#pragma omp parallel for simd reduction(+:temp_l) default(none) \
            firstprivate(A, B, len_) schedule(static)
    for (size_t i = 0; i < len_; ++i) {
      temp_l += *(A + i) * *(B + i);
    }
    MPI_Barrier(comm_world_);
    return AllreduceSum(temp_l);
  }

  inline double L2Norm(const std::vector<double> &a) { return sqrt(InnerProduct(a, a)); }

  inline void init(vector<double> &a, const double val) {
    auto A = a.data() + beg_;
    MPI_Barrier(comm_world_);
#pragma omp parallel for simd default(none) firstprivate(A, val, len_) \
    schedule(static)
    for (size_t i = 0; i < len_; ++i) {
      *(A + i) = val;
    }
    MPI_Barrier(comm_world_);
  }

  inline void init(vector<double> &a, const vector<double> &b) {
    auto A = a.data() + beg_;
    auto B = b.data() + beg_;

#pragma omp parallel for simd default(none) firstprivate(A, len_, B) \
    schedule(static)
    for (size_t i = 0; i < len_; ++i) {
      *(A + i) = *(B + i);
    }
  }

  inline void zero(vector<double> &a) { init(a, 0.0); }

 private:
  MPI_Comm comm_world_;
  int rank_, size_;
  int beg_, end_, len_;

  std::vector<int> table_end_, table_beg_, table_len_;

  void InitLocal(std::vector<int> start, std::vector<int> end) {
    beg_ = start.at(rank_);
    end_ = end.at(rank_);
    len_ = end_ - beg_;
  }

  void InitGlobal(std::vector<int> start, std::vector<int> end) {
    table_beg_ = start;
    table_end_ = end;
    table_len_ = end;
    for (int i = 0; i < end.size(); ++i) {
        table_len_[i] -= start[i];
    }

  }
};



}  // namespace mat
