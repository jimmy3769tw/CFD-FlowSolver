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
  MpiVectorTool(){
    empty_ = true;
  }
  inline int beg() { return beg_; }
  inline int end() { return end_; }
  void Barrier() { MPI_Barrier(comm_world_); }

  void Init(MPI_Comm comm_world, std::vector<int> start,
            std::vector<int> end) {
    empty_ = false;
    comm_world_ = comm_world;
    rank_ = mpi::GetRank(comm_world_);
    size_ = mpi::GetSize(comm_world_);
    InitLocal(start, end);
    InitGlobal(start, end);
  }

  void Init(MPI_Comm &comm_world, int vec_size) {
    size_ = mpi::GetSize(comm_world);
    std::vector<int> len(size_);
    std::vector<int> start(size_);
    std::vector<int> end(size_);
    cout << "size_" << size_;
    for (int i = 0; i < size_; i++) {
      len[i] = vec_size / size_;
    }

    for (int i = 0; i < vec_size % size_; i++) {
        len[i]++;
    }

    start[0] = 0;
    for (int i = 1; i < size_; i++) {
        start[i] = start[i - 1] + len[i - 1];
    }

    for (int i = 0; i < size_; i++) {
        end[i] = start[i] + len[i];
    }
    Init(comm_world, start, end);
  }


  bool Empty() { return empty_;}
  void Allocate(std::vector<double> &x) {
    for (size_t i = 0; i < size_; ++i) {
      MPI_Bcast((void *)&x[table_beg_[i]], table_len_[i], MPI_DOUBLE, i, comm_world_);
    }
    MPI_Barrier(comm_world_);
  }

  double AllreduceSum(double &x) {
    double temp_g = x;
    MPI_Allreduce(&x, &temp_g, 1, MPI_DOUBLE, MPI_SUM, comm_world_);
    return temp_g;
  }

  int AllreduceSum(int &x) {
    int temp_g = x;
    MPI_Allreduce(&x, &temp_g, 1, MPI_INT, MPI_SUM, comm_world_);
    return temp_g;
  }

  inline double InnerProduct(const vector<double> &a, const vector<double> &b) {
    double temp_l{0.0};
    auto A = a.data() + beg_;
    auto B = b.data() + beg_;

#pragma omp parallel for simd reduction(+:temp_l) default(none) \
            firstprivate(A, B, len_) schedule(static)
    for (size_t i = 0; i < len_; ++i) {
      temp_l += *(A + i) * *(B + i);
    }
    return AllreduceSum(temp_l);
  }

  inline double L2Norm(const std::vector<double> &a) { return sqrt(InnerProduct(a, a)); }

  inline void init(vector<double> &a, const double val) {
    auto A = a.data() + beg_;
#pragma omp parallel for simd default(none) firstprivate(A, val, len_) \
    schedule(static)
    for (size_t i = 0; i < len_; ++i) {
      *(A + i) = val;
    }
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
  bool empty_ = true;
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
