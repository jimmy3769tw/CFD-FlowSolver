#pragma once
#include "mpi.h"
#include "mat_mpi_tool.hpp"

namespace mat{
template <typename T>
class SparseMat {
 public:
  using CSR_type = 
      // <ptr, indices, values>
      std::tuple< std::vector<int>, std::vector<int>, std::vector<T> >;

  SparseMat(int rows, int cols) { Construct(rows, cols); }
  SparseMat(int n) { Construct(n, n); }
  SparseMat(){}

  virtual void Multiply(std::vector<T>& x, std::vector<T>& r) = 0;
  void multiply_omp(std::vector<T>& x, std::vector<T>& r) { Multiply(x, r); }

  int row() const { return rows_; }
  int col() const { return cols_; }

  std::vector<T> Multiply(std::vector<T>& x) {
    re_global_.resize(rows_);
    Multiply(x, re_global_);
    return re_global_;
  }
  
  std::vector<T> operator*(std::vector<T>& x) {
    return Multiply(x);
  }

  void Construct(int rows, int cols) {
    if (rows < 0 || cols < 0) {
      throw std::invalid_argument("Matrix dimensions cannot be negative.");
    }
    rows_ = rows;
    cols_ = cols;
  }

  void InitMpi(MPI_Comm comm_world, 
              std::vector<int> beg,
              std::vector<int> end) {
    mpi_.Init(comm_world, beg, end);
  }

  void InitMpi(MPI_Comm comm_world, int vec_size) {
    mpi_.Init(comm_world, vec_size);
  }


  MpiVectorTool mpi_;
 private:
  int rows_, cols_;
  std::vector<T> re_global_;
};
}