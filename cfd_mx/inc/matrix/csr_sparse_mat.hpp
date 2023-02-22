#pragma once
#include <algorithm>  // (until C++11)
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>  // (since C++11) std::pair
#include <vector>
#include "sparse_mat.hpp"

using std::vector;
using std::cout;
using std::endl;
using std::tuple;

namespace mat {

template <typename T>
class CsrMat : public SparseMat<T> {
 public:
  using CsrVectorType = std::tuple<std::vector<int>, std::vector<int>,
                              std::vector<T> >;  // <ptr, indices, values>

  CsrMat() : SparseMat<T>() {}

  CsrMat(int rows, int cols) : SparseMat<T>(rows, cols) { SelfConstruct(); }

  CsrMat(int n) : SparseMat<T>(n) { SelfConstruct(); }

  CsrMat(CsrVectorType c) { Set(c); }

  void resize(int rows, int cols) {
    this->Construct(rows, cols);
    this->SelfConstruct();
  }

  void resize(int n) {
    this->Construct(n, n);
    SelfConstruct();
  }

  void Multiply(std::vector<T>& x, std::vector<T>& r) override {
#pragma omp parallel for shared(r, x)
        for (int i = 0; i < this->row(); ++i) {
      r[i] = T();
      for (int j = ptr_[i]; j < ptr_[i + 1]; ++j) {
        r[i] += val_.at(j) * x.at(idx_.at(j));
      }
    }
  }

  void multiply_mpi(std::vector<T>& x, std::vector<T>& r) {
    this->mpi_.Allocate(x);
#pragma omp parallel for shared(r, x)
    for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i) {
      r[i] = T();
      for (int j = ptr_[i]; j < ptr_[i + 1]; ++j) {
        r[i] += val_[j] * x[idx_[j]];
      }
    }
  }

  void init(int Idx_len, int ptr_len) {
    idx_.clear();
    val_.clear();
    ptr_.clear();

    idx_.reserve(Idx_len);
    val_.reserve(Idx_len);
    ptr_.reserve(ptr_len);

    ptr_.push_back(0);
  }

  void Set(const CsrVectorType& csr) {
    std::tie(ptr_, idx_, val_) = csr;
    auto len = ptr_.size() - 1;
    this->resize(len, len);
    // CSR can't not received the information regarding the number
    // of cols! So, it assume the cols_ = row;
  }

  inline void Show() {
    std::cout << "\n|prt\n";
    for (auto v : ptr_) {
      std::cout << v << ", ";
    }

    std::cout << "\n|Idx\n";
    for (auto v : idx_) {
      std::cout << v << ", ";
    }

    std::cout << "\n|val\n";
    for (auto v : val_) {
      std::cout << v << ", ";
    }

    std::cout << std::endl;
  }

  CsrVectorType GetCsr() {
    return make_tuple(ptr_, idx_, val_);
  }


 private:
  vector<int> ptr_;
  vector<int> idx_;
  vector<T> val_;

  void SelfConstruct() { ptr_.resize(this->row() + 1); }
};
}  // namespace mat




// #endif

// #if defined(CSR0_MPI_ON)
// MpiVectorTool mpi_;
// #endif