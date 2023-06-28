#pragma once 

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <utility>  // (since C++11) std::pair
#include <algorithm> // (until C++11)
#include <tuple>
#include "sparse_mat.hpp"
namespace mat
{

template <typename T>
class EllMat : public SparseMat<T> {

   public:
    EllMat() : SparseMat<T>() {}

    EllMat(int rows, int cols, int nnz) : SparseMat<T>(rows, cols) {
      SelfConstruct(nnz);
    }
    EllMat(int n, int nnz) : SparseMat<T>(n, n) { SelfConstruct(nnz); }

    virtual ~EllMat() { Destruct(); };

    void resize(int rows, int cols, int nnz) {
      this->Construct(rows, cols);
      SelfConstruct(nnz);
    }

    void resize(int n, int nnz) {
      this->Construct(n, n);
      SelfConstruct(nnz);
     }

    void Multiply(std::vector<T> & x, std::vector<T> & r) override {
#pragma omp parallel for
      for (int i = 0; i < this->row(); ++i) {
        const int pt = i * no_no_zero_;
        r[i] = T();
        // for(int j = 0; j < curr_idx_size[i];++j){
        for (int j = 0; j < no_no_zero_; ++j) {
          r[i] += val_[pt + j] * x[idx_[pt + j]];
        }
      }
    }

    void MultiplyMpi(std::vector<T> & x, std::vector<T> & r) {
      this->mpi_.Allocate(x);
      for (auto i = this->mpi_.beg(); i < this->mpi_.end(); ++i) {
        const int pt = i * no_no_zero_;
        r[i] = T();

        for (int j = 0; j < no_no_zero_; ++j) {
          r[i] += val_[pt + j] * x[idx_[pt + j]];
        }
      }
    }

    int MaxRow() { return no_no_zero_; }
    void Show();
    bool analyse();
    EllMat<T>& Set(int row, int col, T val);
    EllMat<T>& Set(typename SparseMat<T>::CSR_type c);

    template <typename X>
    friend std::ostream& operator<<(std::ostream& os, const EllMat<X>& mat);

   private:
    int no_no_zero_, init_idx_;
    std::vector<int> idx_;
    std::vector<double> val_;
    std::vector<int> curr_idx_size;

    inline void SelfConstruct(int max_col_idx) {
      no_no_zero_ = max_col_idx;
      init_idx_ = this->col() - 1;
      idx_.resize(max_col_idx * this->row(), init_idx_);
      val_.resize(max_col_idx * this->row(), 0);
      curr_idx_size.resize(this->row(), 0);
    }
    void Destruct() {}
  };

  template <typename T>
  inline std::ostream& operator<<(std::ostream& os, const EllMat<T>& mat) {
    os << std::endl;
    for (int i = 0, iter; i < mat.row(); ++i) {
      for (int j = 0; j < mat.col(); ++j) {
        const int pt = i * mat.no_no_zero_;
        os << std::setw(4);
        for (iter = 0; iter < mat.curr_idx_size[i]; ++iter) {
          if (j == mat.idx_.at(pt + iter)) {
            os << mat.val_.at(pt + iter);
            break;
          }
        }
        if (iter == mat.curr_idx_size.at(i)) {
          os << "x";
        }
      }
      os << std::endl;
    }
    return os;
  }

  template <typename T>
  inline EllMat<T>& EllMat<T>::Set(typename SparseMat<T>::CSR_type c) {
    auto [prt, idx, val] = c;
    // ! We assume that your mat is a n by n mat for our convenient.
    this->Construct(prt.size() - 1, prt.size() - 1);
    no_no_zero_ = 0;

    // find max index
    for (int row = 0; row < this->row(); ++row) {
      no_no_zero_ = std::max(prt[row + 1] - prt[row], no_no_zero_);
    }

    SelfConstruct(no_no_zero_);
 


    for (int row = 0; row < this->row(); ++row) {
      for (int j = prt[row]; j < prt[row + 1]; ++j) {
        Set(row, idx[j], val[j]);
      }
    }
    return *this;
  }

  template <typename T>
  inline EllMat<T>& EllMat<T>::Set(int row, int col, T val) {
    if (row < 0 || col < 0) {
      throw std::invalid_argument("Matrix dimensions cannot be negative.");
    }

    if (col >= this->col()) {
      std::cout << col << " >= " << this->col() << std::endl;
      throw std::invalid_argument("colIdx >= no_no_zero_ ");
    }

    if (row >= this->row()) throw std::invalid_argument("row >= this->row()");

    const int pt = row * no_no_zero_;

    if (val == T()) {
      // bool is_new_one = true;
      for (int i = 0; i < curr_idx_size[row]; ++i) {
        if (col == idx_[pt + i]) {
          // is_new_one = false; --> remove
          for (int j = i; j < curr_idx_size[row]; ++j) {
            idx_[pt + j] = idx_[pt + j + 1];
            val_[pt + j] = idx_[pt + j + 1];
          }

          idx_[pt + curr_idx_size[row]] = init_idx_;
          val_[pt + curr_idx_size[row]] = T();
          curr_idx_size[row]--;

          return *this;
        }
      }
      return *this;
    } else {
      bool is_new_one = true;
      for (int i = 0; i < curr_idx_size[row]; ++i) {
        if (col == idx_[pt + i]) {
          is_new_one = false;
          idx_[pt + i] = col;
          val_[pt + i] = val;
          return *this;
        }
      }

      if (is_new_one) {
        idx_[pt + curr_idx_size[row]] = col;
        val_[pt + curr_idx_size[row]] = val;
        curr_idx_size[row]++;
      }

      if (curr_idx_size[row] > no_no_zero_) {
        throw std::invalid_argument(
            "[ELL0::Set] -> curr_idx_size[row] >no_no_zero_");
      }

      return *this;
    }
  }

  template <typename T>
  inline bool EllMat<T>::analyse(void) {
    // -----------------------------------
    std::vector<int> a(no_no_zero_ + 1, 0);
    for (auto v : curr_idx_size) {
      ++a[v];
    }
    // -----------------------------------
    int i{0};
    // -----------------------------------
    for (auto v : a) {
      std::cout << "i=" << i++ << ", " << v << "\n";
    }
    // -----------------------------------
    return true;
  }

  template <typename T>
  inline void EllMat<T>::Show() {
    std::cout << "   |val";
    for (int i = 0; i < no_no_zero_; ++i) {
      std::cout << std::setw(4) << "";
    }

    std::cout << "    |Idx\n";

    for (int i = 0; i < this->row(); ++i) {
      // ----------------------------------------

      std::cout << std::setw(3) << i << "|";

      // ----------------------------------------
      for (int j = 0; j < no_no_zero_; ++j) {
        std::cout << std::setw(4) << val_[i * no_no_zero_ + j] << " ";
      }
      // ----------------------------------------

      std::cout << "|";

      // ----------------------------------------
      for (int j = 0; j < no_no_zero_; ++j) {
        std::cout << std::setw(4) << idx_[i * no_no_zero_ + j] << " ";
      }
      // ----------------------------------------

      std::cout << std::endl;
    }
  }



}
