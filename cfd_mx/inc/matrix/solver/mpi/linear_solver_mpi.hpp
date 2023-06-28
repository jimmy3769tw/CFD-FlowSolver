#pragma once

#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "../../mat_mpi_tool.hpp"
#include "../math.hpp"

namespace solver {

using namespace std;
using namespace math;
using namespace mat;

template <typename MatType>
class LinearSolverMpi {
 public:
  LinearSolverMpi() {}

  LinearSolverMpi(MatType& lhs_mat) { Construct(lhs_mat); }

  virtual ~LinearSolverMpi() {}

  void Construct(const MatType& lhs_mat) {
    // if (lhs_mat.mpi_.Empty()) {
    //   throw std::runtime_error("out of dimensions!");
    // }
    lhs_mat_ = lhs_mat;
    mpi_ = lhs_mat.mpi_;
  }

  void SetTolerance(double tolerance) { tolerance_ = tolerance; }

  void SetMaxIters(int niters) {

  }

// return iters, error
  std::tuple<int, double> virtual solve(const std::vector<double>& rhs,
                                        std::vector<double>& x) = 0;

  std::tuple<int, double> operator()(const std::vector<double>& rhs,
                                     std::vector<double>& x) {
    return solve(rhs, x);
  }

 protected:
  double tolerance_{1e-5};
  int iters_max_{3000};
  double InnerProduct(const std::vector<double>& rhs,
                      const std::vector<double>& x) {
    return mpi_.InnerProduct(rhs, x);
  }
  MatType lhs_mat_;
  MpiVectorTool mpi_;
};

}  // namespace solver