#pragma once

#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "math.hpp"
namespace solver {

using namespace std;
using namespace math;

template <typename matrixT>
class LinearSolver {
 public:
  LinearSolver() {}

  LinearSolver(matrixT& lhs_mat) { construct(lhs_mat); }

  virtual ~LinearSolver() {}

  void construct(const matrixT& lhs_mat) { lhs_mat_ = lhs_mat; }

  void virtual SetTolerance(double tolerance) { tolerance_ = tolerance; }

  void SetMaxIters(int niters) {}

  // return iters, error
  std::tuple<int, double> virtual solve(const std::vector<double>& rhs,
                                        std::vector<double>& x) = 0;

  std::tuple<int, double> operator()(const std::vector<double>& rhs,
                                    std::vector<double>& x) {
    return solve(rhs, x);
  }

 protected:
  matrixT lhs_mat_;
  double tolerance_{1e-5};
  int iters_max_{3000};
  double InnerProduct(const std::vector<double>& rhs,
                      const std::vector<double>& x) {
    return InnerProduct(rhs, x);
  }
};

}  // namespace solver