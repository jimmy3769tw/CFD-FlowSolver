#pragma once
#include <vector>
#include <cmath>
#include "math.hpp"
#include "linear_solver.hpp"

namespace solver {
using namespace std;
using namespace math;
template <typename matrixT>
class BicgstabRestart : public LinearSolver<matrixT> {
 public:
  BicgstabRestart() {}

  BicgstabRestart(matrixT& lhs_mat) : LinearSolver<matrixT>(lhs_mat) {
    init(lhs_mat);
  }

  virtual ~BicgstabRestart() {}

  inline void init(matrixT& lhs_mat) {
    length_ = this->lhs_mat_.row();
    r0_.resize(length_, 0.0);
    r1_.resize(length_, 0.0);
    p1_.resize(length_, 0.0);
    s1_.resize(length_, 0.0);
    ap_.resize(length_, 0.0);
    as_.resize(length_, 0.0);
    px1_.resize(length_, 0.0);
    px2_.resize(length_, 0.0);
  }

  void SetTolerance(double tolerance) override {
    this->tolerance_ = tolerance;
    restartTol_ = pRestartFactor_ * this->tolerance_;
  }

  inline std::tuple<int, double> solve(const vector<double>& rhs,
                                       vector<double>& x) override {
    timestep_++;
    double r1r0 = 0, pre_r1r0 = 0, a = 0, w = 0, b = 0,
           norm0 = math::L2Norm(rhs);

    this->lhs_mat_.Multiply(x, r0_);
    // ---------------------------
    math::copy(px1_, px2_);
    math::copy(x, px1_);

#pragma omp parallel for schedule(static) default(none) shared(rhs, r0_)
    for (int i = 0; i < length_; ++i) {
      r0_[i] = rhs[i] - r0_[i];
    }

    math::copy(r0_, r1_);
    math::copy(r0_, p1_);
    // ---------------------------

    r1r0 = math::InnerProduct(r1_, r0_);

    int no_iters_ = 0;

    while (true) {
      this->lhs_mat_.Multiply(p1_, ap_);
      a = r1r0 / math::InnerProduct(ap_, r0_);
#pragma omp parallel for schedule(static)
      for (int i = 0; i < length_; ++i) {
        s1_[i] = r1_[i] - a * ap_[i];
      }

      this->lhs_mat_.Multiply(s1_, as_);

      w = math::InnerProduct(as_, s1_) / math::InnerProduct(as_, as_);

#pragma omp parallel for schedule(static)
      for (int i = 0; i < length_; ++i) {
        x[i] += a * p1_[i] + w * s1_[i];
        r1_[i] = s1_[i] - w * as_[i];
      }

      // Check for terminating conditions
      if ((math::L2Norm(r1_) / norm0) < this->tolerance_) {
        break;
      }

      // ! ## Check for reset of bCGSTAB
      if (no_iters_ >= this->iters_max_) {
        int ResetMembers = 0;

#pragma omp parallel for schedule(static) reduction(+ : ResetMembers)
        for (int i = 0; i < length_; ++i) {
          if (std::isnan(x[i]) || std::isinf(x[i])) {
            x[i] = px1_[i], ++ResetMembers;
          }
        }

        if (ResetMembers > length_ * 0.1) {
          math::copy(px2_, x);

          accumulativeResetMembers_ += length_, ++pResets_;

        } else if (ResetMembers != 0) {
          accumulativeResetMembers_ += ResetMembers, ++pResets_;
        }

        break;
      }
      // *------------------------------------

      pre_r1r0 = r1r0;
      r1r0 = math::InnerProduct(r1_, r0_);
      b = (a / w) * (r1r0 / pre_r1r0);

      // Check rho for restart
      if (timestep_ < 10 || std::abs(pre_r1r0) > restartTol_) {
#pragma omp parallel for schedule(static)
        for (int i = 0; i < length_; ++i) {
          p1_[i] = r1_[i] + b * (p1_[i] - w * ap_[i]);
        }
      } else {
        math::copy(r1_, r0_);
        math::copy(r1_, p1_);
        ++numRestarts_;
      }
      ++no_iters_;
      ++total_no_iters_;
    }

    return std::make_tuple(no_iters_, math::L2Norm(r1_) / norm0);
  }

 private:
  int length_;
  double pRestartFactor_ = 1e-7,
         restartTol_ = pRestartFactor_ * this->tolerance_;

  vector<double> r0_;
  vector<double> r1_;
  vector<double> p1_;
  vector<double> s1_;
  vector<double> ap_;
  vector<double> as_;
  vector<double> px1_;
  vector<double> px2_;

  size_t timestep_ = 0;
  size_t total_no_iters_ = 0, accumulativeResetMembers_ = 0, pResets_ = 0,
         numRestarts_ = 0;
};

}  // namespace solver
