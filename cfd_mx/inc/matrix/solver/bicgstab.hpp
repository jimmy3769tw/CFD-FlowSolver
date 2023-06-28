#pragma once

#include <tuple>
#include <vector>
#include <numeric>
#include "math.hpp"
#include "linear_solver.hpp"

using std::cout;

namespace solver {
using namespace std;
using namespace math;

template <typename MatType>
class Bicgstab : public LinearSolver<MatType> {
 public:
  Bicgstab() {}

  Bicgstab(MatType& lhs_mat) : LinearSolver<MatType>(lhs_mat) { Init(lhs_mat); }

  virtual ~Bicgstab() {}

  inline void Init(const MatType& lhs_mat) {
    this->lhs_mat_ = lhs_mat;
    length_ = lhs_mat.row();
    p_.resize(length_);
    r_.resize(length_);
    r2_.resize(length_);
    v_.resize(length_);
    ss_.resize(length_);
    t_.resize(length_);
  }

  inline std::tuple<int, double> solve(
      const std::vector<double>& rhs, std::vector<double>& x) override {

    int iters{0};
    double norm{0};

    double norm_rhs = L2Norm(rhs);

    this->lhs_mat_.Multiply(x, p_);

    for (int i = 0; i < length_; ++i) r_[i] = rhs[i] - p_[i];

    for (int i = 0; i < length_; ++i) r2_[i] = r_[i];

    rho1_ = 1;
    alpha_ = 1;
    omega_ = 1;

    for (int i = 0; i < length_; ++i) {
      v_[i] = 0.0;
    }

    for (int i = 0; i < length_; ++i) {
      p_[i] = 0.0;
    }

    norm = L2Norm(r_) / norm_rhs;

    iters = 0;
    while (norm > this->tolerance_ && iters < this->iters_max_) {
      ++iters;

      rho2_ = InnerProduct(r2_, r_);

      beta_ = (rho2_ / rho1_) * (alpha_ / omega_);

      for (int i = 0; i < length_; ++i)
        p_[i] = r_[i] + beta_ * (p_[i] - omega_ * v_[i]);

      this->lhs_mat_.Multiply(p_, v_);

      alpha_ = rho2_ / InnerProduct(r2_, v_);

      for (int i = 0; i < length_; ++i) ss_[i] = r_[i] - alpha_ * v_[i];

      this->lhs_mat_.Multiply(ss_, t_);

      omega_ = InnerProduct(t_, ss_) / InnerProduct(t_, t_);

      for (int i = 0; i < length_; ++i)
        x[i] += alpha_ * p_[i] + omega_ * ss_[i];

      for (int i = 0; i < length_; ++i) r_[i] = ss_[i] - omega_ * t_[i];

      rho1_ = rho2_;

      norm = 0;

      norm = InnerProduct(r_, r_);
      norm = sqrt(norm) / norm_rhs;
    }

    return make_tuple(iters, norm);
  }

 private:
  int length_;
  double alpha_, beta_, omega_, rho1_, rho2_, rMath_;
  std::vector<double> p_;
  std::vector<double> r_;
  std::vector<double> r2_;
  std::vector<double> v_;
  std::vector<double> ss_;
  std::vector<double> t_;
};

}  // namespace solver