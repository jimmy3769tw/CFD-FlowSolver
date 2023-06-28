#pragma once

#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "linear_solver_mpi.hpp"
#include "../../mat_mpi_tool.hpp"
#include "../math.hpp"

namespace solver{
    template <typename matrixT>
    class BicgstabMpi : public LinearSolverMpi<matrixT> {
     public:
      BicgstabMpi() {}

      BicgstabMpi(matrixT& lhs_mat) : LinearSolverMpi<matrixT>(lhs_mat) {
        init(lhs_mat);
      }

      virtual ~BicgstabMpi() {}

      inline void init(const matrixT& lhs_mat) {
        this->Construct(lhs_mat);
        length_ = this->lhs_mat_.row();
        p_.resize(length_);
        r_.resize(length_);
        r2_.resize(length_);
        v_.resize(length_);
        ss_.resize(length_);
        t_.resize(length_);
      }

      inline std::tuple<int, double> solve(const std::vector<double>& rhs,
                                           std::vector<double>& x) override {
        int iters{0};
        double norm{0};

        double norm_rhs = this->mpi_.L2Norm(rhs);

        this->lhs_mat_.MultiplyMpi(x, p_);

        for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i)
          r_[i] = rhs[i] - p_[i];

        for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i)
          r2_[i] = r_[i];

        rho1_ = 1;
        alpha_ = 1;
        omega_ = 1;

        for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i) {
          v_[i] = 0.0;
        }

        for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i) {
          p_[i] = 0.0;
        }

        norm = this->mpi_.L2Norm(r_) / norm_rhs;

        iters = 0;
        while (norm > this->tolerance_ && iters < this->iters_max_) {
          ++iters;

          rho2_ = this->InnerProduct(r2_, r_);

          beta_ = (rho2_ / rho1_) * (alpha_ / omega_);

          for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i)
            p_[i] = r_[i] + beta_ * (p_[i] - omega_ * v_[i]);

          this->lhs_mat_.MultiplyMpi(p_, v_);

          alpha_ = rho2_ / this->InnerProduct(r2_, v_);

          for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i)
            ss_[i] = r_[i] - alpha_ * v_[i];

          this->lhs_mat_.MultiplyMpi(ss_, t_);

          omega_ = this->InnerProduct(t_, ss_) / this->InnerProduct(t_, t_);

          for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i)
            x[i] += alpha_ * p_[i] + omega_ * ss_[i];

          for (int i = this->mpi_.beg(); i < this->mpi_.end(); ++i)
            r_[i] = ss_[i] - omega_ * t_[i];

          rho1_ = rho2_;

          norm = 0;

          norm = this->InnerProduct(r_, r_);
          norm = sqrt(norm) / norm_rhs;
        }

        this->mpi_.Allocate(x);

        return std::make_tuple(iters, norm);
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
}