#pragma once

#include <tuple>
#include <vector>
#include <numeric>
#include "math.hpp"

using std::cout;

namespace solver{
    using namespace std;
    using namespace math;
    template<typename matrixT>
    class Bicgstab
    {
        public:

        Bicgstab(){ }

        Bicgstab(matrixT & lhs_mat){ init(lhs_mat);}

        virtual ~Bicgstab(){ }

        // ! init ----------------------
        inline bool init(const matrixT & lhs_mat)
        {
            lhs_mat_ = lhs_mat;
            length_ = lhs_mat.row();
            p_.resize(length_);
            r_.resize(length_);
            r2_.resize(length_);
            v_.resize(length_);
            ss_.resize(length_);
            t_.resize(length_);

            return true;
        }

        void SetTolerance(double tolerance) { zeta_ = tolerance;}

        std::pair<int, double> solve(const std::vector<double> & rhs, std::vector<double> & x); 

        std::pair<int, double> operator ()(const std::vector<double> & rhs, std::vector<double> & x)
        { return solve(rhs, x); }

        private:
            // --------------------
            matrixT lhs_mat_;
            int length_;
            // --------------------

            // --------------------
            double zeta_{1e-5};
            int iters_max_{120};
            // --------------------

            // --------------------
            double  alpha_, beta_, omega_, rho1_, rho2_, rMath_;
            // --------------------

            // --------------------
            std::vector<double> p_;
            std::vector<double> r_;
            std::vector<double> r2_;
            std::vector<double> v_;
            std::vector<double> ss_;
            std::vector<double> t_;
            // --------------------

    };


    // ! main
    template<typename matrixT>
    inline std::pair<int, double> Bicgstab<matrixT>::
    solve( 
        const std::vector<double> & rhs, std::vector<double> & x)
    {
        int iters{0};
        double norm{0};

        // ---------------------
        double norm_rhs = L2Norm(rhs);
        // ---------------------

        lhs_mat_.multiply_omp(x, p_);

        for(int i=0; i<length_; ++i)
            r_[i] = rhs[i] - p_[i];

        for(int i=0; i<length_; ++i) 
            r2_[i] = r_[i];

        rho1_  = 1; alpha_ = 1; omega_ = 1;

        for(int i=0; i<length_; ++i)
        { v_[i] = 0.0;}

        for(int i=0; i<length_; ++i)
        { p_[i] = 0.0;}

        norm = L2Norm(r_) / norm_rhs;

        iters = 0;
        while(norm>zeta_ && iters<iters_max_)
        {
            ++ iters;

            rho2_ = InnerProduct(r2_, r_);

            beta_ = (rho2_/rho1_) * (alpha_/omega_);

            for(int i=0; i < length_; ++i)
                p_[i] = r_[i] + beta_ * (p_[i] - omega_ * v_[i]);

            lhs_mat_.multiply_omp(p_, v_);

            alpha_ = rho2_ / InnerProduct(r2_, v_);

            for(int i=0; i<length_; ++i)
                ss_[i] = r_[i] - alpha_ * v_[i];

            lhs_mat_.multiply_omp(ss_, t_);

            omega_ = InnerProduct(t_, ss_) / InnerProduct(t_, t_);

            for(int i=0; i<length_; ++i)
                x[i] += alpha_ * p_[i] + omega_ * ss_[i];

            for(int i=0; i<length_; ++i)
                r_[i] = ss_[i] - omega_ * t_[i];

            rho1_ = rho2_;

            norm = 0;

            norm = InnerProduct(r_, r_);
            norm = sqrt(norm) / norm_rhs;
        }

        return make_pair(iters, norm);
    }
}