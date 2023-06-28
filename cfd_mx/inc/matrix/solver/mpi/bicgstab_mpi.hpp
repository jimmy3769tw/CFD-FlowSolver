#pragma once

#include <numeric>
#include <tuple>
#include <vector>

#include "../../mat_mpi_tool.hpp"
#include "../math.hpp"

namespace solver{

    using namespace std;
    using namespace math;
    using namespace mat;

    template<typename matrixT>
    class BicgstabMpi
    {
        public:

        BicgstabMpi(){}

        BicgstabMpi(matrixT & lhs_mat){ init(lhs_mat);}

        virtual ~BicgstabMpi(){ }

        // ! init ----------------------
        inline bool init(const matrixT & lhs_mat)
        {
            lhs_mat_ = lhs_mat;
            mpi_ = lhs_mat.mpi_;
            length_ = lhs_mat.row();
            p_.resize(length_);
            r_.resize(length_);
            r2_.resize(length_);
            v_.resize(length_);
            ss_.resize(length_);
            t_.resize(length_);

            return true;
        }
        // * --------------------------
        void setTolerance(double tolerance) { zeta_ = tolerance;}

        std::pair<int, double> solve(const std::vector<double> & rhs, std::vector<double> & x);

        std::pair<int, double> operator()(const std::vector<double>& rhs,
                                          std::vector<double>& x) {
            return solve(rhs, x);
        }

        private:
            // --------------------
            matrixT lhs_mat_;
            MpiVectorTool mpi_;
            int length_;
            // --------------------

            // --------------------
            double zeta_{1e-5};
            int iters_max_{3000};
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


        double InnerProduct(const std::vector<double> & rhs,const  std::vector<double> & x)
        {
            return mpi_.InnerProduct(rhs,x);
        }
    };


    // ! main
    template<typename matrixT>
    inline std::pair<int, double> BicgstabMpi<matrixT>::
    solve( const std::vector<double> & rhs, std::vector<double> & x)
    {

        int iters{0};
        double norm{0};

        // ---------------------
        double norm_rhs = mpi_.L2Norm(rhs);
        // ---------------------

        lhs_mat_.MultiplyMpi(x, p_);

        for(int i=mpi_.beg(); i<mpi_.end(); ++i)
            r_[i] = rhs[i] - p_[i];

        for(int i=mpi_.beg(); i<mpi_.end(); ++i)
            r2_[i] = r_[i];

        rho1_  = 1; alpha_ = 1; omega_ = 1;

        for(int i=mpi_.beg(); i<mpi_.end(); ++i)
        { v_[i] = 0.0;}

        for(int i=mpi_.beg(); i<mpi_.end(); ++i)
        { p_[i] = 0.0;}

        norm = mpi_.L2Norm(r_) / norm_rhs;

        iters = 0;
        while(norm>zeta_ && iters<iters_max_)
        {
            ++iters;

            rho2_ = mpi_.InnerProduct(r2_, r_);

            beta_ = (rho2_/rho1_) * (alpha_/omega_);

            for(int i=mpi_.beg(); i<mpi_.end(); ++i)
                p_[i] = r_[i] + beta_ * (p_[i] - omega_ * v_[i]);

            lhs_mat_.MultiplyMpi(p_, v_);

            alpha_ = rho2_ / mpi_.InnerProduct(r2_, v_);

            for(int i=mpi_.beg(); i<mpi_.end(); ++i)
                ss_[i] = r_[i] - alpha_ * v_[i];

            lhs_mat_.MultiplyMpi(ss_, t_);

            omega_ = mpi_.InnerProduct(t_, ss_) / mpi_.InnerProduct(t_, t_);

            for(int i=mpi_.beg(); i<mpi_.end(); ++i)
                x[i] += alpha_ * p_[i] + omega_ * ss_[i];

            for(int i=mpi_.beg(); i<mpi_.end(); ++i)
                r_[i] = ss_[i] - omega_ * t_[i];

            rho1_ = rho2_;

            norm = 0;

            norm = mpi_.InnerProduct(r_, r_);
            norm = sqrt(norm) / norm_rhs;
        }

        mpi_.Allocate(x);

        return make_pair(iters, norm);
    }
}