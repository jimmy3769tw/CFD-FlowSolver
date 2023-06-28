#pragma once
#include <vector>
#include <cmath>

#include "../math.hpp"
#include "../../mat_mpi_tool.hpp"

namespace solver{
    using namespace std;
    template<typename matrixT>
    class BicgstabRestartMpi
    {
        public:

        BicgstabRestartMpi(){}

        BicgstabRestartMpi(matrixT & lhs_mat){ init(lhs_mat);}

        virtual ~BicgstabRestartMpi(){ }

        // ! init
        inline void init(matrixT lhs_mat)
        {
            mpi_ = lhs_mat.mpi_;
            lhs_mat_ = lhs_mat;
            length_ = lhs_mat.row();
            r0_.resize(length_, 0.0);
            r1_.resize(length_, 0.0);
            p1_.resize(length_, 0.0);
            s1_.resize(length_, 0.0);
            ap_.resize(length_, 0.0);
            as_.resize(length_, 0.0);
            px1_.resize(length_, 0.0);
            px2_.resize(length_, 0.0);
        }

        void SetTolerance(double tolerance) { tol_ = tolerance;}

        std::tuple<int, double> solve(const vector<double>& rhs,
                                      vector<double>& x);

        std::tuple<int, double> operator()(const std::vector<double>& rhs,
                                           std::vector<double>& x) {
            return solve(rhs, x);
        }

        private:
            int length_;
            matrixT lhs_mat_;
            mat::MpiVectorTool mpi_;

            int iters_max_{5000};

            // ----------------------
            double tol_ = 1e-3,
                   pRestartFactor_   = 1e-7,
                   restartTol_ = pRestartFactor_ * tol_;
            // ----------------------

            vector<double> r0_;
            vector<double> r1_;
            vector<double> p1_;
            vector<double> s1_;
            vector<double> ap_;
            vector<double> as_;

            // ----------------------------
            vector<double> px1_;
            vector<double> px2_;
            // ----------------------------

                
            // ----------------------------
            int timestep_ = 0;
            int totalNumIter_ = 0, 
                accumulativeResetMembers_ = 0,
                pResets_ = 0,
                numRestarts_  = 0;
            // ----------------------------

    };



    // ! main
    template <typename matrixT>
    inline std::tuple<int, double> BicgstabRestartMpi<matrixT>::solve(
        const vector<double>& rhs, vector<double>& x) {
            timestep_++;
            double r1r0 = 0, pre_r1r0 = 0, a = 0, w = 0, b = 0,
                   norm0 = mpi_.L2Norm(rhs);

            lhs_mat_.MultiplyMpi(x, r0_);

            // ---------------------------
            math::copy(px1_, px2_);
            math::copy(x, px1_);

#pragma omp parallel for default(none) shared(rhs, r0_)
            for (int i = mpi_.beg(); i < mpi_.end(); ++i) {
            r0_[i] = rhs[i] - r0_[i];
            }

            math::copy(r0_, r1_);
            math::copy(r0_, p1_);
            // ---------------------------

            r1r0 = mpi_.InnerProduct(r1_, r0_);

            int numIter_ = 0;

            while (true) {

            lhs_mat_.MultiplyMpi(p1_, ap_);

            a  = r1r0 / mpi_.InnerProduct(ap_, r0_);

            mpi_.Barrier();
#pragma omp parallel for schedule(static)
                for (int i = mpi_.beg(); i < mpi_.end(); ++i) { s1_[i] = r1_[i] - a*ap_[i]; }
                mpi_.Barrier();

                lhs_mat_.MultiplyMpi(s1_, as_);

                w = mpi_.InnerProduct(as_, s1_) / mpi_.InnerProduct(as_, as_);

                for (int i = mpi_.beg(); i < mpi_.end(); ++i) {
                  x[i] += a * p1_[i] + w * s1_[i];
                  r1_[i] = s1_[i] - w * as_[i];
            }

            // Check for terminating conditions
            if ((mpi_.L2Norm(r1_) / norm0) < tol_) { break; }
            // cout << mpi_.L2Norm(r1_) / norm0 << ", ";

            // ! ## Check for reset of bCGSTAB
            if (numIter_ >= iters_max_) {
                int ResetMembers = 0;

                for(int i = mpi_.beg(); i < mpi_.end(); ++i){
                    if ( std::isnan(x[i]) || std::isinf(x[i]) ) 
                    { x[i] = px1_[i], ++ResetMembers;}
                }

                ResetMembers = mpi_.AllreduceSum(ResetMembers);

                if (ResetMembers > length_*0.1) {

                    for(int i = mpi_.beg(); i<mpi_.end(); ++i)
                    { x[i] = px2_[i]; }
                    accumulativeResetMembers_ += length_, ++pResets_;
                } 
                else if (ResetMembers != 0)
                { accumulativeResetMembers_ += ResetMembers, ++pResets_; }

                break;
            }
            // *------------------------------------

            pre_r1r0 = r1r0;
            r1r0 = mpi_.InnerProduct(r1_, r0_);
            b  = (a / w) * (r1r0 / pre_r1r0);

            // Check rho for restart
            if (timestep_ < 10 || abs(pre_r1r0) > restartTol_)
            {

                #pragma omp parallel for
                for(int i = mpi_.beg(); i<mpi_.end(); ++i)
                { p1_[i] = r1_[i] + b*(p1_[i] - w*ap_[i]); }
            }
            else
            {
                for(int i = mpi_.beg(); i<mpi_.end(); ++i)
                {
                    r0_[i] = r1_[i];
                    p1_[i] = r1_[i];
                }
                ++numRestarts_;
            }
            ++numIter_;
            ++totalNumIter_;
        }
        mpi_.Allocate(x);
        return std::make_tuple(numIter_, mpi_.L2Norm(r1_) / norm0);
    }

} // solver