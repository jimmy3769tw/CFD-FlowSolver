#pragma once
#include <vector>
#include <cmath>
#include "math.hpp"

namespace solver{
    using namespace std;
    using namespace math;

    template<typename matrixT>
    class BicgstabRestart
    {
        public:

        BicgstabRestart(){}

        BicgstabRestart(matrixT & lhs_mat): lhs_mat_(lhs_mat){ init();}

        virtual ~BicgstabRestart(){ }

        inline void init() {
          length_ = lhs_mat_.row();
          r0_.resize(length_, 0.0);
          r1_.resize(length_, 0.0);
          p1_.resize(length_, 0.0);
          s1_.resize(length_, 0.0);
          ap_.resize(length_, 0.0);
          as_.resize(length_, 0.0);
          px1_.resize(length_, 0.0);
          px2_.resize(length_, 0.0);
        }

        void setTolerance(double tolerance) { 
            tol_ = tolerance;
            restartTol_ = pRestartFactor_ * tol_;
        }

        std::pair<int, double> solve(const vector<double> & rhs, vector<double> & x); 

        std::pair<int, double> operator()(const std::vector<double> & rhs, std::vector<double> & x)
        { return solve(rhs, x); }

        private:
            int length_;

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
    
            matrixT lhs_mat_;
    };

    // ! main
    template <typename matrixT>
    inline std::pair<int, double> BicgstabRestart<matrixT>::solve(
        const vector<double>& rhs, vector<double>& x) {

            timestep_++;
            double r1r0 = 0, pre_r1r0 = 0, a = 0, w = 0, b = 0,
                   norm0 = math::L2Norm(rhs);

            lhs_mat_.multiply_omp(x, r0_);
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

            int numIter_ = 0;

            while (true) {
            lhs_mat_.multiply_omp(p1_, ap_);
            a  = r1r0 / math::InnerProduct(ap_, r0_);
            #pragma omp parallel for schedule(static)
            for (int i=0; i<length_; ++i) 
            { s1_[i] = r1_[i] - a*ap_[i]; }

            lhs_mat_.multiply_omp(s1_, as_);

            w = math::InnerProduct(as_, s1_) / math::InnerProduct(as_, as_);

            #pragma omp parallel for schedule(static)
            for (int i=0; i<length_; ++i)
            {
                x[i] += a*p1_[i] + w*s1_[i];
                r1_[i] = s1_[i] - w*as_[i];
            }

            // Check for terminating conditions
            if ((math::L2Norm(r1_) / norm0) < tol_) { break; }

            // ! ## Check for reset of bCGSTAB
            if (numIter_ >= iters_max_) {
                int ResetMembers = 0;

                #pragma omp parallel for schedule(static) reduction(+:ResetMembers)
                for (int i=0; i<length_; ++i){
                    if ( std::isnan(x[i]) || std::isinf(x[i]) ) 
                    { x[i] = px1_[i], ++ResetMembers;}
                }

                if (ResetMembers > length_*0.1) {

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
            b  = (a / w) * (r1r0 / pre_r1r0);

            // Check rho for restart
            if (timestep_ < 10 || std::abs(pre_r1r0) > restartTol_)
            {
                #pragma omp parallel for schedule(static)
                for (int i=0; i<length_; ++i) 
                { p1_[i] = r1_[i] + b*(p1_[i] - w*ap_[i]); }
            }
            else
            {
                math::copy(r1_, r0_);
                math::copy(r1_, p1_);
                ++numRestarts_;
            }
            ++numIter_;
            ++totalNumIter_;
        }
    

        return std::make_pair( numIter_, math::L2Norm(r1_) / norm0);
    }
}
