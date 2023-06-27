#pragma once 

#include <omp.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <utility> // (since C++11) std::pair
#include <algorithm> // (until C++11)
#include <tuple>
// ! row major
#ifdef SPE_MPI_ON
#include "MAT_mpi.hpp"
#endif

#include 

namespace mat
{
    template<typename T>
    class StructGridMat
    {
     public:
      StructGridMat() {}

      StructGridMat(int cal_nx, int cal_ny, int cal_nz) { construct(cal_nx, cal_ny, cal_nz); }

      void resize(int cal_nx, int cal_ny, int cal_nz) { construct(cal_nx, cal_ny, cal_nz); }
      int row() const { return ndim_; }
      int col() const { return ndim_; }

      void Init(int gC, const std::vector<double> &Dx,
                const std::vector<double> &Dy, const std::vector<double> &Dz,
                const std::vector<double> &Dxs, const std::vector<double> &Dys,
                const std::vector<double> &Dzs) {
#pragma omp parallel for schedule(static)
            for (int i=0; i<nx_; ++i)
            {
                for (int j=0; j<ny_; ++j)
                {
                    for (int k=0; k<nz_; ++k)
                    {
                        auto xp = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC);};
                        auto yp = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC);};
                        auto zp = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC);};

                        auto xm = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC-1);};
                        auto ym = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC-1);};
                        auto zm = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC-1);};
                        int ic = k + j*nz_ + i*nynz_;
                        auto ic7 = ic*7;
                        auto ijk =  xm(i) + ym(j) + zm(k) + 
                                    xp(i) + yp(j) + zp(k) ;
                        
                        if (i==0)     { ijk -= xm(i); }
                        if (i==nx_-1) { ijk -= xp(i); }
                        if (j==0)     { ijk -= ym(j); }
                        if (j==ny_-1) { ijk -= yp(j); }
                        if (k==0)     { ijk -= zm(k); }
                        if (k==nz_-1) { ijk -= zp(k); }

                        if (i!=0)     { val_[ic7+ 0] = -xm(i); }else{ val_[ic7+ 0] = T(); }
                        if (j!=0)     { val_[ic7+ 1] = -ym(j); }else{ val_[ic7+ 1] = T(); }
                        if (k!=0)     { val_[ic7+ 2] = -zm(k); }else{ val_[ic7+ 2] = T(); }
                                        val_[ic7+ 3] =  ijk;
                        if (k!=nz_-1) { val_[ic7+ 4] = -zp(k); }else{ val_[ic7+ 4] = T(); }
                        if (j!=ny_-1) { val_[ic7+ 5] = -yp(j); }else{ val_[ic7+ 5] = T(); }
                        if (i!=nx_-1) { val_[ic7+ 6] = -xp(i); }else{ val_[ic7+ 6] = T(); }
                    }
                }
            }
      }

        //  ! ----------------------------------- Spmv -----------------------------------

        inline void Multiply(std::vector<double> &x, std::vector<double> &r) {
#pragma omp parallel for schedule(static)
            for (int i = 0; i < ndim_ ; ++i){
                double sum {0.0};
                auto i7 = i * 7;
                if (i-nynz_ >= 0)       { sum +=  val_[i7+0] * x[i-nynz_]; }
                if (i-nz_ >= 0)         { sum +=  val_[i7+1] * x[i-nz_  ]; }
                if (i-1 >= 0)           { sum +=  val_[i7+2] * x[i-1    ]; }
                                          sum +=  val_[i7+3] * x[i      ];
                if (i+1 <  ndim_)       { sum +=  val_[i7+4] * x[i+1    ]; }
                if (i+nz_ <  ndim_)     { sum +=  val_[i7+5] * x[i+nz_  ]; }
                if (i+nynz_ <  ndim_)   { sum +=  val_[i7+6] * x[i+nynz_]; }
                r[i] = sum;
            }
        }

        void mpi_init(MPI_Comm comm_world, int total)
            { mpi_.init(comm_world, total); }

        auto GetMpi()const{return mpi_;}

        inline void MultiplyMpi(std::vector<T> &x, std::vector<T> &r) {
            mpi_.Allocate(x);
            for (auto i = mpi_.beg(); i < mpi_.end(); ++i) {
                double sum = 0.0;
                auto i7 = i * 7;
                if (i - nynz_ >= 0) {
                    sum += val_[i7 + 0] * x[i - nynz_];
                }
                if (i - nz_ >= 0) {
                    sum += val_[i7 + 1] * x[i - nz_];
                }
                if (i - 1 >= 0) {
                    sum += val_[i7 + 2] * x[i - 1];
                }
                sum += val_[i7 + 3] * x[i];
                if (i + 1 < ndim_) {
                    sum += val_[i7 + 4] * x[i + 1];
                }
                if (i + nz_ < ndim_) {
                    sum += val_[i7 + 5] * x[i + nz_];
                }
                if (i + nynz_ < ndim_) {
                    sum += val_[i7 + 6] * x[i + nynz_];
                }
                r[i] = sum;
            }
        }

        private:
         int ndim_, nx_, ny_, nz_, nynz_;
         std::vector<double> val_;

         void construct(int cal_nx, int cal_ny, int cal_nz) {
            ndim_ = cal_nx * cal_ny * cal_nz;
            nx_ = cal_nx;
            ny_ = cal_ny;
            nz_ = cal_nz;
            nynz_ = cal_ny * cal_nz;
            val_.resize(ndim_ * 7);
         }
    };

}



