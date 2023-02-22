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

#include <cstdlib>

#ifdef SPE1_MPI_ON
#include "MAT_mpi.hpp"
#endif


namespace mat
{
    using namespace std; 
    template<typename T>
    class SPE_matrix1
    {

        public:

        // //* ------ Constructor & Destructor ---------
        SPE_matrix1(){}

        SPE_matrix1(int nx, int ny, int nz){construct(nx, ny, nz);}

        void resize(int nx, int ny, int nz){construct(nx, ny, nz);}
        
        virtual ~SPE_matrix1(){ destruct(); };

        void setupPressure(
            int gC, 
            const std::vector<double> &Dx, 
            const std::vector<double> &Dy, 
            const std::vector<double> &Dz,
            const std::vector<double> &Dxs, 
            const std::vector<double> &Dys, 
            const std::vector<double> &Dzs 
        ){


            // ----------------
            xp_.resize(ndim_);
            yp_.resize(ndim_);
            zp_.resize(ndim_);
            // ----------------
            // ----------------
            xm_.resize(ndim_);
            ym_.resize(ndim_);
            zm_.resize(ndim_);
            // ----------------

            gC_ = gC;

            #pragma omp parallel for schedule(static)
            for (int i=0; i<nx_; ++i)
            for (int j=0; j<ny_; ++j)
            for (int k=0; k<nz_; ++k)
            {
                xp_[i]=  1.0 /Dx[i+gC] / Dxs[i+gC];
                xm_[i]=  1.0 /Dx[i+gC] / Dxs[i+gC-1];

                yp_[j]=  1.0 /Dy[j+gC] / Dys[j+gC];
                ym_[j]=  1.0 /Dy[j+gC] / Dys[j+gC-1];

                zp_[k]=  1.0 /Dz[k+gC] / Dzs[k+gC];
                zm_[k]=  1.0 /Dz[k+gC] / Dzs[k+gC-1];
            }


            #pragma omp parallel for schedule(static)
            for (int i=0; i<nx_; ++i)
            {
                for (int j=0; j<ny_; ++j)
                {
                    for (int k=0; k<nz_; ++k)
                    {
                        int ic = k + j*nz_ + i*nynz_;
                        auto ijk =  xm_[i] + ym_[j] + zm_[k] + 
                                    xp_[i] + yp_[j] + zp_[k] ;

                        if (i==0)     { ijk -= xm_[i]; }
                        if (i==nx_-1) { ijk -= xp_[i]; }
                        if (j==0)     { ijk -= ym_[j]; }
                        if (j==ny_-1) { ijk -= yp_[j]; }
                        if (k==0)     { ijk -= zm_[k]; }
                        if (k==nz_-1) { ijk -= zp_[k]; }

                        val_[ic] =  ijk;
                    }
                }
            }


        }


        std::vector<double> set_Jp_pc(){

            jp_.resize(ndim_);

            jp_on_ = true;

            #pragma omp parallel for schedule(static)
            for (int i=0; i<nx_; ++i)
            for (int j=0; j<ny_; ++j)
            for (int k=0; k<nz_; ++k)
            {
                int ic = k + j*nz_ + i*nynz_;
                jp_[ic] =  1.0 / val_[ic];
            }

            return jp_;
        }


        //  ! ----------------------------------- Spmv -----------------------------------

        inline void spmv_ompJP(std::vector<double> &x , std::vector<double> &r){

            #pragma omp parallel for schedule(static) 
            for (int i=0; i<nx_; ++i)
            for (int j=0; j<ny_; ++j)
            for (int k=0; k<nz_; ++k)
            {
                int ic = k + j*nz_ + i*nynz_;
                double sum {0.0};
                if (i!=0)     { sum += -xm_[i]*jp_[ic] * x[ic-nynz_]; }
                if (j!=0)     { sum += -ym_[j]*jp_[ic] * x[ic-nz_  ]; }
                if (k!=0)     { sum += -zm_[k]*jp_[ic] * x[ic-1    ]; }
                                sum +=  x[ic];
                if (k!=nz_-1) { sum += -zp_[k]*jp_[ic]  * x[ic+1    ]; }
                if (j!=ny_-1) { sum += -yp_[j]*jp_[ic]  * x[ic+nz_  ]; }
                if (i!=nx_-1) { sum += -xp_[i]*jp_[ic]  * x[ic+nynz_]; }
                r[ic] = sum;
            }
        }


        #if defined (SPE1_MPI_ON)

        inline void spmv_mpiJP(std::vector<double> &x , std::vector<double> &r){

             #pragma omp parallel for schedule(static) 
            for (auto ic = mpi_.beg() ; ic < mpi_.end() ; ic++){
                auto [i, rem] = div(ic, nynz_);
                auto [j, k] = div(rem, ny_);

                int ic = k + j*nz_ + i*nynz_;
                if( ic < mpi_.beg() || i >= mpi_.end()) continue;

                double sum {0.0};
                if (i!=0)     { sum += -xm_[i]*jp_[ic] * x[ic-nynz_]; }
                if (j!=0)     { sum += -ym_[j]*jp_[ic] * x[ic-nz_  ]; }
                if (k!=0)     { sum += -zm_[k]*jp_[ic] * x[ic-1    ]; }
                                sum +=  x[ic];
                if (k!=nz_-1) { sum += -zp_[k]*jp_[ic]  * x[ic+1    ]; }
                if (j!=ny_-1) { sum += -yp_[j]*jp_[ic]  * x[ic+nz_  ]; }
                if (i!=nx_-1) { sum += -xp_[i]*jp_[ic]  * x[ic+nynz_]; }
                r[ic] = sum;

            }
        }
        #endif // SPE1_MPI_ON


        inline void spmvJP(std::vector<double> &x , std::vector<double> &r){

            for (int i=0; i<nx_; ++i)
            for (int j=0; j<ny_; ++j)
            for (int k=0; k<nz_; ++k)
            {
                int ic = k + j*nz_ + i*nynz_;
                double sum {0.0};
                if (i!=0)     { sum += -xm_[i]*jp_[ic]  * x[ic-nynz_]; }
                if (j!=0)     { sum += -ym_[j]*jp_[ic]  * x[ic-nz_  ]; }
                if (k!=0)     { sum += -zm_[k]*jp_[ic]  * x[ic-1    ]; }
                                sum +=  x[ic];
                if (k!=nz_-1) { sum += -zp_[k]*jp_[ic]  * x[ic+1    ]; }
                if (j!=ny_-1) { sum += -yp_[j]*jp_[ic]  * x[ic+nz_  ]; }
                if (i!=nx_-1) { sum += -xp_[i]*jp_[ic]  * x[ic+nynz_]; }
                r[ic] = sum;
            }
        }
    

        inline void spmv_omp(std::vector<double> &x , std::vector<double> &r){

            #pragma omp parallel for schedule(static) 
            for (int i=0; i<nx_; ++i)
            for (int j=0; j<ny_; ++j)
            for (int k=0; k<nz_; ++k)
            {
                int ic = k + j*nz_ + i*nynz_;
                double sum {0.0};
                if (i!=0)     { sum += -xm_[i]  * x[ic-nynz_]; }
                if (j!=0)     { sum += -ym_[j]  * x[ic-nz_  ]; }
                if (k!=0)     { sum += -zm_[k]  * x[ic-1    ]; }
                                sum += val_[ic] * x[ic];
                if (k!=nz_-1) { sum += -zp_[k]  * x[ic+1    ]; }
                if (j!=ny_-1) { sum += -yp_[j]  * x[ic+nz_  ]; }
                if (i!=nx_-1) { sum += -xp_[i]  * x[ic+nynz_]; }
                r[ic] = sum;
            }
        }


        #if defined (SPE1_MPI_ON)

        inline void spmv_mpi(std::vector<double> &x , std::vector<double> &r){

            for (auto ic = mpi_.beg() ; ic < mpi_.end() ; ic++){
                auto [i, rem] = div(ic, nynz_);
                auto [j, k]   = div(rem, nz_);

                double sum {0.0};
                if (i!=0)     { sum += -xm_[i]  * x[ic-nynz_]; }
                if (j!=0)     { sum += -ym_[j]  * x[ic-nz_  ]; }
                if (k!=0)     { sum += -zm_[k]  * x[ic-1    ]; }
                                sum += val_[ic] * x[ic];
                if (k!=nz_-1) { sum += -zp_[k]  * x[ic+1    ]; }
                if (j!=ny_-1) { sum += -yp_[j]  * x[ic+nz_  ]; }
                if (i!=nx_-1) { sum += -xp_[i]  * x[ic+nynz_]; }
                r[ic] = sum;
            }
        }
        #endif // SPE1_MPI_ON



        inline void spmv(std::vector<double> &x , std::vector<double> &r){

            for (int i=0; i<nx_; ++i)
            for (int j=0; j<ny_; ++j)
            for (int k=0; k<nz_; ++k)
            {
                int ic = k + j*nz_ + i*nynz_;
                double sum {0.0};
                if (i!=0)     { sum += -xm_[i]  * x[ic-nynz_]; }
                if (j!=0)     { sum += -ym_[j]  * x[ic-nz_  ]; }
                if (k!=0)     { sum += -zm_[k]  * x[ic-1    ]; }
                                sum += val_[ic] * x[ic];
                if (k!=nz_-1) { sum += -zp_[k]  * x[ic+1    ]; }
                if (j!=ny_-1) { sum += -yp_[j]  * x[ic+nz_  ]; }
                if (i!=nx_-1) { sum += -xp_[i]  * x[ic+nynz_]; }
                r[ic] = sum;
            }
        }



        void multiply_omp (std::vector<double> &x , std::vector<double> &r){
            if (jp_on_) {spmv_ompJP(x, r);}
            else { spmv_omp(x, r); }
        }


        
        #if defined (SPE1_MPI_ON)
        void multiply_mpi (std::vector<double> &x , std::vector<double> &r){
            mpi_.Allocate(x);
            if (jp_on_) {spmv_mpiJP(x, r);}
            else { spmv_mpi(x, r); }
        }
        #endif  // SPE1_MPI_ON

        void multiply (std::vector<double> &x , std::vector<double> &r){
            if (jp_on_) {spmvJP(x, r);}
            else { spmv(x, r); }
        }


        //  * ----------------------------------- Spmv -----------------------------------


        #if defined (SPE1_MPI_ON)

        void mpi_init(MPI_Comm comm_world, int total)
            { mpi_.init(comm_world, total); }

        void mpi_init(MPI_Comm comm_world, int total, std::vector<int> beg, std::vector<int> end)
            { mpi_.init(comm_world, total, beg, end); }

        auto GetMpi()const{return mpi_;}

        #endif // SPE1_MPI_ON

        int row() const {return ndim_;}

        private:
            bool jp_on_ = false;

            int ndim_, nx_, ny_, nz_, nynz_, gC_; 
            std::vector<double> xp_, yp_, zp_, 
                                xm_, ym_, zm_;

            std::vector<double> jp_, val_;

            void construct(int nx, int ny, int nz){ 

                ndim_ = nx*ny*nz;
                nx_ = nx;
                ny_ = ny;
                nz_ = nz;
                nynz_ = ny*nz;
                val_.resize(ndim_);
            }




            //---------------------------------------------
            #if defined (SPE1_MPI_ON)
                mat::MpiVectorTool mpi_;
            #endif
            //---------------------------------------------

            void destruct(void){}
    };

}



