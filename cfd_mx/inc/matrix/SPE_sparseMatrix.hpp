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

namespace mat
{
    using namespace std; 
    template<typename T>
    class SPE_matrix
    {

        public:
        using CSR_type = 
                // <ptr, indices, values>
                std::tuple< std::vector<int>,
                std::vector<int>,
                std::vector<T> >;
                // -----------------------

        //* ------ Constructor & Destructor ---------
        SPE_matrix(){}

        SPE_matrix(int nx, int ny, int nz){construct(nx, ny, nz);}

        void resize(int nx, int ny, int nz){construct(nx, ny, nz);}
        
        virtual ~SPE_matrix(){ destruct(); };
    
        void setupPressure(
            int gC, 
            const std::vector<double> &Dx, 
            const std::vector<double> &Dy, 
            const std::vector<double> &Dz,
            const std::vector<double> &Dxs, 
            const std::vector<double> &Dys, 
            const std::vector<double> &Dzs 
        ){

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

                        if (i!=0)     { val_[ic7+ 0] = -xm(i); }else{ val_[ic7+ 0] = double(); }
                        if (j!=0)     { val_[ic7+ 1] = -ym(j); }else{ val_[ic7+ 1] = double(); }
                        if (k!=0)     { val_[ic7+ 2] = -zm(k); }else{ val_[ic7+ 2] = double(); }
                                        val_[ic7+ 3] =  ijk;
                        if (k!=nz_-1) { val_[ic7+ 4] = -zp(k); }else{ val_[ic7+ 4] = double(); }
                        if (j!=ny_-1) { val_[ic7+ 5] = -yp(j); }else{ val_[ic7+ 5] = double(); }
                        if (i!=nx_-1) { val_[ic7+ 6] = -xp(i); }else{ val_[ic7+ 6] = double(); }
                    }
                }
            }
        }






        std::vector<double> set_Jp_pc(){
            std::vector<double> jp;
            #pragma omp parallel for 
            for (int i=0; i<nx_; ++i)
            {
                for (int j=0; j<ny_; ++j)
                {
                    for (int k=0; k<nz_; ++k)
                    {
                        int ic = k + j*nz_ + i*nynz_;
                        jp[ic] = 1.0/ val_[ic*7+ 3] ;
                        if (i!=0)     { val_[ic*7+ 0] *= jp[ic]; }
                        if (j!=0)     { val_[ic*7+ 1] *= jp[ic]; }
                        if (k!=0)     { val_[ic*7+ 2] *= jp[ic]; }
                                        val_[ic*7+ 3] *= jp[ic];
                        if (k!=nz_-1) { val_[ic*7+ 4] *= jp[ic]; }
                        if (j!=ny_-1) { val_[ic*7+ 5] *= jp[ic]; }
                        if (i!=nx_-1) { val_[ic*7+ 6] *= jp[ic]; }
                    }
                }
            }

            // --------------
            val_jp_.resize(ndim_*6);
            jp_on_ = true;
            // --------------
            for (int i=0; i<nx_; ++i)
            {
                for (int j=0; j<ny_; ++j)
                {
                    for (int k=0; k<nz_; ++k)
                    {
                        int ic = k + j*nz_ + i*nynz_;
                        jp[ic] = 1.0/ val_[ic*7+ 3] ;
                        if (i!=0)     { val_jp_[ic*6 + 0] = val_[ic*7+ 0] *= jp[ic]; }
                        if (j!=0)     { val_jp_[ic*6 + 1] = val_[ic*7+ 1] *= jp[ic]; }
                        if (k!=0)     { val_jp_[ic*6 + 2] = val_[ic*7+ 2] *= jp[ic]; }
                        if (k!=nz_-1) { val_jp_[ic*6 + 3] *= jp[ic]; }
                        if (j!=ny_-1) { val_jp_[ic*6 + 4] *= jp[ic]; }
                        if (i!=nx_-1) { val_jp_[ic*6 + 5] *= jp[ic]; }
                    }
                }
            }

            return jp;
        }


        //  ! ----------------------------------- Spmv -----------------------------------

        void spmvPointerAB(std::vector<double> &x , std::vector<double> &r){
            auto VAL = val_.data();
            auto R = r.data();
            auto X = x.data();
            #pragma omp parallel for simd default(none) \
            firstprivate(R, X, VAL)
            for (int i = 0; i < ndim_ ; ++i){
                double sumB {0.0};
                double sumA {0.0};
                auto i7 = i * 7;
                if (i-nynz_ >= 0)       { sumA +=  *(VAL+i7+0) *  *(X+i-nynz_); }
                if (i-nz_ >= 0)         { sumB +=   *(VAL+i7+1) *  *(X+i-nz_  ); }
                if (i-1 >= 0)           { sumA +=  *(VAL+i7+2) *  *(X+i-1    ); }
                                          sumB  +=   *(VAL+i7+3) *  *(X+i      );
                if (i+1 <  ndim_)       { sumA +=  *(VAL+i7+4) *  *(X+i+1    ); }
                if (i+nz_ <  ndim_)     { sumB  +=  *(VAL+i7+5) *  *(X+i+nz_  ); }
                if (i+nynz_ <  ndim_)   { sumA +=  *(VAL+i7+6) *  *(X+i+nynz_); }
                *(R+i) = sumB + sumA;
            }
        }



        inline void 
        spmvJp_omp(std::vector<double> &x , std::vector<double> &r){



            for (int i = 0; i < ndim_ ; ++i){
                double sum = 0.0;
                auto i6 = i * 6;
                if (i-nynz_ >= 0)       { sum +=  val_jp_[i6+ 0] * x[i-nynz_]; }
                if (i-nz_ >= 0)         { sum +=  val_jp_[i6+ 1] * x[i-nz_  ]; }
                if (i-1 >= 0)           { sum +=  val_jp_[i6+ 2] * x[i-1    ]; }
                                          sum +=  x[i];
                if (i+1 <  ndim_)       { sum +=  val_jp_[i6+ 3] * x[i+1    ]; }
                if (i+nz_ <  ndim_)     { sum +=  val_jp_[i6+ 4] * x[i+nz_  ]; }
                if (i+nynz_ <  ndim_)   { sum +=  val_jp_[i6+ 5] * x[i+nynz_]; }
                r[i] = sum;
            }
        }




        inline void 
        multiply_omp (std::vector<double> &x , std::vector<double> &r){
            if (jp_on_) {spmvJp_omp(x, r);}
            else {spmv_omp(x, r);}
        }




        inline void spmv_omp(std::vector<double> &x , std::vector<double> &r){
        
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


        inline void 
        multiply(std::vector<double> &x , std::vector<double> &r){

            for (int i = 0; i < ndim_ ; ++i){
                double sum = 0.0;
                auto i7 = i * 7;
                if (i-nynz_ >= 0)       { sum +=  val_[i7+ 0] * x[i-nynz_]; }
                if (i-nz_ >= 0)         { sum +=  val_[i7+ 1] * x[i-nz_  ]; }
                if (i-1 >= 0)           { sum +=  val_[i7+ 2] * x[i-1    ]; }
                                          sum +=  val_[i7+ 3] * x[i      ];
                if (i+1 <  ndim_)       { sum +=  val_[i7+ 4] * x[i+1    ]; }
                if (i+nz_ <  ndim_)     { sum +=  val_[i7+ 5] * x[i+nz_  ]; }
                if (i+nynz_ <  ndim_)   { sum +=  val_[i7+ 6] * x[i+nynz_]; }
                r[i] = sum;
            }
        }

        inline void 
        multiplyJp(std::vector<double> &x , std::vector<double> &r){

            for (int i = 0; i < ndim_ ; ++i){
                double sum = 0.0;
                auto i6 = i * 6;
                if (i-nynz_ >= 0)       { sum +=  val_jp_[i6+ 0] * x[i-nynz_]; }
                if (i-nz_ >= 0)         { sum +=  val_jp_[i6+ 1] * x[i-nz_  ]; }
                if (i-1 >= 0)           { sum +=  val_jp_[i6+ 2] * x[i-1    ]; }
                                          sum +=  x[i];
                if (i+1 <  ndim_)       { sum +=  val_jp_[i6+ 3] * x[i+1    ]; }
                if (i+nz_ <  ndim_)     { sum +=  val_jp_[i6+ 4] * x[i+nz_  ]; }
                if (i+nynz_ <  ndim_)   { sum +=  val_jp_[i6+ 5] * x[i+nynz_]; }
                r[i] = sum;
            }
        }


        #if defined (SPE_MPI_ON)


        void mpi_init(MPI_Comm comm_world, int total)
            { mpi_.init(comm_world, total); }

        auto GetMpi()const{return mpi_;}


        inline bool
        spmvJp_mpi( std::vector<T> & x, std::vector<T> & r) 
        {
            mpi_.Allocate(x);
            for(auto i = mpi_.beg(); i < mpi_.end();++i){
                double sum = 0.0;
                auto i6 = i * 6;
                if (i-nynz_ >= 0)       { sum +=  val_jp_[i6+ 0] * x[i-nynz_]; }
                if (i-nz_ >= 0)         { sum +=  val_jp_[i6+ 1] * x[i-nz_  ]; }
                if (i-1 >= 0)           { sum +=  val_jp_[i6+ 2] * x[i-1    ]; }
                                          sum +=  x[i];
                if (i+1 <  ndim_)       { sum +=  val_jp_[i6+ 3] * x[i+1    ]; }
                if (i+nz_ <  ndim_)     { sum +=  val_jp_[i6+ 4] * x[i+nz_  ]; }
                if (i+nynz_ <  ndim_)   { sum +=  val_jp_[i6+ 5] * x[i+nynz_]; }
                r[i] = sum;
            }
            return true;
        }



        inline bool
        spmv_mpi( std::vector<T> & x, std::vector<T> & r) 
        {
            mpi_.Allocate(x);
            for(auto i = mpi_.beg(); i < mpi_.end();++i){
                double sum = 0.0;
                auto i7 = i * 7;
                if (i-nynz_ >= 0)       { sum +=  val_[i7+ 0] * x[i-nynz_]; }
                if (i-nz_ >= 0)         { sum +=  val_[i7+ 1] * x[i-nz_  ]; }
                if (i-1 >= 0)           { sum +=  val_[i7+ 2] * x[i-1    ]; }
                                          sum +=  val_[i7+ 3] * x[i      ];
                if (i+1 <  ndim_)       { sum +=  val_[i7+ 4] * x[i+1    ]; }
                if (i+nz_ <  ndim_)     { sum +=  val_[i7+ 5] * x[i+nz_  ]; }
                if (i+nynz_ <  ndim_)   { sum +=  val_[i7+ 6] * x[i+nynz_]; }
                r[i] = sum;
            }
            return true;
        }


        inline bool 
        multiply_mpi( std::vector<T> & x, std::vector<T> & r) 
        {
            if (jp_on_) {spmvJp_mpi(x, r);}
            else {spmv_mpi(x, r);}
        }

        #endif
      


        int row() const {return ndim_;}

        private:



            //---------------------------------------------
            #if defined (SPE_MPI_ON)
                MpiVectorTool mpi_;
            #endif
            //---------------------------------------------

            bool jp_on_ = false;
            int ndim_, nx_, ny_, nz_, nynz_; 
            std::vector<double> val_;
            std::vector<double> val_jp_;

            void construct(int nx, int ny, int nz){ 

                ndim_ = nx*ny*nz;
                nx_ = nx;
                ny_ = ny;
                nz_ = nz;
                nynz_ = ny*nz;
                val_.resize(ndim_*7);
            }

            void destruct(void){}
    };

}



