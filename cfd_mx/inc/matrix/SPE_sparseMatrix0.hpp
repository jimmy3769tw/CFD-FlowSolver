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



#ifdef SPE0_MPI_ON
    #include "MAT_mpi.hpp"
#endif 

namespace mat
{

    template<typename T>
    class SPE_matrix0
    {

        public:

        using CSR_type = 
                // <ptr, indices, values>
                std::tuple< std::vector<int>,
                std::vector<int>,
                std::vector<T> >;
                // -----------------------

        //* ------ Constructor & Destructor ---------
        SPE_matrix0(){}

        SPE_matrix0(int nx, int ny, int nz){construct(nx, ny, nz);}

        void resize(int nx, int ny, int nz){construct(nx, ny, nz);}
        
        virtual ~SPE_matrix0(){ destruct(); };
    
        void setupPressure(
            int gC, 
            const std::vector<double> &Dx, 
            const std::vector<double> &Dy, 
            const std::vector<double> &Dz,
            const std::vector<double> &Dxs, 
            const std::vector<double> &Dys, 
            const std::vector<double> &Dzs 
        ){


            #pragma omp parallel for 
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
                        auto ijk =  xm(i) + ym(j) + zm(k) + 
                                    xp(i) + yp(j) + zp(k) ;
                        
                        if (i==0)     { ijk -= xm(i); }
                        if (i==nx_-1) { ijk -= xp(i); }
                        if (j==0)     { ijk -= ym(j); }
                        if (j==ny_-1) { ijk -= yp(j); }
                        if (k==0)     { ijk -= zm(k); }
                        if (k==nz_-1) { ijk -= zp(k); }
                        if (i!=0)     { val_[ic+0*ndim_] = -xm(i); }else{val_[ic+0*ndim_] = double();}
                        if (j!=0)     { val_[ic+1*ndim_] = -ym(j); }else{val_[ic+1*ndim_] = double();}
                        if (k!=0)     { val_[ic+2*ndim_] = -zm(k); }else{val_[ic+2*ndim_] = double();}
                                        val_[ic+3*ndim_] =  ijk;
                        if (k!=nz_-1) { val_[ic+4*ndim_] = -zp(k); }else{val_[ic+4*ndim_] = double();}
                        if (j!=ny_-1) { val_[ic+5*ndim_] = -yp(j); }else{val_[ic+5*ndim_] = double();}
                        if (i!=nx_-1) { val_[ic+6*ndim_] = -xp(i); }else{val_[ic+6*ndim_] = double();}
                    }
                }
            }
        }



        std::vector<double> set_Jp_pc(){

            std::vector<double> jp(ndim_);

            #pragma omp parallel for 
            for (int i=0; i<nx_; ++i)
            {
                for (int j=0; j<ny_; ++j)
                {
                    for (int k=0; k<nz_; ++k)
                    {
                        int ic = k + j*nz_ + i*nynz_;
                        jp[ic] = 1. / val_[ic+3*ndim_] ;
                        if (i!=0)     { val_[ic+0*ndim_] *= jp[ic]; }
                        if (j!=0)     { val_[ic+1*ndim_] *= jp[ic]; }
                        if (k!=0)     { val_[ic+2*ndim_] *= jp[ic]; }
                                        val_[ic+3*ndim_] *= jp[ic];
                        if (k!=nz_-1) { val_[ic+4*ndim_] *= jp[ic]; }
                        if (j!=ny_-1) { val_[ic+5*ndim_] *= jp[ic]; }
                        if (i!=nx_-1) { val_[ic+6*ndim_] *= jp[ic]; }
                    }
                }
            }
            return jp;
        }


        void multiply_omp(std::vector<double> &x , std::vector<double> &r){

            #pragma omp parallel for schedule(static) \
            firstprivate(ndimV_, ndim_)
            for (int ic = 0; ic < ndim_ ; ++ic){
                double sum = 0.0;
                if (ic-nynz_ >= 0)       { sum +=  val_[ic+0*ndim_] * x[ic-nynz_]; }
                if (ic-nz_ >= 0)         { sum +=  val_[ic+1*ndim_] * x[ic-nz_  ]; }
                if (ic-1 >= 0)           { sum +=  val_[ic+2*ndim_] * x[ic-1    ]; }
                                           sum +=  val_[ic+3*ndim_] * x[ic      ];
                if (ic+1 <  ndim_)       { sum +=  val_[ic+4*ndim_] * x[ic+1    ]; }
                if (ic+nz_ <  ndim_)     { sum +=  val_[ic+5*ndim_] * x[ic+nz_  ]; }
                if (ic+nynz_ <  ndim_)   { sum +=  val_[ic+6*ndim_] * x[ic+nynz_]; }
                r[ic] = sum;
            }

        }


        #if defined (SPE0_MPI_ON)

        void multiply_mpi(std::vector<double> &x , std::vector<double> &r){

            mpi_.Allocate(x);

            #pragma omp parallel for schedule(static) \
            firstprivate(ndimV_, ndim_)
            for (int ic = mpi_.beg(); ic < mpi_.end() ; ++ic){
                double sum = 0.0;
                if (ic-nynz_ >= 0)       { sum +=  val_[ic+0*ndim_] * x[ic-nynz_]; }
                if (ic-nz_ >= 0)         { sum +=  val_[ic+1*ndim_] * x[ic-nz_  ]; }
                if (ic-1 >= 0)           { sum +=  val_[ic+2*ndim_] * x[ic-1    ]; }
                                           sum +=  val_[ic+3*ndim_] * x[ic      ];
                if (ic+1 <  ndim_)       { sum +=  val_[ic+4*ndim_] * x[ic+1    ]; }
                if (ic+nz_ <  ndim_)     { sum +=  val_[ic+5*ndim_] * x[ic+nz_  ]; }
                if (ic+nynz_ <  ndim_)   { sum +=  val_[ic+6*ndim_] * x[ic+nynz_]; }
                r[ic] = sum;
            }

        }
        #endif 




        void multiply(std::vector<double> &x , std::vector<double> &r){

            for (int ic = 0; ic < ndim_ ; ++ic){
                double sum = 0.0;
                if (ic-nynz_ >= 0)       { sum +=  val_[ic+0*ndim_] * x[ic-nynz_]; }
                if (ic-nz_ >= 0)         { sum +=  val_[ic+1*ndim_] * x[ic-nz_  ]; }
                if (ic-1 >= 0)           { sum +=  val_[ic+2*ndim_] * x[ic-1    ]; }
                                           sum +=  val_[ic+3*ndim_] * x[ic      ];
                if (ic+1 <  ndim_)       { sum +=  val_[ic+4*ndim_] * x[ic+1    ]; }
                if (ic+nz_ <  ndim_)     { sum +=  val_[ic+5*ndim_] * x[ic+nz_  ]; }
                if (ic+nynz_ <  ndim_)   { sum +=  val_[ic+6*ndim_] * x[ic+nynz_]; }
                r[ic] = sum;
            }
        }

        int row() const {return ndim_;}



        #if defined (SPE0_MPI_ON)
            void mpi_init(MPI_Comm comm_world, int total)
                { mpi_.init(comm_world, total); }

            void mpi_init(MPI_Comm comm_world, int total, std::vector<int> beg, std::vector<int> end)
                { mpi_.init(comm_world, total, beg, end); }

            auto GetMpi()const{return mpi_;}
        #endif // SPE0_MPI_ON


        private:

            int ndim_, nx_, ny_, nz_, nynz_;
            int ndimV_[7] ;



            std::vector<double> val_;

            void construct(int nx, int ny, int nz){ 

                ndim_ = nx*ny*nz;
                nx_ = nx;
                ny_ = ny;
                nz_ = nz;
                nynz_ = ny*nz;
                val_.resize(ndim_*7);

                for (int i = 0; i < 7; i++)
                {
                    ndimV_[i] = i*ndim_;
                }
            }

            void destruct(void){}

         
            //---------------------------------------------
            #if defined (SPE0_MPI_ON)
                mat::MpiVectorTool mpi_;
            #endif
            //---------------------------------------------


    };

}



