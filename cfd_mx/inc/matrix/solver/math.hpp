
#pragma once 
#include <complex>
#include <vector>
#include <tuple>

#include<iostream>
namespace math{

    using namespace std;

    inline void init(vector<double> &a, const double val) {
      auto A = a.data();
      auto len = a.size();

#pragma omp parallel for simd default(none) firstprivate(A, val, len) \
    schedule(static)
      for (size_t i = 0; i < len; ++i) {
        *(A + i) = val;
      }
    }
    inline void init(vector<double> &a, const vector<double> &b) {
      size_t len = a.size();
      auto A = a.data();
      auto B = b.data();

#pragma omp parallel for simd default(none) firstprivate(A, len, B) \
    schedule(static)
      for (size_t i = 0; i < len; ++i) {
        *(A + i) = *(B + i);
      }
    }
    inline void copy(const vector<double> &a, vector<double> &b) { init(b, a); }



    inline double 
    InnerProduct(const vector<double> & a, const vector<double> & b)
    {
        double r{0.0};
        auto len = a.size();
        auto A = a.data();
        auto B = b.data();

        #pragma omp parallel for simd reduction(+:r) default(none) \
        firstprivate(A, B, len) schedule(static)
        for(size_t i=0; i<len ; ++i)
        {  r += *(A+i) * *(B+i); }

        return r;
    }


    inline double 
    L2Norm(const std::vector<double> & a)
    { return sqrt( InnerProduct(a,a) ); }


    inline void 
    zero(vector<double> & a){  init(a, 0.0); }
}