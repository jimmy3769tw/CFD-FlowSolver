#pragma once
#include "grid/structured_grid.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

class PhysicalVal{
  public:

  PhysicalVal(StructuredGrid &grid): grid_(&grid){}
 protected:
  StructuredGrid* grid_;


  template <typename T>
  inline size_t icelCal(T i, T j, T k) const;

  template <typename T>
  inline size_t icel(T i, T j, T k) const;


  inline void Fill(std::vector<double> &v, double val) {
    #pragma omp for schedule(static)
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] = val;
    }
  }
};

template <typename T>
size_t PhysicalVal::icel(T i, T j, T k) const {
    return grid_->icel(i, j, k);
}


template <typename T>
size_t PhysicalVal::icelCal(T i, T j, T k) const {
    return grid_->icelCal(i, j, k);
}


// undefined reference to `unsigned long PhysicalVal::icel<int>(int, int, int) const'

class StaggeredVelocity : public PhysicalVal{
  std::vector<double> viseff_;

 public:
  StaggeredVelocity(StructuredGrid &grid);

  double &viseff(int i, int j, int k){
    return viseff_[icelCal(i, j, k)];
  }

  std::vector<double> u, v, w;

  template <typename T>
  double &U(T i, T j, T k) {
    return u[icel(i, j, k)];
  }

  template <typename T>
  auto &V(T i, T j, T k) {
    return v[icel(i, j, k)];
  }

  template <typename T>
  auto &W(T i, T j, T k) {
    return w[icel(i, j, k)];
  }

  template <typename T>
  auto U(T i, T j, T k) const {
    return u[icel(i, j, k)];
  }

  template <typename T>
  auto V(T i, T j, T k) const {
    return v[icel(i, j, k)];
  }

  template <typename T>
  double W(T i, T j, T k) const {
    return w[icel(i, j, k)];
  }

  template <typename T>
  std::vector<double> &U(T whichDim) {
    if (whichDim == 0) {
      return u;
    } else if (whichDim == 1) {
      return v;
    } else if (whichDim == 2) {
      return w;
    }
    throw std::runtime_error("out of dimensions!");
  }

  void FillVel(double ui, double vi, double wi) {
    Fill(u, ui);
    Fill(v, vi);
    Fill(w, wi);
  }

  auto GetGrid() {
    return grid_;
  }

 private:
  void resize(size_t n) {
    u.resize(n);
    v.resize(n);
    w.resize(n);
  }
};

class Pressure : public PhysicalVal{
 public:
  std::vector<double> p;
Pressure(StructuredGrid& grid);

  template <typename T>
  double& operator()(T i, T j, T k) {
    return p[icel(i, j, k)];
  }

  template <typename T>
  double operator()(T i, T j, T k) const {
    return p[icel(i, j, k)];
  }

  void InitPressure(double pi) {
    std::fill(p.begin(), p.end(), pi);
  }
};


class ImmersedBoundary : public PhysicalVal{
 private:
  size_t yShift_;
  size_t zShift_;

  public:
  ImmersedBoundary(StructuredGrid& grid);

  template <typename T>
  auto Fx(T i, T j, T k) const {
    return f[icelCal(i, j, k)];
  }

  template <typename T>
  auto Fy(T i, T j, T k) const {
    return f[icelCal(i, j, k) + yShift_];
  }

  template <typename T>
  auto Fz(T i, T j, T k) const {
    return f[icelCal(i, j, k) + zShift_];
  }

  template <typename T>
  auto Eta(T i, T j, T k) const {
    return eta[icel(i, j, k)];
  }

  template <typename T>
  auto &Fx(T i, T j, T k) {
    return f[icelCal(i, j, k)];
  }

  template <typename T>
  auto &Fy(T i, T j, T k) {
    return f[icelCal(i, j, k) + yShift_];
  }

  template <typename T>
  auto &Fz(T i, T j, T k) {
    return f[icelCal(i, j, k) + zShift_];
  }

  template <typename T>
  auto &Eta(T i, T j, T k) {
    return eta[icel(i, j, k)];
  }

  std::vector<double> f, eta;
  std::vector<double> Val_sumz;
  std::vector<double> Val_sumz_sumy;
  std::vector<double> ValSum;
  std::vector<double> cylinderCenter;
};

