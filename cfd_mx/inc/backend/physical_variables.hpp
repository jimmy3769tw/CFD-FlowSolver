#pragma once
#include "../grid/structured_grid.hpp"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
inline void Fill(std::vector<double> &v, double val) {
  #pragma omp for schedule(static)
  for (size_t i = 0; i < v.size(); ++i) {
    v[i] = val;
  }
}

class StaggeredVelocity {
  std::vector<double> viseff_;
 public:
  StaggeredVelocity(StructuredGrid &grid) : grid_(&grid) {
    resize(grid_->no_grid);
    viseff_.resize(grid_->cal_no_grid);
  }

  double &viseff(int i, int j, int k){
    return viseff_[grid_->icelCal(i, j, k)];
  }

  std::vector<double> u, v, w;
  template <typename T>
  auto &U(T i, T j, T k) {
    return u[grid_->icel(i, j, k)];
  }

  template <typename T>
  auto &V(T i, T j, T k) {
    return v[grid_->icel(i, j, k)];
  }

  template <typename T>
  auto &W(T i, T j, T k) {
    return w[grid_->icel(i, j, k)];
  }

  template <typename T>
  auto U(T i, T j, T k) const {
    return u[grid_->icel(i, j, k)];
  }

  template <typename T>
  auto V(T i, T j, T k) const {
    return v[grid_->icel(i, j, k)];
  }

  template <typename T>
  auto W(T i, T j, T k) const {
    return w[grid_->icel(i, j, k)];
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
  StructuredGrid *grid_;
  void resize(size_t n) {
    u.resize(n);
    v.resize(n);
    w.resize(n);
  }
};

class Pressure {
 public:
  std::vector<double> p;
Pressure(StructuredGrid& grid) : grid_(&grid){
      p.resize(grid.no_grid);
  }

  template <typename T>
  double& operator()(T i, T j, T k) {
    return p[grid_->icel(i, j, k)];
  }

  template <typename T>
  double operator()(T i, T j, T k) const {
    return p[grid_->icel(i, j, k)];
  }

  void InitPressure(double pi) {
    std::fill(p.begin(), p.end(), pi);
  }

  private:
  StructuredGrid* grid_;
};


class ImmersedBoundary {
  public:
  ImmersedBoundary(StructuredGrid& grid): grid_(&grid) {
    eta.resize(grid_->no_grid);
    f.resize(grid_->cal_no_grid * 3);
  }

  template <typename T>
  auto Fx(T i, T j, T k) const {
    return f[grid_->icelCal(i, j, k)];
  }

  template <typename T>
  auto Fy(T i, T j, T k) const {
    return f[grid_->icelCal(i, j, k) + yShift_];
  }

  template <typename T>
  auto Fz(T i, T j, T k) const {
    return f[grid_->icelCal(i, j, k) + zShift_];
  }

  template <typename T>
  auto Eta(T i, T j, T k) const {
    return eta[grid_->icel(i, j, k)];
  }

  template <typename T>
  auto &Fx(T i, T j, T k) {
    return f[grid_->icelCal(i, j, k)];
  }

  template <typename T>
  auto &Fy(T i, T j, T k) {
    return f[grid_->icelCal(i, j, k) + yShift_];
  }

  template <typename T>
  auto &Fz(T i, T j, T k) {
    return f[grid_->icelCal(i, j, k) + zShift_];
  }

  template <typename T>
  auto &Eta(T i, T j, T k) {
    return eta[grid_->icel(i, j, k)];
  }

  std::vector<double> f, eta;
  std::vector<double> Val_sumz;
  std::vector<double> Val_sumz_sumy;
  std::vector<double> ValSum;
  std::vector<double> cylinderCenter;

  double cylinderDimension;

 private:
  StructuredGrid* grid_;
  int yShift_ = grid_->cal_no_grid;
  int zShift_ = 2 * grid_->cal_no_grid;
};

