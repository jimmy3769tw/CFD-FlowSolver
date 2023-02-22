#pragma once
#include <utility>  // std::pair, std::make_pair

// ! ======================== Ud ========================
inline std::pair<double, double> CalUpwind(
    const double &Up, const double &Un, const int &p, const int &ic,
    const int &n, std::vector<double> &Di, std::vector<double> &Ds,
    const int &Idx, std::vector<double> &phi) {
  double phi_P, phi_N;

  if (Up > 0.0) {
    phi_P = phi[ic];
  } else {
    phi_P = phi[p];
  }

  if (Un > 0.0) {
    phi_N = phi[n];
  } else {
    phi_N = phi[ic];
  }

  return make_pair(phi_P, phi_N);
}
// * ======================== Ud ========================

// ! ======================== Lud ========================
inline std::pair<double, double> CalLudSpanwise(
    const double &Up, const double &Un, const int &pp, const int &p,
    const int &ic, const int &n, const int &nn, const std::vector<double> &Di,
    const std::vector<double> &Ds, const int &Idx, std::vector<double> &phi) {
  double phi_P, phi_N;

  if (Up > 0.0) {
    phi_P = (phi[ic] + 0.5 * Ds[Idx] * (phi[ic] - phi[n]) / Di[Idx]);
  } else {
    phi_P = (phi[p] + 0.5 * Ds[Idx + 1] * (phi[p] - phi[pp]) / Di[Idx + 1]);
  }

  if (Un > 0.0) {
    phi_N = (phi[n] + 0.5 * Ds[Idx - 1] * (phi[n] - phi[nn]) / Di[Idx - 1]);
  } else {
    phi_N = (phi[ic] + 0.5 * Ds[Idx] * (phi[ic] - phi[p]) / Di[Idx + 1]);
  }

  return make_pair(phi_P, phi_N);
}

inline std::pair<double, double> CalLudChordwise(
    const double &Up, const double &Un, const int &pp, const int &p,
    const int &ic, const int &n, const int &nn, const std::vector<double> &Di,
    const std::vector<double> &Ds, const int &Idx, std::vector<double> &phi) {
  double phi_P, phi_N;

  if (Up > 0.0) {
    phi_P = (phi[ic] + 0.5 * Ds[Idx] * (phi[ic] - phi[n]) / Di[Idx]);
  } else {
    phi_P = (phi[p] + 0.5 * Ds[Idx + 1] * (phi[p] - phi[pp]) / Di[Idx + 1]);
  }

  if (Un > 0.0) {
    phi_N = (phi[n] + 0.5 * Ds[Idx - 1] * (phi[n] - phi[nn]) / Di[Idx - 1]);
  } else {
    phi_N = (phi[ic] + 0.5 * Ds[Idx] * (phi[ic] - phi[p]) / Di[Idx + 1]);
  }

  return make_pair(phi_P, phi_N);
}
