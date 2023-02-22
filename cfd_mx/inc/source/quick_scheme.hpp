#pragma once
#include <tuple>
#include <utility>  // std::pair, std::make_pair
#include <vector>

// ! ======================== Quick ========================
std::pair<double, double> CalQuickSpanwise(
    const double &Up, const double &Un, const int &pp, const int &p,
    const int &ic, const int &n, const int &nn, const std::vector<double> &Di,
    const std::vector<double> &Ds, const int &Idx, std::vector<double> &phi) {
  double phi_P, phi_N;

  if (Up > 0.0)
    phi_P =
        0.5 * (phi[ic] + phi[p]) -
        0.125 * Di[Idx + 1] * Di[Idx + 1] / Ds[Idx] *
            ((phi[p] - phi[ic]) / Di[Idx + 1] - (phi[ic] - phi[n]) / Di[Idx]);
  else
    phi_P = 0.5 * (phi[ic] + phi[p]) - 0.125 * Di[Idx + 1] * Di[Idx + 1] /
                                           Ds[Idx + 1] *
                                           ((phi[pp] - phi[p]) / Di[Idx + 2] -
                                            (phi[p] - phi[ic]) / Di[Idx + 1]);

  if (Un > 0.0)
    phi_N = 0.5 * (phi[n] + phi[ic]) - 0.125 * Di[Idx] * Di[Idx] / Ds[Idx - 1] *
                                           ((phi[ic] - phi[n]) / Di[Idx] -
                                            (phi[n] - phi[nn]) / Di[Idx - 1]);
  else
    phi_N = 0.5 * (phi[n] + phi[ic]) - 0.125 * Di[Idx] * Di[Idx] / Ds[Idx] *
                                           ((phi[p] - phi[ic]) / Di[Idx + 1] -
                                            (phi[ic] - phi[n]) / Di[Idx]);

  return make_pair(phi_P, phi_N);
}

inline std::pair<double, double> CalQuickChordwise(
    const double &Up, const double &Un, const int &pp, const int &p,
    const int &ic, const int &n, const int &nn, const std::vector<double> &Di,
    const std::vector<double> &Ds, const int &Idx, std::vector<double> &phi) {
  double phi_P, phi_N;
  // ----------------------------------
  if (Up > 0.0)
    phi_P = 0.5 * (phi[ic] + phi[p]) - 0.125 * Ds[Idx] * Ds[Idx] / Di[Idx] *
                                           ((phi[p] - phi[ic]) / Ds[Idx] -
                                            (phi[ic] - phi[n]) / Ds[Idx - 1]);
  else
    phi_P = 0.5 * (phi[ic] + phi[p]) - 0.125 * Ds[Idx] * Ds[Idx] / Di[Idx + 1] *
                                           ((phi[pp] - phi[p]) / Ds[Idx + 1] -
                                            (phi[p] - phi[ic]) / Ds[Idx]);
  // ----------------------------------

  // ----------------------------------
  if (Un > 0.0)
    phi_N = 0.5 * (phi[n] + phi[ic]) - 0.125 * Ds[Idx - 1] * Ds[Idx - 1] /
                                           Di[Idx - 1] *
                                           ((phi[ic] - phi[n]) / Ds[Idx - 1] -
                                            (phi[n] - phi[nn]) / Ds[Idx - 2]);
  else
    phi_N =
        0.5 * (phi[n] + phi[ic]) -
        0.125 * Ds[Idx - 1] * Ds[Idx - 1] / Di[Idx] *
            ((phi[p] - phi[ic]) / Ds[Idx] - (phi[ic] - phi[n]) / Ds[Idx - 1]);
  // ----------------------------------
  return make_pair(phi_P, phi_N);
}

// * ======================== Quick ========================
