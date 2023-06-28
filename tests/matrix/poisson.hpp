
#pragma once

#include <vector>
class HeatConduction{
 public:
  HeatConduction(int n) : nxy(n) {
    int n2 = nx * ny;  // Number of points in the grid.
    x.resize(n2, 0.0);

    ptr.clear();
    ptr.reserve(n2 + 1);
    ptr.push_back(0);
    idx.clear();
    idx.reserve(n2 * 5);  // We use 5-point stencil, so the matrix
    val.clear();
    val.reserve(n2 * 5);  // will have at most n2 * 5 nonzero elements.

    rhs.resize(n2);

    for (int j = 0, k = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i, ++k) {
        if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1) {
          // Boundary point. Use Dirichlet condition.
          idx.push_back(k);
          val.push_back(1.0);

          rhs[k] = 0.0;
        } else {
          // Interior point. Use 5-point finite difference stencil.
          idx.push_back(k - ny);
          val.push_back(-1.0 / (dx * dx));

          idx.push_back(k - 1);
          val.push_back(-1.0 / (dx * dx));

          idx.push_back(k);
          val.push_back(2.0 / (dx * dx) + 2.0 / (dy * dy));

          idx.push_back(k + 1);
          val.push_back(-1.0 / (dy * dy));

          idx.push_back(k + ny);
          val.push_back(-1.0 / (dy * dy));

          rhs[k] = 1.0;
        }

        ptr.push_back(idx.size());
      }
    }
  }

  auto GetCsr() { return std::make_tuple(ptr, idx, val); }
  std::vector<int> ptr;
  std::vector<int> idx;
  std::vector<double> val;
  std::vector<double> rhs;
  std::vector<double> x;

 private:
  int nxy;
  int nx = nxy, ny = nxy;

  const int gc = 1;
  const double ly = 1;
  const double lx = 1;

  int ndim = nx * ny;

  double dx = lx/nx;
  double dy = ly/ny;



  int icel(int i, int j) {
    return i*ny+j;
  }

  int icelShift(int i, int j) {
    return (i-gc) + (j-gc)*nx;
  }
};
