#pragma once
// #include <fstream>  // file output
#include <string>
#include <vector>

// void file(std::string fileN, std::vector<double> &result, int nx, int ny) {
//   auto dx = 1/nx;
//   auto dy = 1/ny;

//   std::vector<double> X(nx + 1, 0), Y(ny + 1, 0);

//   for (int i = 1; i < nx + 1; ++i) X[i] = X[i - 1] + dx;

//   for (int i = 1; i < ny + 1; ++i) Y[i] = Y[i - 1] + dy;

//   std::ofstream fileQ;

//   fileQ.precision(std::numeric_limits<double>::digits10 + 2);

//   fileQ.open(fileN);

//   fileQ << "TITLE= \"2D Navier-Stokes Solution\""
//         << "\n"
//         << "VARIABLES= x, y, z, T"
//         << "\n"
//         << "ZONE T= \"Single Zone\""
//         << ", I=" << std::to_string(nx) << ", J=" << std::to_string(ny) << ", K=1"
//         << ", F=POINT"
//         << "\n";

//   for (int j = 0; j < ny; ++j) {
//     for (int i = 0; i < nx; ++i) {
//       fileQ << X.at(i) << "\t" << Y.at(j) << "\t" << 0.0 << "\t"
//             << result.at(i * ny + j) << "\n";
//     }
//   }
// }
