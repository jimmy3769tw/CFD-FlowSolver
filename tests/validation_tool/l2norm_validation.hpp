#pragma once

template <typename T>
auto L2Norm(std::vector<T>& a, std::vector<T>& b) {
  T sum = 0;
  for (int i = 0; i < a.size(); i++) {
    sum += pow(a[i] - b[i], 2);
  }
  return sqrt(sum);
}
