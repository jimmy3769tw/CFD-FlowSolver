#pragma once

std::vector<int> GetRandomVectorInt(int size, int range_begin = 0) {
  std::vector<int> result(size);
  for (auto& x : result) {
    x = rand() % size + range_begin;
  }
  return result;
}

std::vector<double> GetRandomVectorDouble(int size, int range_begin = 0) {
  std::vector<double> result(size);
  for (auto& x : result) {
    x = rand() % size + range_begin + rand() / size;
  }
  return result;
}