#pragma once

vector<int> GetRandomVectorInt(int size, int range_begin = 0) {
  vector<int> result(size);
  for (auto& x : result) {
    x = rand() % size + range_begin;
  }
  return result;
}

vector<double> GetRandomVectorDouble(int size, int range_begin = 0) {
  vector<double> result(size);
  for (auto& x : result) {
    x = rand() % size + range_begin + rand() / size;
  }
  return result;
}