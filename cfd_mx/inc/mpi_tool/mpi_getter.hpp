#pragma once
#include "mpi.h"
namespace mpi {

int GetRank(MPI_Comm word) {
  int rank;
  MPI_Comm_rank(word, &rank);
  return rank;
}

int GetRank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int GetSize() {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

int GetSize(MPI_Comm word) {
  int size;
  MPI_Comm_size(word, &size);
  return size;
}
}