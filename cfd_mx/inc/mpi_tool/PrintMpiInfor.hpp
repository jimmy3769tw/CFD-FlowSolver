#pragma once
#include "mpi.h"
#include <iostream>
namespace mpi {
void PrintMpiInfo() {
  char name[1024];
  int length = 1024, minor;
  MPI_Get_processor_name(name, &length);

  int major;
  MPI_Get_version(&major, &minor);
  // MPI_Request req[10];

  std::cout << "\nMPI Version " << major << "." << minor << std::endl;
  std::cout << "This Project is from " << name << std::endl;
}
}  // namespace mpi