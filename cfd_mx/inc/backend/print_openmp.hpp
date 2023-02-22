#pragma once
#include <unistd.h>  // for print pid

void PrintIsOpenmpExist() {
  printf("pid: %d\n", getpid());
#ifdef _OPENMP
  std::cout << "OpenMP: " << _OPENMP << std::endl;
#else
  std::cout << "Your compiler does not support OpenMP." << std::endl;
#endif
  printf("pid: %d\n", getpid());
}