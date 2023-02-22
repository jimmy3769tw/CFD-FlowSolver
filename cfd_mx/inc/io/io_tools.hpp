#pragma once
#include <dirent.h>
#include <sys/types.h>
#include "../import/stl.hpp"

void CreatOutputFile(){
  if (opendir("Information") == NULL) {
    if (system("mkdir Information") != 0) {
      throw std::invalid_argument("Failed: Creat a output file");
    }
  }

  if (opendir("mx_out") == NULL) {
    if (system("mkdir mx_out") != 0) {
      throw std::invalid_argument("Failed: Creat a output file");
    }
  }

  if (opendir("Information/Chronograph") == NULL) {
    if (system("mkdir Information/Chronograph") != 0) {
      throw std::invalid_argument("Failed: Creat a output file");
    }
  }
}