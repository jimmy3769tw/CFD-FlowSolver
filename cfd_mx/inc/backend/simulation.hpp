#pragma once
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "time_variance_authority.hpp"


struct IniCondition{
  double u = double();
  double v = double();
  double w = double();
};

class Simulation {
 public:
  std::vector<std::string> io_string;
  int tid = 0;
  int pid = 0;
  TimeVarianceAuthority tva;
  Simulation& SetReynoldsNumber(double re) {
    reynolds_no_ = re;
    nu_ = 1.0 / re;
    isSetRe_ = true;
    return *this;
  }
  auto GetReynoldsNumber() { return reynolds_no_; }
  auto GetNu() { return nu_; }

  Simulation& SetBoundaryCondition(std::string str) {
    boundaryConditionName_ = str;
    isSetBoundaryCondition_ = false;
    return *this;
  }
  std::string GetBoundaryCondition() { return boundaryConditionName_; }

  size_t poisson_eq_iter_max_ = 5000;
  double poisson_criteria = 1.e-3;
  // amgcl def 1e-8
  // ahmad def 1e-3

  void SetArgcArgv(int argc, char** argv) {
    argc_ = argc;
    argv_ = argv;
  }

  auto GetArgc() { return argc_; }

  auto GetArgv() { return argv_; }


  IniCondition ini_condition;

  Simulation& SetVelocityInitialCondition(double u, double v, double w) {
    ini_condition.u = u;
    ini_condition.v = v;
    ini_condition.w = w;
    return *this;
  }


 private:
  bool isSetRe_ = false;
  bool isSetBoundaryCondition_ = false;
  std::string boundaryConditionName_;
  double dtMin = 1e-12,
         /* Timstep adjustment: CFL       */ cflFactor = 0.05,
         /* Timstep adjustment: Grid Fo   */ gridFoFactor = 0.15;
  double reynolds_no_;
  double nu_;
  double terminal_time_;
  int argc_;
  char **argv_;
};

