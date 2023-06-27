#pragma once
#include <cmath>

class CalEtaInterface {
 public:
  virtual void solver(phi::ImmersedBoundary& dfib, CalDomain& domain) = 0;
};


class CalEtaFactory {
  public:
  CalEtaFactory(){};

  CalEtaFactory(StructuredGrid& grid) { Init(grid); }
  ~CalEtaFactory() {
    delete cal_eta_imp_;
  }

  CalEtaFactory& Init(StructuredGrid& grid) { grid_ = &grid; return *this; }

  CylinderZ& SetCylinderZ() {
    cal_eta_imp_ = new CylinderZ(*grid_);
    return *(CylinderZ*)cal_eta_imp_;
  }
  
  CalEtaInterface* cal_eta_imp_ = nullptr;

  StructuredGrid* grid_;
};
