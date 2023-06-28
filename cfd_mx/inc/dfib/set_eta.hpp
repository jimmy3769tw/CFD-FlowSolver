// #pragma once
// #include <cmath>

// #include "backend/physical_variables.hpp" // immersedBoundary
// #include "grid/structured_grid.hpp"
// #include "backend/domain.hpp"
// #include "set_eta_cylinder.hpp"

// class CylinderZ;

// class CalEtaInterface {
//  public:
//   virtual void solver(immersedBoundary& dfib, CalDomain& domain) = 0;
// };
// class CalEtaFactory {
//   public:
//   CalEtaFactory(){}


//   CalEtaFactory(StructuredGrid& grid) { 
//     Init(grid); 
//   }

//   ~CalEtaFactory() {
//     delete cal_eta_imp_;
//   }

//   CalEtaFactory& Init(StructuredGrid& grid) { 
//     grid_ = &grid; 
//     return *this; 
//   }

//   CylinderZ& SetCylinderZ() {
//     cal_eta_imp_ = new CylinderZ(*grid_);
//     return *(CylinderZ*)cal_eta_imp_;
//   }

// private:
//   CalEtaInterface* cal_eta_imp_ = nullptr;
//   StructuredGrid* grid_;
// };
