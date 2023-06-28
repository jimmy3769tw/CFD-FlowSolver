#pragma once
#include <cmath>

#include "backend/physical_variables.hpp"
#include "grid/structured_grid.hpp"
#include "backend/domain.hpp"

class ImmersedBoundary;
class CalEtaInterface {
 public:
  virtual void solver(ImmersedBoundary& dfib, LocalDomain& domain) = 0;
};
