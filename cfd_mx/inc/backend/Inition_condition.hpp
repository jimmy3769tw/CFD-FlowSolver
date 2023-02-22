#pragma once
#include "controlPanel.hpp"


void Locality(
    Simulation& simu,
    CalDomain& Lo,
    velocity& T0,
    velocity& T1,
    velocity& T3,
    pressure& t1,
    grid& gridA
)
{
   const auto [nx, ny, nz, gC] = gridA.nxyzgC; 

   for (size_t i = Lo.x_start  ; i < Lo.x_end ; ++i )
   for (size_t j = Lo.y_start  ; j < Lo.y_end ; ++j )
   for (size_t k = Lo.z_start  ; k < Lo.z_end ; ++k )
   {
      const int icel = i*nz*ny + j*nz + k;
      T0.u[icel] = 0.;         T0.v[icel] = 0.;         T0.w[icel] = 0.;
      T1.u[icel] = 0.;         T1.v[icel] = 0.;         T1.w[icel] = 0.;
      T3.u[icel] = 0.;         T3.v[icel] = 0.;         T3.w[icel] = 0.;
      t1.p[icel] = 0.;

      if (i == 2) {
         T0.u[icel - nz*ny] = 0.;
         T0.u[icel - 2*nz*ny] = 0.;
      }
      if (i == nx -3 ) {
         T0.u[icel + nz*ny] = 0.;
         T0.u[icel + 2*nz*ny] = 0.;
      }
      if (i == 2) {
         T0.u[icel-nz] = 0.;
         T0.u[icel-2*nz] = 0.;
      }
      if (i == nx -3 ) {
         T0.u[icel + nz] = 0.;
         T0.u[icel + 2*nz] = 0.;
      }
      if (k == 2) {
         T0.u[icel - 1] = 0.;
         T0.u[icel - 2] = 0.;
      }
      if (k == nz -3 ) {
         T0.u[ icel + 1] = 0.;
         T0.u[ icel + 2] = 0.;
      }
   }
}

