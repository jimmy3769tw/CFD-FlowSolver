#pragma once
#define CONVECTION_DIFUSSION_QUICK
// #define CONVECTION_DIFUSSION_LUD
// #define CONVECTION_DIFUSSION_UD

#include "grid/structured_grid.hpp"
#include "backend/physical_variables.hpp"
/*
* A present main , p is positve, n is negative, m is minus
* In the quick shcemem and non-nuiform case, there are some diferent between A and B.
*/
#include "quick_scheme.hpp"
#include "upwind_scheme.hpp"

void CalConvectionAndDiffusion(Simulation &simu, StaggeredVelocity &old_vel,
                               StaggeredVelocity &curr_vel,
                               const LocalDomain &location,
                               StructuredGrid &grid) {
  auto ii = [&](auto i, auto j, auto k) { return grid.icel(i, j, k); };
  const int nx = grid.nx, ny = grid.ny, nz = grid.nz;
  int d = 0, one = 1;
  auto gc = grid.no_ghost_cell;
  auto dt = simu.tva.GetDt();

  // ! ---------------------------------- x ----------------------------------

  if (location.x_end == nx - gc) {
    d = one;
  } else {
    d = 0;
  }

    #pragma omp parallel for
    for (auto i = location.x_start; i < location.x_end-d ; ++i )
    for (auto j = location.y_start; j < location.y_end   ; ++j )
    for (auto k = location.z_start; k < location.z_end   ; ++k )
    { //  =================== u =================== //

        int ic = grid.icel(i, j, k);
        // ! -------------- convection term --------------
        const auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] = grid.Get12Nb(ic);
        const auto interpolation = grid.GetInterpolationX(i);

        const double Up = 0.5 * (old_vel.u[xp]+old_vel.u[ic]); //main 
        const double Un = 0.5 * (old_vel.u[xm]+old_vel.u[ic]); //main

        const double Vp = old_vel.v[ic]+(old_vel.v[xp]-old_vel.v[ic])*interpolation;
        const double Vn = old_vel.v[ym]+(old_vel.v[ym+xp-ic]-old_vel.v[ym])*interpolation;

        const double Wp = old_vel.w[ic]+(old_vel.w[xp]-old_vel.w[ic])*interpolation;
        const double Wn = old_vel.w[zm]+(old_vel.w[zm+xp-ic]-old_vel.w[zm])*interpolation;

        // * ----------------------------------------------
        
#if defined(CONVECTION_DIFUSSION_QUICK)

        auto [uPx, uNx] = CalQuickSpanwise(Up, Un, xpp, xp, ic, xm, xmm, grid.dx,
                                           grid.staggered_dx, i, old_vel.u);
        auto [uPy, uNy] = CalQuickChordwise(Vp, Vn, ypp, yp, ic, ym, ymm, grid.dy,
                                            grid.staggered_dy, j, old_vel.u);
        auto [uPz, uNz] = CalQuickChordwise(Wp, Wn, zpp, zp, ic, zm, zmm, grid.dz,
                                            grid.staggered_dz, k, old_vel.u);

#elif defined(CONVECTION_DIFUSSION_LUD)

        auto [uPx, uNx] = CalLudSpanwise(Up, Un, xpp, xp, ic, xm, xmm, grid.dx,
                                         grid.staggered_dx, i, old_vel.u);
        auto [uPy, uNy] = CalLudChordwise(Vp, Vn, ypp, yp, ic, ym, ymm, grid.dy,
                                          grid.staggered_dy, j, old_vel.u);
        auto [uPz, uNz] = CalLudChordwise(Wp, Wn, zpp, zp, ic, zm, zmm, grid.dz,
                                          grid.staggered_dz, k, old_vel.u);

#elif defined(CONVECTION_DIFUSSION_UD)

        auto [uPx, uNx] =
            CalUpwind(Up, Un, xp, ic, xm, grid.dx, grid.staggered_dx, i, old_vel.u);
        auto [uPy, uNy] =
            CalUpwind(Vp, Vn, yp, ic, ym, grid.dy, grid.staggered_dy, j, old_vel.u);
        auto [uPz, uNz] =
            CalUpwind(Wp, Wn, zp, ic, zm, grid.dz, grid.staggered_dz, k, old_vel.u);

#endif

        auto convection  = - dt *((uPx * uPx - uNx * uNx) / grid.staggered_dx[i] +
                               (uPy * Vp - uNy * Vn) / grid.dy[j] +
                               (uPz * Wp - uNz * Wn) / grid.dz[k]);
        // * -------------- convection term --------------

        // ! -------------- diffusion term --------------

        // #ifdef TERBULENCE_SMAGORINSKY
        //     auto nu =  simu.GetNu() + old_vel.viseff(i, j, k);
        // #else
            auto nu = simu.GetNu();
        // #endif

         auto diffusion = nu*dt*(
                    ( (old_vel.u[xp]-old_vel.u[ic]) / grid.dx[i+1] 
                    - (old_vel.u[ic]-old_vel.u[xm]) / grid.dx[i]   ) / grid.staggered_dx[i]+
                    ( (old_vel.u[yp]-old_vel.u[ic]) / grid.staggered_dy[j]   
                    - (old_vel.u[ic]-old_vel.u[ym]) / grid.staggered_dy[j-1] ) / grid.dy[j]+
                    ( (old_vel.u[zp]-old_vel.u[ic]) / grid.staggered_dz[k]   
                    - (old_vel.u[ic]-old_vel.u[zm]) / grid.staggered_dz[k-1] ) / grid.dz[k]);
         // * -------------- diffusion term --------------

         curr_vel.u[ic] = old_vel.u[ic] + convection + diffusion;
    }

// ! ---------------------------------- y  ----------------------------------

    if(location.y_end == ny-gc){ d = one;}
    else{ d = 0;}

#pragma omp parallel for
        for (auto i = location.x_start; i < location.x_end ; ++i )
        for (auto j = location.y_start; j < location.y_end-d ; ++j )
        for (auto k = location.z_start; k < location.z_end ; ++k )
        { //  =================== v =================== //
            int ic = grid.icel(i,j,k);
            // ! -------------- convection term --------------
            const auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] 
                          = grid.Get12Nb(ic);
    
            const auto interpolation = grid.GetInterpolationY(j);

            const double Vp = 0.5 * (old_vel.v[yp]+old_vel.v[ic]); //main
            const double Vn = 0.5 * (old_vel.v[ym]+old_vel.v[ic]); //main

            const double Up = old_vel.u[ic]+(old_vel.u[yp]-old_vel.u[ic])*interpolation;
            const double Un = old_vel.u[xm]+(old_vel.u[xm+yp-ic]-old_vel.u[xm])*interpolation;

            const double Wp = old_vel.w[ic]+(old_vel.w[yp]-old_vel.w[ic])*interpolation;
            const double Wn = old_vel.w[zm]+(old_vel.w[zm+yp-ic]-old_vel.w[zm])*interpolation;

#if defined(CONVECTION_DIFUSSION_QUICK)

            auto [vPx, vNx] =
                CalQuickChordwise(Up, Un, xpp, xp, ic, xm, xmm, grid.dx,
                                  grid.staggered_dx, i, old_vel.v);
            auto [vPy, vNy] =
                CalQuickSpanwise(Vp, Vn, ypp, yp, ic, ym, ymm, grid.dy,
                                 grid.staggered_dy, j, old_vel.v);
            auto [vPz, vNz] =
                CalQuickChordwise(Wp, Wn, zpp, zp, ic, zm, zmm, grid.dz,
                                  grid.staggered_dz, k, old_vel.v);

#elif defined(CONVECTION_DIFUSSION_LUD)

            auto [vPx, vNx] = CalLudChordwise(Up, Un, xpp, xp, ic, xm, xmm,
                                              grid.dx, grid.staggered_dx, i, old_vel.v);
            auto [vPy, vNy] = CalLudSpanwise(Vp, Vn, ypp, yp, ic, ym, ymm,
                                             grid.dy, grid.staggered_dy, j, old_vel.v);
            auto [vPz, vNz] = CalLudChordwise(Wp, Wn, zpp, zp, ic, zm, zmm,
                                              grid.dz, grid.staggered_dz, k, old_vel.v);

#elif defined (CONVECTION_DIFUSSION_UD)

            auto [vPx, vNx] =
                CalUpwind(Up, Un, xp, ic, xm, grid.dx, grid.staggered_dx, i, old_vel.v);
            auto [vPy, vNy] =
                CalUpwind(Vp, Vn, yp, ic, ym, grid.dy, grid.staggered_dy, j, old_vel.v);
            auto [vPz, vNz] =
                CalUpwind(Wp, Wn, zp, ic, zm, grid.dz, grid.staggered_dz, k, old_vel.v);

#endif

            double convection =
                -dt * ((Up * vPx - Un * vNx) / grid.dx[i] +
                       (vPy * vPy - vNy * vNy) / grid.staggered_dy[j] +
                       (Wp * vPz - Wn * vNz) / grid.dz[k]);
            // * -------------- convection term --------------

            // ! -------------- diffusion term --------------

#ifdef TERBULENCE_SMAGORINSKY
            const auto nu = simu.GetNu() + old_vel.viseff[grid.icelCal(i, j, k)];
#else
            const auto nu = simu.GetNu();
#endif
            const double diffusion =
                nu * dt *
                (((old_vel.v[xp] - old_vel.v[ic]) / grid.staggered_dx[i] -
                  (old_vel.v[ic] - old_vel.v[xm]) / grid.staggered_dx[i - 1]) /
                     grid.dx[i] +
                 ((old_vel.v[yp] - old_vel.v[ic]) / grid.dy[j + 1] -
                  (old_vel.v[ic] - old_vel.v[ym]) / grid.dy[j]) /
                     grid.staggered_dy[j] +
                 ((old_vel.v[zp] - old_vel.v[ic]) / grid.staggered_dz[k] -
                  (old_vel.v[ic] - old_vel.v[zm]) / grid.staggered_dz[k - 1]) /
                     grid.dz[k]);
            // * -------------- diffusion term --------------

            curr_vel.v[ic] = old_vel.v[ic] + convection + diffusion;
        }

        if (location.z_end == nz - gc) {
        d = one;
        } else {
        d = 0;
        }

#pragma omp parallel for
        for (auto i = location.x_start; i < location.x_end; ++i)
        for (auto j = location.y_start; j < location.y_end; ++j)
        for (auto k = location.z_start; k < location.z_end - d;
             ++k) {  //  =================== w =================== //
            const int ic = grid.icel(i, j, k);
            // ! -------------- convection term --------------
            auto [xm, xp, ym, yp, zm, zp, xmm, xpp, ymm, ypp, zmm, zpp] =
                grid.Get12Nb(ic);

            const auto interpolation = grid.GetInterpolationZ(k);

            const double Wp = 0.5 * (old_vel.w[zp] + old_vel.w[ic]);  // main
            const double Wn = 0.5 * (old_vel.w[zm] + old_vel.w[ic]);  // main

            const double Up =
                old_vel.u[ic] + (old_vel.u[zp] - old_vel.u[ic]) * interpolation;
            const double Un =
                old_vel.u[xm] + (old_vel.u[xm + zp - ic] - old_vel.u[xm]) * interpolation;

            const double Vp =
                old_vel.v[ic] + (old_vel.v[zp] - old_vel.v[ic]) * interpolation;
            const double Vn =
                old_vel.v[ym] + (old_vel.v[ym + zp - ic] - old_vel.v[ym]) * interpolation;
            // *----------------------------------------------

#if defined(CONVECTION_DIFUSSION_QUICK)

            auto [wPx, wNx] = CalQuickChordwise(Up, Un, xpp, xp, ic, xm, xmm,
                                                grid.dx, grid.staggered_dx, i, old_vel.w);
            auto [wPy, wNy] = CalQuickChordwise(Vp, Vn, ypp, yp, ic, ym, ymm,
                                                grid.dy, grid.staggered_dy, j, old_vel.w);
            auto [wPz, wNz] = CalQuickSpanwise(Wp, Wn, zpp, zp, ic, zm, zmm,
                                               grid.dz, grid.staggered_dz, k, old_vel.w);

#elif defined(CONVECTION_DIFUSSION_LUD)

            auto [wPx, wNx] = CalLudChordwise(Up, Un, xpp, xp, ic, xm, xmm,
                                              grid.dx, grid.staggered_dx, i, old_vel.w);
            auto [wPy, wNy] = CalLudChordwise(Vp, Vn, ypp, yp, ic, ym, ymm,
                                              grid.dy, grid.staggered_dy, j, old_vel.w);
            auto [wPz, wNz] = CalLudSpanwise(Wp, Wn, zpp, zp, ic, zm, zmm,
                                             grid.dz, grid.staggered_dz, k, old_vel.w);

#elif defined(CONVECTION_DIFUSSION_UD)

            auto [wPx, wNx] =
                CalUpwind(Up, Un, xp, ic, xm, grid.dx, grid.staggered_dx, i, old_vel.w);
            auto [wPy, wNy] =
                CalUpwind(Vp, Vn, yp, ic, ym, grid.dy, grid.staggered_dy, j, old_vel.w);
            auto [wPz, wNz] =
                CalUpwind(Wp, Wn, zp, ic, zm, grid.dz, grid.staggered_dz, k, old_vel.w);

#endif

            const double convection =
                -dt * ((Up * wPx - Un * wNx) / grid.dx[i] +
                            (Vp * wPy - Vn * wNy) / grid.dy[j] +
                            (wPz * wPz - wNz * wNz) / grid.staggered_dz[k]);
            // * -------------- convection term --------------

            // * -------------- diffusion term --------------
#ifdef TERBULENCE_SMAGORINSKY
            auto nu = simu.GetNu() + old_vel.viseff[grid.icelCal(i, j, k)];
#else
            auto nu = simu.GetNu();
#endif

            const double diffusion =
                nu * dt *
                (((old_vel.w[xp] - old_vel.w[ic]) / grid.staggered_dx[i] -
                  (old_vel.w[ic] - old_vel.w[xm]) / grid.staggered_dx[i - 1]) /
                     grid.dx[i] +
                 ((old_vel.w[yp] - old_vel.w[ic]) / grid.staggered_dy[j] -
                  (old_vel.w[ic] - old_vel.w[ym]) / grid.staggered_dy[j - 1]) /
                     grid.dy[j] +
                 ((old_vel.w[zp] - old_vel.w[ic]) / grid.dz[k + 1] -
                  (old_vel.w[ic] - old_vel.w[zm]) / grid.dz[k]) /
                     grid.staggered_dz[k]);
            // * -------------- diffusion term --------------
            curr_vel.w[ic] = old_vel.w[ic] + convection + diffusion;
        }
}