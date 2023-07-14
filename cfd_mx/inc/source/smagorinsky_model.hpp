#pragma once
#include "backend/simulation.hpp"
#include "backend/physical_variables.hpp"
#include "grid/structured_grid.hpp"
#include <vector>



void CalSmagorinskyModel(
    Simulation& simu,
    StaggeredVelocity& vel,
    LocalDomain& Lo,
    StructuredGrid& grid
){
    constexpr double CS = 0.18;
    const int cal_size = grid.cal_no_grid;

    std::vector<double>  du_dx(3*cal_size);
    std::vector<double>  dv_dx(3*cal_size);
    std::vector<double>  dw_dx(3*cal_size);

    const static auto square = [&](double a){return std::pow(a, 2);};

    for (size_t i = Lo.x_start; i < Lo.x_end ; ++i )
    for (size_t j = Lo.y_start; j < Lo.y_end ; ++j )
    for (size_t k = Lo.z_start; k < Lo.z_end ; ++k )
    {
        const int ii = grid.icelCal(i,j,k)*3;
        du_dx[ii+0] = (vel.U(i-1, j, k) - vel.U(i, j, k)) / grid.staggered_dx[i];
        du_dx[ii+1] = (vel.U(i, j-1, k) - vel.U(i, j, k)) / grid.dy[j] ;
        du_dx[ii+2] = (vel.U(i, j, k-1) - vel.U(i, j, k)) / grid.dz[j] ;
        
        dv_dx[ii+0] = (vel.V(i-1, j, k) - vel.V(i, j, k)) / grid.dx[j] ;
        dv_dx[ii+1] = (vel.V(i, j-1, k) - vel.V(i, j, k)) / grid.staggered_dy[j] ;
        dv_dx[ii+2] = (vel.V(i, j, k-1) - vel.V(i, j, k)) / grid.dz[j] ;
        
        dw_dx[ii+0] = (vel.W(i-1, j, k) - vel.W(i, j, k)) / grid.dx[j] ;
        dw_dx[ii+1] = (vel.W(i, j-1, k) - vel.W(i, j, k)) / grid.dy[j] ;
        dw_dx[ii+2] = (vel.W(i, j, k-1) - vel.W(i, j, k)) / grid.staggered_dz[k] ;
    }

    for (size_t i = Lo.x_start; i < Lo.x_end ; ++i )
    for (size_t j = Lo.y_start; j < Lo.y_end ; ++j )
    for (size_t k = Lo.z_start; k < Lo.z_end ; ++k )
    {
        const int t  = grid.icelCal(i,j,k)*3;
        const double temp = std::abs(
                2.0*(square(du_dx[t+0])) + square(dv_dx[t+1]) + square(dw_dx[t+2])
                + square(du_dx[t+1] + dv_dx[t+0])
                + square(du_dx[t+2] + dw_dx[t+0])
                + square(dv_dx[t+2] + dw_dx[t+1])
            );
        const double delta = std::pow(grid.dx[i] * grid.dy[j] * grid.dz[k], 1.0/3.0);
        vel.VisEff(i, j, k) = square(CS*delta) * std::sqrt(temp);
    }

}
