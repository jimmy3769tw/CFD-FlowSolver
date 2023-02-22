#pragma once 
#include"controlPanel.hpp"


double CheckL2Norm(
    Simulation& simu,
    const velocity& A,
    const velocity& B,
    CalDomain& Lo,
    grid& gA
){
    auto [nx, ny, nz , gC] = gA.nxyzgC;

    double temp = 0.0;
    double L2norm = 0.0;
    int icel;

    for (size_t i = Lo.x_start; i < Lo.x_end ; ++i )
    for (size_t j = Lo.y_start; j < Lo.y_end ; ++j )
    for (size_t k = Lo.z_start; k < Lo.z_end ; ++k )
    {
        icel = gA.icel(i,j,k);
        temp += pow((A.u[icel] - B.u[icel]),2);
    }

    for (size_t i = Lo.x_start; i < Lo.x_end ; ++i )
    for (size_t j = Lo.y_start; j < Lo.y_end ; ++j )
    for (size_t k = Lo.z_start; k < Lo.z_end ; ++k )
    {
        icel = gA.icel(i,j,k);
        temp += pow((A.v[icel] - B.v[icel]),2);
    }

    for (size_t i = Lo.x_start; i < Lo.x_end ; ++i )
    for (size_t j = Lo.y_start; j < Lo.y_end ; ++j )
    for (size_t k = Lo.z_start; k < Lo.z_end ; ++k )
    {
        icel = gA.icel(i,j,k);
        temp += pow((A.w[icel] - B.w[icel]),2);
    }
    L2norm = sqrt( temp / ((gA.nx - 4) * (gA.ny - 4) * (gA.nz - 4)));

    return L2norm;
}