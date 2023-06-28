#pragma once


#include "omp.h"

bool SorPipeLine_omp(
    SORcoefficient& Sor,
    Simulation& simu,
    velocity& T1,
    pressure& t1,
    LocalDomain& Lo,
    grid& gA
)
{
    const auto [nx, ny, nz, gC] = gA.nxyzgC;
    double omega = simu.p_sor_omega;
    double mChangeMax;
    int itmax = simu.p_sor_iter_max;
    double rdt = 1.0 / simu.dt;
    if (simu.firstSOR) {
        simu.firstSOR == false;
        Sor.cf.resize(gA.iceltotCal * 8);
        #pragma omp parallel for firstprivate(Lo)
        for (auto i = Lo.x_start ; i < Lo.x_end ; ++i )
        for (auto j = Lo.y_start ; j < Lo.y_end ; ++j )
        for (auto k = Lo.z_start ; k < Lo.z_end ; ++k )
        {
          const int ii = 8 * gA.icelCal(i, j, k);
          auto coef = gA.dy[j] * gA.dz[k] / gA.staggered_dx[i] +
                      gA.dy[j] * gA.dz[k] / gA.staggered_dx[i - 1] +
                      gA.dx[i] * gA.dz[k] / gA.staggered_dy[j] +
                      gA.dx[i] * gA.dz[k] / gA.staggered_dy[j - 1] +
                      gA.dx[i] * gA.dy[j] / gA.staggered_dz[k] +
                      gA.dx[i] * gA.dy[j] / gA.staggered_dz[k - 1];
          Sor.cf[ii] = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i];
          Sor.cf[ii + 1] = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i - 1];
          Sor.cf[ii + 2] = gA.Dx[i] * gA.Dz[k] / gA.Dys[j];
          Sor.cf[ii + 3] = gA.Dx[i] * gA.Dz[k] / gA.Dys[j - 1];
          Sor.cf[ii + 4] = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k];
          Sor.cf[ii + 5] = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k - 1];
          Sor.cf[ii + 7] = -1.0 / coef;
        }
    }  

    double norm_rhs = 0.0;
    #pragma omp parallel for firstprivate(Lo, rdt) reduction(+ : norm_rhs)
    for (auto i = Lo.x_start ; i < Lo.x_end ; ++i )
    for (auto j = Lo.y_start ; j < Lo.y_end ; ++j )
    for (auto k = Lo.z_start ; k < Lo.z_end ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const int ii8 = gA.icelCal(i,j,k)*8;
        auto mChange = ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] ) * gA.Dy[j] * gA.Dz[k] 
                     + ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] ) * gA.Dx[i] * gA.Dz[k] 
                     + ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] ) * gA.Dx[i] * gA.Dy[j] ;

        Sor.cf[ii8+6] = mChange * rdt;
        norm_rhs += std::pow( mChange * rdt, 2);
    }

    norm_rhs = sqrt(norm_rhs);

    double interTime = omp_get_wtime();

    // * --------- setting Peridic and Neunann --------- 

    double residual;
    for (simu.iters = 0;  simu.iters < simu.p_sor_iter_max; ++simu.iters)
    {
        residual = 0;
        #pragma omp parallel for reduction(+:residual) firstprivate(Lo, norm_rhs)
        for (auto i = Lo.x_start ; i < Lo.x_end ; ++i )
        for (auto j = Lo.y_start ; j < Lo.y_end ; ++j )
        for (auto k = Lo.z_start ; k < Lo.z_end ; ++k )
        {
            const int icel = gA.icel(i,j,k);
            const int ii8 = gA.icelCal(i,j,k)*8;

            auto pNEW =(- t1.p[gA.icel(i+1,j,k)] * Sor.cf[ii8  ]
                        - t1.p[gA.icel(i-1,j,k)] * Sor.cf[ii8+1]
                        - t1.p[gA.icel(i,j+1,k)] * Sor.cf[ii8+2]
                        - t1.p[gA.icel(i,j-1,k)] * Sor.cf[ii8+3]
                        - t1.p[gA.icel(i,j,k+1)] * Sor.cf[ii8+4]
                        - t1.p[gA.icel(i,j,k-1)] * Sor.cf[ii8+5]
                + Sor.cf[ii8+6]) * Sor.cf[ii8+7];

            double pChange = std::abs(pNEW - t1.p[icel]);
            t1.p[icel] += omega * (pNEW - t1.p[icel]);
            residual += pChange/norm_rhs;
        }

        if (sqrt(residual) < simu.p_criteria){break;}
    }

    simu.error = sqrt(residual);
    return true;
}
