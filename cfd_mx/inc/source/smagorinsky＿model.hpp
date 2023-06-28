#pragma once

#ifdef TERBULENCE_SMAGORINSKY

inline void GradientGauss(
    std::vector<double>& phi,
    std::vector<double>& dPhidX,
    std::vector<double>& Dxm,
    std::vector<double>& Dym,
    std::vector<double>& Dzm,
    LocalDomain& Lo,
    grid& gridA
){
    const auto [nx, ny, nz, gC] = gridA.nxyzgC;

    for (size_t i = Lo.x_start; i < Lo.x_end ; ++i )
    for (size_t j = Lo.y_start; j < Lo.y_end ; ++j )
    for (size_t k = Lo.z_start; k < Lo.z_end ; ++k )
    {
        const auto [xm, xp, ym, yp, zm, zp] = gridA.getNb6(i,j,k);
        const int ii = gridA.icelCal(i,j,k) *3;
        dPhidX[ii+0] = 0.5/Dxm[i] * (phi[xm] - phi[xp]);
        dPhidX[ii+1] = 0.5/Dym[j] * (phi[ym] - phi[yp]);
        dPhidX[ii+2] = 0.5/Dzm[k] * (phi[zm] - phi[zp]);
    }
}


void CalSmagorinskyModel(
    Simulation& simu,
    velocity& TX,
    pressure& t1,
    LocalDomain& Lo,
    grid& gridA
){
    const double CS = 0.18;

    if (TX.viseff.size() < gridA.iceltotCal){
        TX.viseff.resize(gridA.iceltotCal);
    }

    const auto [nx, ny, nz, gC] = gridA.nxyzgC;

    std::vector<double>  dudx(3*gridA.iceltotCal);
    std::vector<double>  dvdx(3*gridA.iceltotCal);
    std::vector<double>  dwdx(3*gridA.iceltotCal);

    const static auto square = [&](double a){return std::pow(a, 2);};


    GradientGauss(TX.u, dudx, gridA.Dxs, gridA.Dy, gridA.Dz, Lo, gridA);
    GradientGauss(TX.v, dvdx, gridA.Dx, gridA.Dys, gridA.Dz, Lo, gridA);
    GradientGauss(TX.w, dwdx, gridA.Dx, gridA.Dy, gridA.Dzs, Lo, gridA);

    for (size_t i = Lo.x_start; i < Lo.x_end ; ++i )
    for (size_t j = Lo.y_start; j < Lo.y_end ; ++j )
    for (size_t k = Lo.z_start; k < Lo.z_end ; ++k )
    {
        const auot icelCal = gridA.icelCal(i,j,k);
        const int t  = gridA.icelCal(i,j,k)*3;
        const double temp = std::abs(

                            2.0*(square(dudx[t+0])) + square(dvdx[t+1]) + square(dwdx[t+2])
                            + square(dudx[t+1] + dvdx[t+0])
                            + square(dudx[t+2] + dwdx[t+0])
                            + square(dvdx[t+2] + dwdx[t+1])

            );



        const double delta = std::pow(gridA.Dx[i] * gridA.Dy[j] * gridA.Dz[k], 1.0/3.0);

        TX.viseff[icelCal] = square(CS*delta) * std::sqrt(temp);
    }

}

#endif
