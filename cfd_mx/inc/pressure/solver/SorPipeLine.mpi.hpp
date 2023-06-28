#pragma once
# include "mpi_tool/init_mpi.hpp"

bool SorPipeLine_mpi_hybrid_omp(
    MPI_Comm comm_world,
    std::vector<int> &mpi_neighborhood,
    SORcoefficient& Sor,
    Simulation& simu,
    velocity& T1,
    pressure& t1,
    LocalDomain& Lo,
    grid& gA
){
    const int word_size = Lo.word_size;
    const auto [nx, ny, nz, gC] = gA.nxyzgC;
    double omega = simu.p_sor_omega;
    double mChangeMax;
    int itmax = simu.p_sor_iter_max;
    double rdt = 1.0 / simu.dt;
    if (simu.firstSOR) {
        simu.firstSOR = false;
        Sor.cf.resize(gA.iceltotCal * 8);
        #pragma omp parallel for firstprivate(Lo, rdt)
        for (auto i = Lo.x_start ; i < Lo.x_end ; ++i )
        for (auto j = Lo.y_start ; j < Lo.y_end ; ++j )
        for (auto k = Lo.z_start ; k < Lo.z_end ; ++k )
        {
            const int ii =  8*gA.icelCal(i,j,k);
            auto coef = gA.dy[j] * gA.dz[k] / gA.staggered_dx[i] +
                        gA.dy[j] * gA.dz[k] / gA.staggered_dx[i - 1] +
                        gA.dx[i] * gA.dz[k] / gA.staggered_dy[j] +
                        gA.dx[i] * gA.dz[k] / gA.staggered_dy[j - 1] +
                        gA.dx[i] * gA.dy[j] / gA.staggered_dz[k] +
                        gA.dx[i] * gA.dy[j] / gA.staggered_dz[k - 1];
            Sor.cf[ii  ] = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i]  ;
            Sor.cf[ii+1] = gA.Dy[j] * gA.Dz[k] / gA.Dxs[i-1];
            Sor.cf[ii+2] = gA.Dx[i] * gA.Dz[k] / gA.Dys[j]  ;
            Sor.cf[ii+3] = gA.Dx[i] * gA.Dz[k] / gA.Dys[j-1];
            Sor.cf[ii+4] = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k]  ;
            Sor.cf[ii+5] = gA.Dx[i] * gA.Dy[j] / gA.Dzs[k-1];
            Sor.cf[ii+7] = - 1.0 / coef;
        }
    }  


    double norm_rhs = 0;
    #pragma omp parallel for firstprivate(Lo, rdt) reduction(+ : norm_rhs)
    for (auto i = Lo.x_start ; i < Lo.x_end ; ++i )
    for (auto j = Lo.y_start ; j < Lo.y_end ; ++j )
    for (auto k = Lo.z_start ; k < Lo.z_end ; ++k )
    {
        const int icel = gA.icel(i,j,k);
        const int ii = gA.icelCal(i,j,k)*8;
        auto mChange = ( T1.u[icel] - T1.u[gA.icel(i-1,j,k)] ) * gA.Dy[j] * gA.Dz[k] 
                     + ( T1.v[icel] - T1.v[gA.icel(i,j-1,k)] ) * gA.Dx[i] * gA.Dz[k] 
                     + ( T1.w[icel] - T1.w[gA.icel(i,j,k-1)] ) * gA.Dx[i] * gA.Dy[j];
        Sor.cf[ii+6] = mChange * rdt;
        norm_rhs += std::pow( mChange * rdt, 2);
    }

    auto snetBuf_norm_rhs = sqrt(norm_rhs);
    MPI_Allreduce(&snetBuf_norm_rhs, &norm_rhs, 1, MPI_DOUBLE, MPI_SUM, comm_world);
    double interTime = omp_get_wtime();
    double residual;
    auto recvbuf = residual;
    for (simu.iters = 0;  simu.iters < simu.p_sor_iter_max; ++simu.iters) {
        // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~  MPI no-blocking send & recv ~~~~~~~~~~~~~~~~~~
        #ifdef MPI_DEBUG
        mpi_iSR_double_x_debugger(word_size, simu.TID, comm_world, mpi_neighborhood, t1.p, Lo, gA);
        #else
        mpi_iSR_double_x(1, comm_world, mpi_neighborhood, t1.p, Lo, gA );
        #endif
        MPI_Barrier(comm_world);
        // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~  MPI no-blocking send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~~
        residual = 0.0;
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
        auto sendbuf = residual;
        MPI_Allreduce(&sendbuf, &recvbuf,1, MPI_DOUBLE, MPI_SUM, comm_world);
        if (sqrt(recvbuf) < simu.p_criteria){break;}
    }

    simu.error =  sqrt(recvbuf);

    return true;
}
