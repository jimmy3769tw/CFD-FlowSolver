#pragma once

// Pressure_transform_X_result(t1, Mx, gridA);
void Pressure_transform_X_result_Dir(
    pressure &t1, 
    PressureMat &Mx ,
    LocalDomain& Lo,
    grid & gA
){
    #pragma omp parallel firstprivate(Lo)
    {
        #pragma omp for
        for (size_t i = Lo.x_start; i < Lo.x_end; ++i)
        for (size_t j = Lo.y_start; j < Lo.y_end; ++j)
        for (size_t k = Lo.z_start; k < Lo.z_end; ++k)
        {
            t1.p[gA.icel(i,j,k)] = Mx.x_result[gA.icelDir(i,j,k)];
        }
    }
}


std::tuple<double , double> getMax(
    const std::vector<double> &x,
    const LocalDomain& Lo,
    grid & gA
){

    double maxVal = - 1.0e10;
    double minVal = 1.0e10;

    #pragma omp parallel firstprivate(Lo) reduction(max : maxVal) reduction (min : minVal)
    {
        #pragma omp for
        for (auto i = Lo.x_start; i < Lo.x_end; ++i)
        for (auto j = Lo.y_start; j < Lo.y_end; ++j)
        for (auto k = Lo.z_start; k < Lo.z_end; ++k)
        {
            if (maxVal < x[gA.icelCal(i,j,k)]) 
            { maxVal =  x[gA.icelCal(i,j,k)];}


            if (minVal > x[gA.icelCal(i,j,k)]) 
            { minVal =  x[gA.icelCal(i,j,k)];}
        }
    }

    return std::make_pair(maxVal, minVal);

}



#ifdef EIGEN_ON

void Pressure_transform_x_Eigen(
    pressure &t1, 
    PressureMat &Mx ,
    LocalDomain& Lo,
    grid & gridA
){
    const auto [nx, ny , nz] = gridA.nxyz;

    for (size_t i = Lo.x_start; i < Lo.x_end; ++i)
    for (size_t j = Lo.y_start; j < Lo.y_end; ++j)
    for (size_t k = Lo.z_start; k < Lo.z_end; ++k)
    {
        t1.p[gridA.icel(i,j,k)] = Mx.x_Eigen[gridA.icelCal(i,j,k)];
    }
}

#endif