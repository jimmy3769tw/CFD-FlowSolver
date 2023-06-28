#pragma once

#include "virtualForceIntergrator.hpp"

auto getCD_CL_mpi(
    MPI_Comm & comm,
    Simulation& simu,
    LocalDomain& localDomain,
    ImmersedBoundary& Dfib,
    grid& gA
)
{
    virtualF_Int(Dfib, localDomain, gA, 0);

    virtualF_Int(Dfib, localDomain, gA, 1);

    auto area_cD    = gA.lz;
    auto area_cL    = gA.lz;

    auto cD = -2.0 * Dfib.ValSum[0] / area_cD;
 
    auto cL = -2.0 * Dfib.ValSum[1] / area_cL;


    double temp_g[2]{0.0, 0.0}, temp_l[2] = {cD, cL};
    MPI_Allreduce(temp_l, temp_g, 2, MPI_DOUBLE, MPI_SUM, comm);

    cD = temp_g [0];
    cL = temp_g [1];


    if (simu.PID == 0)
    {

    cout << "[cD, cL] : " << cD << ", " << cL << std::endl; 

    std::ofstream file;

    std::string name = "Information/Time_cDcL";

    name += ".dat";

    file.open (name, std::ios::out|ios::app);

    std::string tab = " ";

    if (simu.loop == 1 ){

        std::vector<std::string> variables;
        variables.push_back("simulation time");
        variables.push_back("C<sub>D</sub>");
        variables.push_back("C<sub>L</sub>");

            file 
            << "TITLE     = \"\"\n"
            << "VARIABLES = \""
            << variables.at(0)
            << "\",\""
            << variables.at(1)
            << "\",\""
            << variables.at(2)
            << "\"\n"
            << "ZONE T=\""
            << simu.ZONE()
            << "\"";
    }

    file << "\n"<<  simu.getSimuTime() << tab << cD << tab  << cL ;

    file.close();
    }


    return true;
}
