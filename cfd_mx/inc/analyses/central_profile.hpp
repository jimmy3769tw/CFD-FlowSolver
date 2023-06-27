#pragma once

void getCentProfile(
    Simulation& simu,
    velocity& T3,
    CalDomain& Lo,
    grid& gridA
)
{
    auto [nx, ny, nz , gC] = gridA.nxyzgC;

	int io = simu.currentFile;

        vector<double> uc(ny-4);
        vector<double> vc(nx-4);


        const int k1 = (gridA.nz - 4)/2;
        const int k2 = (gridA.nz - 4)/2 - 1;
		for (size_t j = 0;  j < ny-4; ++j) 
		{
            const int i = (gridA.nx - 4)/2;

            const auto icelk1 = gridA.icel(i+2, j+2, k1+2);

            const auto icelk2 = gridA.icel(i+2, j+2, k2+2);

            const double buffer1 = (T3.u[icelk1-(ny*nz)] + T3.u[icelk1]);

            const double buffer2 = (T3.u[icelk2-(ny*nz)] + T3.u[icelk2]);
			uc[j] = 0.125 * (buffer1 + buffer2);
		}

		for (size_t i = 0; i < nx-4; ++i)
		{
			// const int icelShift = (i+2)*nz*ny + (j+2)*nz + (k+2);
            const int j = (gridA.ny - 4)/2;
            const auto icelk1 = gridA.icel(i+2, j+2, k1+2);

            const auto icelk2 = gridA.icel(i+2, j+2, k2+2);

            const double bufferk1 = (T3.v[icelk1-(nz)] + T3.v[icelk1]);

            const double bufferk2 = (T3.v[icelk2-(nz)] + T3.v[icelk2]);
			vc[i] = 0.25 * (bufferk1 + bufferk2);
		}

        std::fstream file;
        std::string filename = "Information/Uprofile/"; // Information/Uprofile/
		filename += std::to_string(io);
        filename += "_Uprofile";
        int Reint = simu.Re;
		filename += std::to_string(Reint);
		filename += ".dat";
        file.open(filename,std::ios::out);
        {
            file << "TITLE =\"present (u)\"" << std::endl;
            file << "VARIABLES = \"Location(y)\",\"u-Velocity\""<< std::endl;
            file << "ZONE T=\"present\"\t i= "<< nx-4 << std::endl;
            for (size_t i =0 ; i < nx -4;i++)
            {
                file << gridA.Xc[i]  << '\t'<< uc[i] << std::endl;
            }
        }
        file.close();

        std::fstream file2;
        std::string filename2 = "Information/Vprofile/";

		filename2 += std::to_string(io);
        filename2 += "_Vprofile";
		filename2 += "Re_";
		filename2 += std::to_string(Reint);
		filename2 += ".dat";
        file2.open(filename2,std::ios::out);
        {
            file2 << "TITLE =\"present (v)\"" << std::endl;
            file2 << "VARIABLES = \"Location(x)\",\"v-Velocity\""<< std::endl;
            file2 << "ZONE T=\"present\"\t i= "<< ny-4 << std::endl;
            for (size_t j =0 ; j < ny -4;j++)
            {
                file2 << gridA.Yc[j] << '\t' << vc[j] << std::endl;
            }
        }
        file2.close();
}