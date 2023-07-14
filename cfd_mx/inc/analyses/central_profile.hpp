#pragma once

#include <vector>

#include "backend/physical_variables.hpp"

void getCentProfile(
    Simulation& simu,
    StaggeredVelocity& vel,
    LocalDomain& Lo,
    StructuredGrid& grid
)
{
	static int io = 0;
    
    std::vector<double> y_plot(grid.cal_ny);
    std::vector<double> x_plot(grid.cal_nx);

    std::vector<double> y_plot_pos(grid.cal_ny);
    std::vector<double> x_plot_pos(grid.cal_nx);


    const int k1 = (grid.cal_nz)/2  + grid.no_ghost_cell;
    const int k2 = (grid.cal_nz)/2 - 1 + grid.no_ghost_cell;

    const int j2 = (grid.cal_ny)/2 - 1 + grid.no_ghost_cell;
    const int i2 = (grid.cal_nx)/2 - 1 + grid.no_ghost_cell;


    // y plot
    for (int j =  grid.no_ghost_cell;  j < grid.cal_ny+ grid.no_ghost_cell; ++j) {
        y_plot[j-grid.no_ghost_cell] = ( vel.U(i2, j, k1) + vel.U(i2, j, k2)) / 2.0;
        y_plot_pos[j-grid.no_ghost_cell] = grid.y_pos.at(j);
    }

    for (int i = grid.no_ghost_cell;  i < grid.cal_nx+ grid.no_ghost_cell; ++i) {
        x_plot[i-grid.no_ghost_cell] = (vel.V(i, j2, k1) + vel.V(i, j2, k2)) / 2.0;
        x_plot_pos[i-grid.no_ghost_cell] = grid.x_pos.at(i);
    }


    std::fstream file;
    std::string filename = "Information/"; // Information/Uprofile/
    // filename += std::to_string(io);
    filename += "x_plot_";
    filename += "re_";
    filename += std::to_string(int(simu.GetReynoldsNumber()));
    filename += ".dat";
    file.open(filename,std::ios::out);
    {
        for (size_t i = 0 ; i < x_plot.size(); i++) {
            file << x_plot_pos.at(i)  << '\t'<< x_plot[i] << std::endl;
        }
    }
    file.close();

    std::fstream file2;
    std::string filename2 = "Information/";

    // filename2 += std::to_string(io);
    filename2 += "y_plot_";
    filename2 += "re_";
    filename2 += std::to_string(int(simu.GetReynoldsNumber()));
    filename2 += ".dat";
    file2.open(filename2,std::ios::out);
    {
        for (size_t j = 0 ; j < y_plot.size(); j++) {
            file2 << y_plot_pos.at(j) << '\t' << y_plot[j] << std::endl;
        }
    }
    file2.close();
    io++;
}