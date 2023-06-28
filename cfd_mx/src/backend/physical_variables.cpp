#include "backend/physical_variables.hpp"

Pressure::Pressure(StructuredGrid& grid) : PhysicalVal(grid){
    p.resize(grid.no_grid);
}

StaggeredVelocity::StaggeredVelocity(StructuredGrid &grid) :
    PhysicalVal(grid) {
    resize(grid_->no_grid);
    viseff_.resize(grid_->cal_no_grid);
}


ImmersedBoundary::ImmersedBoundary(StructuredGrid& grid):
    PhysicalVal(grid) {
    eta.resize(grid_->no_grid);
    f.resize(grid_->cal_no_grid * 3);
    yShift_ = grid_->cal_no_grid;
    zShift_ = 2 * grid_->cal_no_grid;
}

