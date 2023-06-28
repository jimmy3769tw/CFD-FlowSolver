#include "backend/domain.hpp"
#include "grid/structured_grid.hpp"

void GlobalDomain::DivideDomain(int dims, int n) {
   if (dims == 0) {
     throw std::invalid_argument("DivideDomain !!!!\n");
   }
   resize(dims);
   DivideLen(dims, n);
   start[0] = grid_->no_ghost_cell;
   for (int i = 1; i < dims; i++) {
     start[i] = start[i - 1] + len[i - 1];
   }

   for (int i = 0; i < dims; i++) {
     end[i] = start[i] + len[i];
   }
}


void LocalDomain::Init(std::vector<int> dims, std::vector<int> rank) {
    global_x_.DivideDomain(dims[0], grid_->cal_nx);
    global_y_.DivideDomain(dims[1], grid_->cal_ny);
    global_z_.DivideDomain(dims[2], grid_->cal_nz);
    InitLocal(rank[0], rank[1], rank[2]);
}


LocalDomain::LocalDomain(StructuredGrid& grid): grid_(&grid),
        global_x_(GlobalDomain{grid}),
        global_y_(GlobalDomain{grid}),
        global_z_(GlobalDomain{grid}) {}


std::pair<std::vector<int>, std::vector<int> > LocalDomain::GetOneDimStartEnd() {
int slid = grid_->cal_ny_times_nz;
auto beg = global_x_.start;
auto end = global_x_.end;

for (auto &x : beg) {
    x -= grid_->no_ghost_cell;
}
for (auto &x : end) {
    x -= grid_->no_ghost_cell;
}
for (auto &x : beg) {
    x *= slid;
}
for (auto &x : end) {
    x *= slid;
}

return std::make_pair(beg, end);
}