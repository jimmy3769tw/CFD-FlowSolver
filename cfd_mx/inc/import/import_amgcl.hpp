#pragma once

#include <amgcl/adapter/crs_tuple.hpp>  // for build crs_tuple
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>  // for build OpenMP
#include <amgcl/io/mm.hpp>

// ! make_solver
// 1.1
#include <amgcl/make_solver.hpp>
// 1.2
#include <amgcl/coarsening/plain_aggregates.hpp>
#include <amgcl/coarsening/rigid_body_modes.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/lgmres.hpp>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>


//   the solver backend:
typedef amgcl::backend::builtin<double> SBackend;
//   the preconditioner backend:
typedef amgcl::backend::builtin<double> PBackend;

//=========== Compose the solver type ===========//

typedef amgcl::make_solver<
    amgcl::amg<PBackend, amgcl::coarsening::smoothed_aggregation,
               amgcl::relaxation::gauss_seidel>,
    amgcl::solver::bicgstab<SBackend> >
    AmgclSolver;
