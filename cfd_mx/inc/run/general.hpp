#pragma once

#ifndef RUN_GENERAL
#define RUN_GENERAL

#include "profiling/chronograph.hpp"

#include "resize_and_simuSetting.hpp"

#include "grid/Generator.hpp"

// --------------------------------
#include "analyses/cfl.hpp"

#include "analyses/CheckL2norm.hpp"

#include "analyses/CheckSteadystate.hpp"

#include "analyses/CentralProfile.hpp"
// --------------------------------

// --------------------------------
#include "source/engeryEquation.hpp"

#include "source/convection_and_difussion.hpp"

#include "source/convection_and_difussion.uniformGrid.hpp"

#include "source/smagorinskyModel.hpp"

// --------------------------------

#include "InitionCondition.hpp"

#include "pressure/Transform.hpp"

// --------------------------------
#include "io/Plot3D/qfile.hpp"

#include "io/Plot3D/xfileWrite.hpp"

#include "io/Plot3D/qfileRead.hpp"
// --------------------------------

// --------------------------------
#include "matrix/PressureMatrix.hpp"

#include "matrix/ELL_sparseMatrix.hpp"

#include "matrix/SPE_sparseMatrix0.hpp"

#include "matrix/SPE_sparseMatrix1.hpp"

#include "matrix/SPE_sparseMatrix.hpp"

#include "matrix/b_matrix.hpp"
// --------------------------------

// --------------------------------
#include "dfib/setEta.hpp"

#include "dfib/virtualForceIntergrator.hpp"

#include "dfib/updateUandF.hpp"
// --------------------------------

// --------------------------------
#include "solver/bicgstab.hpp"

#include "solver/bicgstab0.hpp"

#include "solver/bicgstabRe.hpp"

#include "solver/bicgstab_restart.hpp"
// --------------------------------


// --------------------------------
#include "BoundaryCondition/Boundary_Condition.hpp"

#include "BoundaryCondition/updateBC.hpp"
// --------------------------------

#include "pressure/solver/SorPipeLine.omp.hpp"

#include "io/ioTool.hpp"
#endif
