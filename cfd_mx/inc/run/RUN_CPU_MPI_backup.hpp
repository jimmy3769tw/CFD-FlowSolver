# include "controlPanel.hpp"

# include "run/general.hpp"
# include "backend/print_openmp.hpp"
# include <mpi.h>

# include "run/RUN_CPU.hpp"
# include "mpi_tool/init_mpi.hpp"
# include "dfib/virtualForceIntergrator_mpi.hpp"
# include "mpi_tool/pointTOpoint/blocking_p2p.hpp"
# include "mpi_tool/pointTOpoint/nonBlocking_p2p.hpp"
# include "pressure/solver/SorPipeLine.mpi.hpp"
# include "solver/mpi/bicgstab_restart.hpp"

// #define MPI_DEBUG

bool Run_Hy_MPI_OpenMP(
    grid &gA,
    ClockStruct &timer,
    Simulation &simu,
    int argc, char **argv
){

  timer.beginNew.start();

  // Global variable
  // ImmersedBoundary Dfib;
  // pressure t1;
  // velocity T0, T1, T3;
  // SORcoefficient Sor;
  // PressureMat Mx;
  // Local variable
  ImmersedBoundary Dfib_loc;
  pressure t1_loc;
  velocity T0_loc, T1_loc, T3_loc;
  SORcoefficient Sor_loc;
  PressureMat Mx_loc;

  PrintIsOpenmpExist();
  int ompThreads = omp_get_max_threads();  // Using Inv OMP_NUM_THREADS
  omp_set_num_threads(ompThreads);
  // ~~~~~~~~~~~~~~~~~~~~~~~~ openMP ~~~~~~~~~~~~~~~~~~~~~~~~

  // ! ============================  divid Domain ============================
  CalDomain globalDomain;
  CalDomain localDomain;
  std::vector<int> grid_size{gA.nx, gA.ny, gA.nz};
  std::vector<int> dims_loc{2};   //! change here to use more process: ex: mpirun -np 2 ./mx
                                  //! NOW ONLY SUPORT ONE DIMENTION DIVISUN
  globalDomain.InitTable(grid_size, {1, 1, 1});
  globalDomain.InitLocal(0, 0, 0);

  bool reorder = true;

  auto [mx_comm_world, 
        mpi_word_size, 
        mpi_word_rank, 
        mpi_coord, 
        mpi_neighborhood] = mpi_init(reorder, dims_loc, argc, argv);

  for (size_t i = 0 ; i < (4 - dims_loc.size()) ; ++i){
    dims_loc.push_back(1);
    mpi_coord.push_back(0);
  }

  localDomain.InitTable(grid_size, dims_loc);
  localDomain.InitLocal(mpi_coord.at(0), mpi_coord.at(1), mpi_coord.at(2));

  // * ========================================================================

  // check CalDomain
  simu.PID = mpi_word_rank;

  if (simu.PID == 0)
  {
    if(opendir("Information") == NULL)
    { if (system("mkdir Information") != 0){ return 1; } }

    if(opendir("mx_out") == NULL)
    { if (system("mkdir mx_out") != 0){ return 1; } }

    if(opendir("Information/Chronograph") == NULL)
    { if (system("mkdir Information/Chronograph") != 0){ return 1; } }
  }


  for(size_t i = 0; i < mpi_word_size ; ++i){
    MPI_Barrier(mx_comm_world);
    if (mpi_word_rank == i){
      cout  << "[Rank: {begin,endof}] = "<< mpi_word_rank 
            << ": [{ "  << localDomain.x_start
            << ", "     << localDomain.x_end
            << "},{ "   << localDomain.y_start
            << ", "     << localDomain.y_end
            << "},{ "   << localDomain.z_start
            << ", "     << localDomain.z_end
            << "}]"     << endl;
    }
  }

  MPI_Barrier(mx_comm_world);
  
  if (simu.PID == 0) cout << "\n<i_begin_table>,<x_endof_table>";
  
  for(size_t j = 0; j < mpi_word_size ; ++j){
    cout << std::flush;
    MPI_Barrier(mx_comm_world);

    if (mpi_word_rank == j){
      cout << "\n===== " << j << " =====" << endl;

      for(size_t i = 0; i < mpi_word_size ; ++i){
        cout << "{" << localDomain.i_begin_table.at(i) << ", ";
        cout << localDomain.x_endof_table.at(i) << "}, ";
      }

    }
  }

  cout << std::flush;
  MPI_Barrier(mx_comm_world);
  if (simu.PID == 0) cout << "\n<i_length_table>";
  for(size_t j = 0; j < mpi_word_size ; ++j){
    MPI_Barrier(mx_comm_world);
    if (mpi_word_rank == j){
      cout << "\n===== " << j << " =====" << endl;
      for(size_t i = 0; i < mpi_word_size; ++i){
        cout << localDomain.i_length_table.at(i) << ", ";
      }
    }
  }

  // for (auto &x:T0.U(0)){
  //   x = simu.PID;
  // }


  // for(size_t j = 0; j < mpi_word_size ; ++j)
  // {
  //   MPI_Barrier(mx_comm_world);
  //   if (mpi_word_rank == j){
  //     cout << "\n===== " << j << " =====" << endl;
  //      T0.showVelt(0, gA);
  //   }
  // }

  // mpi_Bcast(mx_comm_world, gA, localDomain, T0.U(0));
  // mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, T0.U(0), localDomain, gA );
  // mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, T0.U(0), localDomain, gA );

  // for(size_t j = 0; j < mpi_word_size ; ++j)
  // {
  //   MPI_Barrier(mx_comm_world);
  //   if (mpi_word_rank == j){
  //     cout << "\n===== " << j << " =====" << endl;
  //      T0.showVelt(0, gA);
  //   }
  // }


  // return 0;

  // resize_variable(gA, t1, T0, T1, T3, Dfib); 
  resize_variable(gA, t1_loc, T0_loc, T1_loc, T3_loc, Dfib_loc); 
  {
    double ui = 0.0, vi = 0.0, wi = 0.0;

    // T0.iniU_omp(ui, vi, wi);
    // T1.iniU_omp(ui, vi, wi);
    // T3.iniU_omp(ui, vi, wi);
    // t1.init_p(0.0);

    T0_loc.iniU_omp(ui, vi, wi);
    T1_loc.iniU_omp(ui, vi, wi);
    T3_loc.iniU_omp(ui, vi, wi);
    t1_loc.init_p(0.0);
  }


  generateGride(simu, ShareM, globalDomain, gA);

  if (simu.PID == 0) {
    gA.io_csv("Information/gA.csv");
    OutputPlot3D_Xfile(simu, gA);
  } 

  // ! ## The first BC -------------
  // T1.u = T3.u = T0.u ;
  // T1.v = T3.v = T0.v ;
  // T1.w = T3.w = T0.w ;

  T1_loc.u = T3_loc.u = T0_loc.u;
  T1_loc.v = T3_loc.v = T0_loc.v;
  T1_loc.w = T3_loc.w = T0_loc.w;

  // UpdateAllVelocityOnBoundary( globalDomain, T0, t1, gA );
  // BC_staggered_copy( globalDomain, T0, T1, gA );

  UpdateAllVelocityOnBoundary( localDomain, T0_loc, t1_loc, gA);
  BC_staggered_copy( localDomain, T0_loc, T1_loc, gA);

  // double accumlate = 0;
  // for (int i = 0; i < T1.u.size(); i++) {
  //   accumlate += std::abs(T1.u[i]- T1.lo)
  // }
  // !## prepare for poisson equation --------------------------
  timer.pressure.start();


  #if  defined (P_SOLVER_ELL) \
    || defined (P_SOLVER_CSR) \
    || defined (P_SOLVER_SPE) \
    || defined (P_SOLVER_AMGCL_BUILTIN) \
    // Mx.x_result.resize(gA.iceltotCal, 0.0);
    Mx_loc.x_result.resize(gA.iceltotCal, 0.0);
  #endif



  // !### ---------------------------------------- ELL (BICG)
  #if defined (P_SOLVER_ELL)

  // ---------------------------
  // createPressureMatrix(Mx, simu, localDomain, gA);
  // ---------------------------

  // ---------------------------
  // solver::BicgstabRestart<ELL_type> pSolver(Mx.matA_ell);
  // ---------------------------


  // ---------------------------
  createPressureMatrix(Mx_loc, simu, globalDomain, gA);
  // ---------------------------
  // ---------------------------
  auto [beg, end] = localDomain.GetStartEnd(gA.nyCal * gA.nzCal);
  Mx_loc.matA_ell.mpi_init(mx_comm_world, gA.iceltotCal, beg, end);
  // ---------------------------

  // ---------------------------
  solver::bicgstab_restart<ELL_type> pSolver_loc(Mx_loc.matA_ell);
  // ---------------------------
  // !### ---------------------------------------- CSR
  #elif defined (P_SOLVER_CSR)

  //  ----------------------------------------------------------------
  // Mx.matA_csr.set( gA.createPressureMatrixCSR() );
  //  ----------------------------------------------------------------

  //  ----------------------------------------------------------------
  // solver::BicgstabRestart<CSR_type> pSolver(Mx.matA_csr);
  //  ----------------------------------------------------------------


  //  ----------------------------------------------------------------
  Mx_loc.matA_csr.set( gA.createPressureMatrixCSR() );
  //  ----------------------------------------------------------------
  //  ----------------------------------------------------------------
  {
  auto [beg, end] = localDomain.GetStartEnd(gA.nyCal * gA.nzCal);
  Mx_loc.matA_csr.mpi_init(mx_comm_world, gA.iceltotCal, beg,  end);
  }
  //  ----------------------------------------------------------------
  //  ----------------------------------------------------------------
  solver::bicgstab_restart<CSR_type> pSolver_loc(Mx_loc.matA_csr);
  //  ----------------------------------------------------------------


  // !### ---------------------------------------- SPE
  #elif defined(P_SOLVER_SPE)

  // --------------------------------------
  // Mx.matA_spe.resize(gA.nxCal, gA.nyCal, gA.nzCal);
  // --------------------------------------

  // --------------------------------------
  // Mx.matA_spe.setupPressure( gA.gC, gA.Dx, gA.Dy, gA.Dz,gA.Dxs, gA.Dys, gA.Dzs);
  // --------------------------------------

  // #if defined(JACOBI_PC)
  // gA.jp =  Mx.matA_spe.set_Jp_pc(); 
  // #endif


  // --------------------------------------
  // solver::BicgstabRestart<SPE_type> pSolver(Mx.matA_spe);
  // --------------------------------------


  // --------------------------------------
  Mx_loc.matA_spe.resize(gA.nxCal, gA.nyCal, gA.nzCal);
  // --------------------------------------
  Mx_loc.matA_spe.setupPressure(
    gA.gC, gA.Dx, gA.Dy, gA.Dz, gA.Dxs, gA.Dys, gA.Dzs);
  // --------------------------------------
  {
  auto [beg, end] = localDomain.GetStartEnd(gA.nyCal * gA.nzCal);
  Mx_loc.matA_spe.mpi_init(mx_comm_world, gA.iceltotCal, beg, end);
  }
  // --------------------------------------

  // --------------------------------------
  #if defined(JACOBI_PC)
  gA.jp =  Mx_loc.matA_spe.set_Jp_pc();
  #endif
  // --------------------------------------

  //  ----------------------------------------------------------------
  solver::bicgstab_restart<SPE_type> pSolver_loc(Mx_loc.matA_spe);
  //  ----------------------------------------------------------------
  #endif


  #ifndef JACOBI_PC
  // gA.jp.clear();
  #endif


  #if  defined(P_SOLVER_ELL) \
    || defined(P_SOLVER_CSR) || defined(P_SOLVER_SPE)

    // pSolver.setTolerance(simu.p_criteria);
    pSolver_loc.setTolerance(simu.p_criteria);

  #endif  // P_SOLVER_ELL  P_SOLVER_CSR  P_SOLVER_SPE

  timer.pressure.stop();
  // * ----------------------------------------------------

  std::fill(Dfib_loc.eta.begin(), Dfib_loc.eta.end(), 0.0);
  if (simu.DfibMethod == "OFF") {}
  else if (simu.DfibMethod == "DFIB_Cylinder-X")
    CalEtaCylinderXdirection(Dfib_loc, globalDomain, gA);
  else if (simu.DfibMethod == "DFIB_Cylinder-Z")
    DFIB_CylinderZ(Dfib_loc, globalDomain, gA);
  else
    throw std::invalid_argument("DFIB method??");
  // * ------------------------------------------------------
    WriteQfile(Dfib_loc, simu, t1_loc, T0_loc, gA);

    #ifdef TERBULENCE_SMAGORINSKY
    // T0.viseff.resize(gA.iceltotCal, 0.0);
    T0_loc.viseff.resize(gA.iceltotCal, 0.0);
    T1_loc.viseff.resize(gA.iceltotCal);
    #endif


  auto perLoopWTime = timer.beginNew.elapsedTime();
  // ! ==================  main loop begin ==================
  for (simu.loop = 1; simu.get_finishloop(); simu.AddLoop()) {
    // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (int whichDim = 0; whichDim < 3; whichDim++) {
      mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood,
      T0_loc.U(whichDim), localDomain, gA);
    }
    // * ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~

    // !## 1. get the T1/u* --------------------
    timer.convectionDifussion.start();
    // ---------------------------------
    #if defined (TERBULENCE_SMAGORINSKY)
    // CalSmagorinskyModel(ShareM, simu, T0, t1, globalDomain, gA);
    CalSmagorinskyModel(ShareM, simu, T0_loc, t1_loc, localDomain, gA);
    #endif
    // ---------------------------------

    // ---------------------------------
    // ConvectionDifussion(simu, T0, T1, globalDomain, gA);
    ConvectionDifussion(simu, T0_loc, T1_loc, localDomain, gA);
    // ---------------------------------

    // ------- the error could be happen if you call BC_updateSlid().
    // BC_updateSlid(globalDomain, T1, gA);
    // BC_updateSlid(localDomain, T1_loc, gA);
    // ---------------------------------
    timer.convectionDifussion.stop();
    //*---------------------------------------------------

    // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
    for (int whichDim = 0; whichDim < 3 ; whichDim++) {
      mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood,
      T1_loc.U(whichDim), localDomain, gA);
    }
    // * ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
    // -----------------------
    // get_Cfl(T1, globalDomain, gA, simu);
    // -----------------------

    // -----------------------
    // get_Cfl(T1_loc, localDomain, gA, simu);
    // -----------------------



    // !## 2. Get the pressure poisson equation's result.
    timer.pressure.start();


    // !### --------------------- SOR method
    int pLoopIni  = 1;
    #if defined (P_SOLVER_SOR)
    if (simu.loop > pLoopIni) {
    // ---------------------------------
    // SorPipeLine_omp(Sor, simu, T1, t1, globalDomain, gA);
    SorPipeLine_mpi_hybrid_omp(mx_comm_world,
    mpi_neighborhood, Sor_loc, simu, T1_loc, t1_loc, localDomain, gA);
    // ---------------------------------
    }

    // !### ----------------------- CSR/SPE/ELL (BICG)(CG) / AMGCL
    #elif defined(P_SOLVER_BICG_CSR) \
      ||  defined(P_SOLVER_BICG_SPE) \
      ||  defined(P_SOLVER_BICG_ELL) \
      ||  defined(P_SOLVER_AMGCL_BUILTIN)


    if (simu.loop > pLoopIni) {
      // ---------------------------------
      // CalMatB( T1, Mx, simu, globalDomain, gA);
      // ---------------------------------

      // ---------------------------------
      // std::tie(simu.iters, simu.error) = pSolver(Mx.mat_b, Mx.x_result);
      // ---------------------------------

      // ---------------------------------
      // Pressure_transform_X_result(t1, Mx, globalDomain, gA);
      // ---------------------------------

      // ---------------------------------
      CalMatB(T1_loc, Mx_loc, simu, localDomain, gA);
      mpi_Bcast(mx_comm_world, gA, localDomain, Mx_loc.mat_b);
      // ---------------------------------

      // ---------------------------------
      std::tie(simu.iters, simu.error) = pSolver_loc(Mx_loc.mat_b, Mx_loc.x_result);
      // ---------------------------------

      // ---------------------------------
      Pressure_transform_X_result(t1_loc, Mx_loc, localDomain, gA);
      // ---------------------------------
    }

    // ------------------------------------------------------------------
    auto [min_p , max_p] = getMax(Mx_loc.x_result, localDomain, gA);
    {
      auto min_pTEMP = min_p;
      auto max_pTEMP = max_p;
      MPI_Allreduce(&min_pTEMP, &min_p, 1, MPI_DOUBLE, MPI_MIN, mx_comm_world);
      MPI_Allreduce(&max_pTEMP, &max_p, 1, MPI_DOUBLE, MPI_MAX, mx_comm_world);
    }
    if (simu.PID == 0) std::cout
    << "[Min_p, Max_p]:"<< min_p  << ", "<< max_p << std::endl;
    // ------------------------------------------------------------------

    #endif
    if (simu.PID == 0) simu.printInfo();

    timer.pressure.stop();
    // * -----------------------------------------------------------


    // ! ~~~~~~~~~~~~~~~~~~~~ send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood,
      t1_loc.p, localDomain, gA);
    // * ~~~~~~~~~~~~~~~~~~~ send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // !## 3. Implement the DFIB method,
    // ! the velocity (T1) be update to velocity (T3).
    timer.updateT1toT3.start();
    // ----------------------
    // UpdateForceAndVelocity(Dfib, simu, t1, T1, T3, globalDomain, gA);
    // ----------------------

    // ----------------------
    UpdateForceAndVelocity(Dfib_loc, simu, t1_loc, T1_loc, T3_loc, localDomain, gA); 
    // ----------------------
    timer.updateT1toT3.stop();
    // * -----------------------------------------------------------


    // ! Check the staedy state (L2 norm) -------------------------------
    timer.checkL2norm.start();
    // CheckSteadyStatebyMaxVal_omp(simu, ShareM, T0, T3, globalDomain, gA);
    // CheckSteadyStatebyMaxVal_mpi(mx_comm_world, simu, ShareM, T0, T3, localDomain, gA);
    // CheckSteadyStatebyL2norm(simu, T0_loc, T3_loc, localDomain, gA);
    timer.checkL2norm.stop();
    // *  ----------------------------------------------------------------

    // ! ## 4. copy the vel.
    timer.updateT3toT0.start();
    // T0.copy_omp(T3);
    T0_loc.copy_omp(T3_loc);
    timer.updateT3toT0.stop();
    // * --------------------


    // !## 5. Update vel and pressure at ghost cell.
    timer.BC.start();
    // UpdateAllVelocityOnBoundary( globalDomain, T0, t1, gA );
    // BC_staggered_copy( globalDomain, T0, T1, gA );
    UpdateAllVelocityOnBoundary(localDomain, T0_loc, t1_loc, gA);
    BC_staggered_copy(localDomain, T0_loc, T1_loc, gA);

    timer.BC.stop();
    // * -----------------------------------------

    // !## 6 Get both Cd and Cl.  -------------
    if (simu.DfibMethod != "OFF") {
      timer.getCdCl.start();
      getCD_CL_mpi(mx_comm_world, simu, localDomain, Dfib_loc, gA);
      timer.getCdCl.stop();
    }
    // * --------------------------------

    // !## 7 write plot3D formart .-------------
    timer.beginNew.stop();
    if (simu.get_writefile()) {
      // ! ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~
      for (int whichDim = 0; whichDim < 3 ; whichDim++) {
        mpi_iSR_double_x_Collect_to_Master(mx_comm_world,
        simu.PID, mpi_word_size, T3_loc.U(whichDim),
        localDomain, gA);
      }
      mpi_iSR_double_x_Collect_to_Master(
        mx_comm_world, simu.PID, mpi_word_size, t1_loc.p, localDomain, gA);

      MPI_Barrier(mx_comm_world);

      // * ~~~~~~~~~~~~~~~~~~~~~~~~~~  send & recv ~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (simu.PID == 0) {
        // WriteQfile(Dfib, simu, t1, T3, gA);
        WriteQfile(Dfib_loc, simu, t1_loc, T3_loc, gA);
      }
    }
    // * -----------------------------------------

    // !## 8. IO ------------------------------------------------
    MPI_Barrier(mx_comm_world);
    if (simu.PID == 0)
    cout
      <<  "[simu]{time, dt} ="  << simu.getSimuTime()
      <<  ", " << simu.dt       << endl
      <<  "[wall time] < "      << timer.beginNew.elapsedTime() - perLoopWTime
      <<  " OF "                << timer.beginNew.elapsedTime()
      <<  " >" << endl

    // ! =============== NEXT time step ===============
    <<     "===================================================" << endl
    <<     "LOOP :"  << simu.loop + 1
    <<     gA.show() << ", file :" << simu.get_file() << endl;
    // * ----------------------------------------------------------
    perLoopWTime = timer.beginNew.elapsedTime();
    timer.beginNew.start();

    recordingTime(timer, simu, ShareM);
  }
  timer.beginNew.stop();
  // ! ==================  main loop begin ==================


  MPI_Barrier(mx_comm_world);
  MPI_Finalize();
  return true;
}

