    #elif defined (P_SOLVER_AMGCL_EIGEN)

      createBMatrix_Eigen(T1, Mx, simu, Lo, gridA);

      createBMatrix_seq(T1, Mx, simu, Lo, gridA);


      std::tie(simu.iters, simu.error) = solve(Mx.matB_Eigen, Mx.x_Eigen);

      Pressure_transform_x_Eigen(t1, Mx, Lo, gridA);

      simu.iters = iters;
