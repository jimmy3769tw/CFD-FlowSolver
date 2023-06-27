#include <gtest/gtest.h>
#include "backend/simulation.hpp"
#include "grid/structured_grid.hpp"
#include "grid/uniform_structured_grid.hpp"
#include "io/plot3d/read_qfile.hpp"
#include "run/run_cpu_mpi.hpp"
#include "../validation_tool/l2norm_validation.hpp"

TEST(re100, test1) {
  Simulation simu;
  UniformStructuredGrid grid(40, 40, 40);
  grid.setLen(1, 1, 1).Init();
  simu.SetReynoldsNumber(100.0)
      .tva.SetDt(0.001)
      .SetWritingDt(1)
      .SetTerminalTime(19);
    
  grid.bc_selector.CavityFlow();
    
  
  projection_method::CpuOpenMpMpi projection_method(simu, grid);
  projection_method.solve();
  auto validation = get<5>(ReadQfile(grid, "P3D19.q"));
  auto result = get<5>(ReadQfile(grid, "mx_out/P3D19.q"));

  for (int i = 0; i < result.size(); i++) {
    auto l2norm = L2Norm(validation[i], result[i]);
    std::cout << "l2norm: " << l2norm << std::endl;
    EXPECT_TRUE(l2norm < 1.0e-6);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
