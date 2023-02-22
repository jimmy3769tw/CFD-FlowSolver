#include <gtest/gtest.h>

#include "../../cfd_mx/inc/backend/simulation.hpp"
#include "../../cfd_mx/inc/grid/structured_grid.hpp"
#include "../../cfd_mx/inc/grid/uniform_structured_grid.hpp"
#include "../../cfd_mx/inc/io/plot3d/read_qfile.hpp"
#include "../../cfd_mx/inc/run/run_cpu_mpi.hpp"
#include "../validation_tool/l2norm_validation.hpp"

TEST(re100, test1) {
  Simulation simu;
  UniformStructuredGrid grid(80, 80, 80);
  grid.setLen(1, 1, 1).Init();
  simu.SetReynoldsNumber(1000.0)
      .tva.SetDt(0.001)
      .SetWritingDt(1)
      .SetTerminalTime(2);
  projection_method::CpuOpenMpMpi projection_method(simu, grid);
  projection_method.solve();
  auto validation = get<5>(ReadQfile(grid, "P3D40.q"));
  auto result = get<5>(ReadQfile(grid, "mx_out/P3D02.q"));

  for (int i = 0; i < result.size(); i++) {
    EXPECT_TRUE(L2Norm(validation[i], result[i]) < 1.0e-15);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
