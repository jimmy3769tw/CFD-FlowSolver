#include <gtest/gtest.h>

#include "../../cfd_mx/inc/backend/simulation.hpp"
#include "../../cfd_mx/inc/grid/structured_grid.hpp"
#include "../../cfd_mx/inc/grid/uniform_structured_grid.hpp"
#include "../../cfd_mx/inc/io/plot3d/read_qfile.hpp"
#include "../validation_tool/l2norm_validation.hpp"

TEST(re100, test1) {
  Simulation simu;
  UniformStructuredGrid grid(80, 80, 80);
  grid.setLen(1, 1, 1).Init();

  auto validation = get<5>(ReadQfile(grid, "P3D40.q"));
  auto result = get<5>(ReadQfile(grid, "mx_out/P3D01.q"));

  for (int i = 0; i < result.size(); i++) {
    EXPECT_TRUE(L2Norm(validation[i], result[i]) < 0.1e-30);
  }
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
