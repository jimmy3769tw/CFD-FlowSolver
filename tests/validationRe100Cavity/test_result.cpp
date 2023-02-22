#include <gtest/gtest.h>

#include "../../cfd_mx/inc/backend/simulation.hpp"
#include "../../cfd_mx/inc/grid/structured_grid.hpp"
#include "../../cfd_mx/inc/grid/uniform_structured_grid.hpp"
#include "../../cfd_mx/inc/io/plot3d/read_qfile.hpp"
#include "../../cfd_mx/inc/run/run_cpu.hpp"
#include "../validation_tool/l2norm_validation.hpp"

TEST(re100, test1) {
  Simulation simu;
  UniformStructuredGrid grid(40, 40, 40);
  grid.setLen(1, 1, 1).Init();


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
