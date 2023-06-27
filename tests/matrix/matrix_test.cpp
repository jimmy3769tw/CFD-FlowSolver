#include <gtest/gtest.h>

#include "../../cfd_mx/inc/matrix/csr_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/ell_sparse_mat.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab.hpp"
#include "../../cfd_mx/inc/matrix/solver/bicgstab_restart.hpp"
#include "../validation_tool/l2norm_validation.hpp"
#include "GetRandomVector.hpp"
#include "io_2d_plot3d.hpp"
#include "poisson.hpp"

namespace {  
class CsrMatTest : public ::testing::Test {
};

TEST(CsrMatTest, resize_row_col) {
    mat::CsrMat<double> matA(5, 4);
    EXPECT_EQ(matA.row(), 5);
    EXPECT_EQ(matA.col(), 4);
    matA.resize(4, 3);
    EXPECT_EQ(matA.row(), 4);
    EXPECT_EQ(matA.col(), 3);
}

TEST(CsrMatTest, get_csr) {
  mat::CsrMat<double> matA;
  HeatConduction heat_conduction(10);
  matA.Set(heat_conduction.GetCsr());
  EXPECT_EQ(heat_conduction.GetCsr(), matA.GetCsr());
}

} // namespace