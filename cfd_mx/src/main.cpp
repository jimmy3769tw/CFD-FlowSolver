#include "../inc/backend/physical_variables.hpp"
#include "../inc/backend/simulation.hpp"
#include "../inc/grid/structured_grid.hpp"
#include "../inc/grid/uniform_structured_grid.hpp"
#include "../inc/run/run_cpu_mpi.hpp"
#include "../inc/run/run_cpu.hpp"
#include "../inc/run/run_cpu_mpi.hpp"
int main(int argc, char **argv) {  
  Simulation simu;
  simu.SetArgcArgv(argc, argv);
      UniformStructuredGrid grid(40, 40, 40);
  grid.setLen(1, 1, 1).Init();
  simu.SetReynoldsNumber(100.0)
      .tva.SetDt(0.001)
      .SetWritingDt(1)
      .SetTerminalTime(40);
  projection_method::CpuOpenMp projection_method(simu, grid);
  projection_method.solve();
  std::cout << "Finish !\n";
  return 0;
}