#include "backend/physical_variables.hpp"
#include "backend/simulation.hpp"
#include "grid/structured_grid.hpp"
#include "grid/parital_uniform_structured_gird.hpp"
#include "run/run_cpu.hpp"

int main(int argc, char **argv) {  
  Simulation simu;
  simu.SetArgcArgv(argc, argv);

  PartialUniformStructuredGrid grid(256, 80, 30);
  grid.setLen(17, 12, 3).Init();
  grid.bc_selector.FlowPassing();

  // grid.cal_eta_factory.SetCylinderZ()
          // .SetCylinderCenter()
          // .SetRadius(0.5)
          // .SetSubGrid(100, 100);


  simu.SetReynoldsNumber(40.0)
      .SetVelocityInitialCondition(1, 0, 0)
      .tva.SetDt(0.001)
      .SetWritingDt(1)
      .SetTerminalTime(40);

  projection_method::CpuOpenMp projection_method(simu, grid);
  projection_method.solve();
  std::cout << "Finish !\n";
  return 0;
}

