#pragma once
#include <string>
#include <vector>

inline auto virtualF_Int(ImmersedBoundary& Dfib, CalDomain& Lo, grid& gA,
                         const int whichDirection) {
  if (Dfib.ValSum.size() < 3) Dfib.ValSum.resize(3);

  const int ptr = gA.iceltotCal * whichDirection;
  double sum = 0.0;

#pragma omp parallel for firstprivate(Lo, ptr) reduction(+ : sum)
  for (auto i = Lo.x_start; i < Lo.x_end; ++i)
    for (auto j = Lo.y_start; j < Lo.y_end; ++j)
      for (auto k = Lo.z_start; k < Lo.z_end; ++k) {
        sum += Dfib.f[ptr + gA.icelCal(i, j, k)] *
               (gA.Dx[i] * gA.Dy[j] * gA.Dz[k]);
      }

  Dfib.ValSum[whichDirection] = sum;

  return 0;
}

std::pair<double, double> noDelay_IO_CD_CL(Simulation& simu,
                                           CalDomain& Lo,
                                           ImmersedBoundary& Dfib,
                                           grid& gA) {
  virtualF_Int(Dfib, Lo, gA, 0);
  virtualF_Int(Dfib, Lo, gA, 1);

  auto area_cD = gA.lz;
  auto area_cL = gA.lz;

  auto cD = -2.0 * Dfib.ValSum[0] / area_cD;
  auto cL = -2.0 * Dfib.ValSum[1] / area_cL;

  if (simu.PID == 0) {
    std::cout << "[cD, cL] : " << cD << ", " << cL << std::endl;

    std::ofstream file;

    std::string name = "Information/Time_cDcL";

    name += ".dat";

    file.open(name, std::ios::out | std::ios::app);

    std::string tab = " ";

    if (simu.loop == 1) {
      std::vector<std::string> variables;
      variables.push_back("simulation time");
      variables.push_back("C<sub>D</sub>");
      variables.push_back("C<sub>L</sub>");

      file << "TITLE     = \"\"\n"
           << "VARIABLES = \"" << variables.at(0) << "\",\"" << variables.at(1)
           << "\",\"" << variables.at(2) << "\"\n"
           << "ZONE T=\"" << simu.ZONE() << "\"";
    }

    file << "\n" << simu.getSimuTime() << tab << cD << tab << cL;

    file.close();
  }

  return std::make_pair(cD, cL);
}

std::pair<double, double> Delay_IO_CD_CL(Simulation& simu, CalDomain& Lo,
                                         ImmersedBoundary& Dfib,
                                         grid& gA) {
  std::string name = "Information/Time_cDcL";
  name += ".dat";

  virtualF_Int(Dfib, Lo, gA, 0);
  virtualF_Int(Dfib, Lo, gA, 1);

  auto area_cD = gA.lz;
  auto area_cL = gA.lz;

  auto dragCoefficient = -2.0 * Dfib.ValSum[0] / area_cD;
  auto liftCoefficient = -2.0 * Dfib.ValSum[1] / area_cL;

  simu.dragCoefficient[simu.loop % simu.delayIO] = dragCoefficient;
  simu.liftCoefficient[simu.loop % simu.delayIO] = liftCoefficient;

  bool IOfale = simu.loop % simu.delayIO;

  if (simu.PID != 0) return make_pair(dragCoefficient, liftCoefficient);

  if (simu.loop == 1) {
    std::ofstream file;
    file.open(name, std::ios::out | ios::app);
    std::vector<std::string> variables;
    variables.push_back("simulation time");
    variables.push_back("C<sub>D</sub>");
    variables.push_back("C<sub>L</sub>");

    file << "TITLE     = \"\"\n"
         << "VARIABLES = \"" << variables.at(0) << "\",\"" << variables.at(1)
         << "\",\"" << variables.at(2) << "\"\n"
         << "ZONE T=\"" << simu.ZONE() << "\"";
  }

  if (!(IOfale)) {
    std::ofstream file;
    file.open(name, std::ios::out | std::ios::app);
    std::string tab = " ";
    for (int i = 1; i < simu.delayIO + 1; ++i)
      file << "\n"
           << simu.getSimuTime() << tab
           << simu.dragCoefficient[simu.loop % simu.delayIO] << tab
           << simu.liftCoefficient[simu.loop % simu.delayIO];

    file.close();
  }

  return make_pair(dragCoefficient, liftCoefficient);
}
