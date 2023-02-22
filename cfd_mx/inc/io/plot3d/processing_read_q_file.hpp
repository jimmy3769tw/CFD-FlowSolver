
auto ProcessReadQfileint(int Nblock, float mach, float alpha, float reyn,
                         float time, std::vector<std::vector<float>> data) {
  std::vector<std::vector<double>> data_new(
      nodata, std::vector<double>(grid.no_grid, 0.0));

  for (size_t k = gc; k < cal_nz + gc; ++k)
    for (size_t j = gc; j < cal_ny + gc; ++j)
      for (size_t i = gc; i < cal_nx + gc; ++i) {
        data_new[0][grid.icel(i, j, k)] = data[0][grid.icelCal(i, j, k)];
        data_new[4][grid.icel(i, j, k)] = data[4][grid.icelCal(i, j, k)];
      }

  for (size_t k = gc; k < cal_nz + gc; ++k)
    for (size_t j = gc; j < cal_ny + gc; ++j)
      for (size_t i = gc; i < cal_nx + gc; ++i) {
        if (nxCal + gc - 1 == i)
          data_new[1][grid.icel(i, j, k)] = data[1][grid.icelCal(i, j, k)];
        else
          data_new[1][grid.icel(i, j, k)] =
              (data[1][grid.icelCal(i, j, k)] +
               data[1][grid.icelCal(i + 1, j, k)]) /
              2.0;

        if (nyCal + gc - 1 == j)
          data_new[2][grid.icel(i, j, k)] = data[2][grid.icelCal(i, j, k)];
        else
          data_new[2][grid.icel(i, j, k)] =
              (data[2][grid.icelCal(i, j, k)] +
               data[2][grid.icelCal(i, j + 1, k)]) /
              2.0;

        if (nzCal + gc - 1 == k)
          data_new[3][grid.icel(i, j, k)] = data[3][grid.icelCal(i, j, k)];
        else
          data_new[3][grid.icel(i, j, k)] =
              (data[3][grid.icelCal(i, j, k)] +
               data[3][grid.icelCal(i, j, k + 1)]) /
              2.0;
      }

  return std::make_tuple(Nblock, mach, alpha, reyn, time, data_new[0],
                         data_new[1], data_new[2], data_new[3], data_new[4]);
}
