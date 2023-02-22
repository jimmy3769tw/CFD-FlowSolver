

auto unit_SML(vector<double> &a, const double lSml, int nS, int off_nS,
              double l, int n) {
  cout << "\n-------------unit_SML -------------\n";

  // * | goast Cell | Part[0](L) | Part[1].lx | Part[2].lx| Part[3].lx |
  // Part[4].lx|goast Cell|
  // * | 0 1        | Part[0].nx | dySml | nx-2, nx-1|
  auto nCal = n - 2 * no_ghost_cell;

  const double dS = lSml / nS, d = (l - lSml) / (nCal - nS);

  const double SBgn = 0.5 * (l - lSml), SEnd = 0.5 * (l + lSml);

  a.at(no_ghost_cell) = 0.0;

  if (nS != 0) {
    // * Calculation the grid
    if ((n - nS) % 2 != 0) throw std::invalid_argument("(n-nSml)%2 != 0");

    int nlEnd[2];
    int nsEnd[1];

    nlEnd[0] = (nCal - nS) / 2;
    nsEnd[0] = nlEnd[0] + nS;
    nlEnd[1] = nsEnd[0] + nlEnd[0];

    nlEnd[0] += off_nS;
    nsEnd[0] += off_nS;

    cout << "nlEnd[0] " << nlEnd[0] << "\n";
    cout << "nsEnd[0] " << nsEnd[0] << "\n";
    cout << "nlEnd[1] " << nlEnd[1] << "\n";

    if (nlEnd[0] <= 0) throw std::invalid_argument("nL[0]");

    for (size_t i = 0; i < nlEnd[0]; i++)
      a.at(i + no_ghost_cell + 1) = a.at(i + no_ghost_cell) + d;

    for (size_t i = nlEnd[0]; i < nsEnd[0]; ++i)
      a.at(i + no_ghost_cell + 1) = a.at(i + no_ghost_cell) + dS;

    for (size_t i = nsEnd[0]; i < nlEnd[1]; i++)
      a.at(i + no_ghost_cell + 1) = a.at(i + no_ghost_cell) + d;

  } else {
    for (size_t i = no_ghost_cell; i < n - no_ghost_cell; i++)
      a.at(i + 1) = a.at(i) + l / nCal;
  }

  init_gC_dx();

  return true;
}


bool init_nonuniform_SML() {
  // * | goast Cell | Part[0](L) | Part[1].lx | Part[2].lx| Part[3].lx |
  // Part[4].lx|goast Cell|
  // * | 0 1        | Part[0].nx | dySml| nx-2, nx-1|
  // *Setting
  const double lxS = 0.3 * lx, lyS = 0.3 * ly, lzS = 0.5 * lz;

  const int nxS = 100, nyS = 100, nzS = 0;

  const int off_nS_x = -4, off_nS_y = 0, off_nS_z = 0;

  unit_SML(x_pos, lxS, nxS, off_nS_x, lx, nx);
  unit_SML(y_pos, lyS, nyS, off_nS_y, ly, ny);
  unit_SML(z_pos, lzS, nzS, off_nS_z, lz, nz);

  return true;
}

void init_uniform() {
  auto dx = lx / cal_nx;
  auto dy = ly / cal_ny;
  auto dz = lz / cal_nz;
  double Xi = 0.0, Yi = 0.0, Zi = 0.0;

  x_pos.at(0) = Xi - 2.0 * dx;
  for (size_t i = 1; i < nx + 1; ++i) {
    x_pos.at(i) = x_pos.at(i - 1) + dx;
  }

  y_pos[0] = Yi - 2.0 * dy;
  for (size_t j = 1; j < ny + 1; ++j) {
    y_pos.at(j) = y_pos.at(j - 1) + dy;
  }

  z_pos[0] = Zi - 2.0 * dz;
  for (size_t k = 1; k < nz + 1; ++k) {
    z_pos.at(k) = z_pos.at(k - 1) + dz;
  }

  fill(dx.begin(), dx.end(), dx);
  fill(dy.begin(), dy.end(), dy);
  fill(dz.begin(), dz.end(), dz);

  fill(staggered_dx.begin(), staggered_dx.end(), dx);
  fill(staggered_dy.begin(), staggered_dy.end(), dy);
  fill(staggered_dz.begin(), staggered_dz.end(), dz);
}

