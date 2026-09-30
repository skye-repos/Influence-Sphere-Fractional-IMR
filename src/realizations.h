#pragma once
#include "influence.h"
#include "simulate.h"

struct AvgSimResult {
  std::vector<std::vector<double>> brnch_matrix;
  // [sweep][F0][θ]: network branching factor at each sweep, averaged over runs
  std::vector<std::vector<std::vector<double>>> brnch_step_matrix;
};

AvgSimResult avg_simulate_IMR(const int N, const double p,
                              const int num_realizations, const int &seed,
                              double θ_min, double θ_max, int Nθ, double F0_min,
                              double F0_max, int NF0, int max_sweeps);
