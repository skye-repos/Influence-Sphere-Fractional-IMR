#pragma once
#include "influence.h"
#include "simulate.h"
#include <vector>

struct AvgSimResult {
  std::vector<std::vector<double>> flips_matrix;
  std::vector<std::vector<double>> ffrac_matrix;
  std::vector<std::vector<double>> ctime_matrix;
  std::vector<std::vector<double>> exin0_matrix;
  std::vector<std::vector<double>> exin1_matrix;
  std::vector<std::vector<double>> brnch_matrix;
  std::vector<std::vector<double>> outdg_matrix;
};

AvgSimResult avg_simulate_IMR(const int N, const double p,
                              const int num_realizations, const int &seed,
                              double θ_min, double θ_max, int Nθ, double F0_min,
                              double F0_max, int NF0);
