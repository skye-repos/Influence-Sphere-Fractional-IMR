#pragma once
#include "graph.h"
#include <cstddef>
#include <random>
#include <vector>

struct SimResult {
  std::vector<std::vector<int>> tracking;
  int sweep_count;
};

SimResult simulate_IMR(const Graph &g, std::mt19937 &rng, double θ, double F0,
                       int max_sweeps);

int new_state(const Graph &g, double θ, int node, const std::vector<int> &state);

double count_flips(const SimResult &result);
double count_ffrac(const SimResult &result);
