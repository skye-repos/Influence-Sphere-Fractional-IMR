#pragma once
#include "graph.h"
#include "simulate.h"
#include <random>
#include <tuple>
#include <vector>

struct InfluenceDiGraph {
  double branching_factor;
  DiGraph influence;
};

InfluenceDiGraph influence_graph(const Graph &g, const std::vector<int> &state,
                                 double θ);
