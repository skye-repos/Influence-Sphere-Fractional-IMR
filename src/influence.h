#pragma once
#include "graph.h"
#include "simulate.h"
#include <random>
#include <tuple>
#include <vector>

struct InfluenceResult {
  std::vector<int> dist;
  std::vector<int> affected;
  int branching;
};

InfluenceResult influence_sphere(const Graph &g, const std::vector<int> &state,
                                 int node, double θ);

std::tuple<std::vector<double>, std::vector<double>>
total_instability(const Graph &g, const std::vector<int> &state, double θ);

struct OpinionInstability {
  std::vector<double> I0;
  std::vector<double> I1;
};

OpinionInstability opinion_instability(const std::vector<double> &I_total,
                                       std::vector<int> &state);

double expectation(const std::vector<double> &I_vec);
