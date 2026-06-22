#pragma once
#include "graph.h"
#include "simulate.h"
#include <random>
#include <tuple>
#include <vector>

DiGraph influence_graph(const Graph &g, const std::vector<int> &state,
                        const double &θ);
