#include "influence.h"
#include "simulate.h"
#include <numeric>
#include <queue>
#include <vector>

InfluenceDiGraph influence_graph(const Graph &g, const std::vector<int> &state,
                                 double θ) {
  const int N = g.nv();

  InfluenceDiGraph result{0.0, DiGraph(N)};

  for (int node = 0; node < N; ++node) {
    std::vector<int> flipped = state;
    flipped[node] = 1 - flipped[node];

    for (int nbr : g.neighbors(node)) {
      int orig_op = new_state(g, θ, nbr, state);
      int flip_op = new_state(g, θ, nbr, flipped);
      double bf = 0.0;

      if (flip_op != state[nbr] && orig_op == state[nbr]) {
        result.influence.add_edge(node, nbr);
        ++bf;
      }
      result.branching_factor += bf / static_cast<double>(N);
    }
  }

  return result;
}
