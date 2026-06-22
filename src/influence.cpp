#include "influence.h"
#include "simulate.h"
#include <vector>

DiGraph influence_graph(const Graph &g, const std::vector<int> &state,
                        const double &θ) {
  const int N = g.nv();

  DiGraph result(N);

  for (int node = 0; node < N; ++node) {
    std::vector<int> flipped = state;
    flipped[node] = 1 - flipped[node];

    for (int nbr : g.neighbors(node)) {
      int orig_op = new_state(g, θ, nbr, state);
      int flip_op = new_state(g, θ, nbr, flipped);

      if (flip_op != state[nbr] && orig_op == state[nbr]) {
        result.add_edge(node, nbr);
      }
    }
  }

  return result;
}
