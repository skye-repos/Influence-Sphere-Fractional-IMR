#include "simulate.h"
#include "graph.h"
#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <random>
#include <vector>

// Logic to calculate the new state
int new_state(const Graph &g, double θ, int node,
              const std::vector<int> &state) {
  double m0 = 0;
  double m1 = 0;

  for (const int &nbr : g.neighbors(node)) {
    (state[nbr] == 1) ? m1 += g.deg[nbr] : m0 += g.deg[nbr];
  }

  double frac = (m0 + m1 > 0) ? m1 / (m0 + m1) : 0.0;

  int new_node_state = (frac >= θ) ? 1 : 0;

  return new_node_state;
}

// Actual Simulations
SimResult simulate_IMR(const Graph &g, std::mt19937 &rng, double θ, double F0,
                       int max_sweeps) {
  const int N = g.nv();

  std::vector<int> state(N);
  std::bernoulli_distribution coin(F0);

  for (int i = 0; i < N; ++i) {
    state[i] = coin(rng) ? 0 : 1;
  }

  SimResult results;
  results.tracking.push_back(state);
  results.sweep_count = 0;

  std::vector<int> order(N);
  std::iota(order.begin(), order.end(), 0);

  for (int sweep = 0; sweep < max_sweeps; ++sweep) {
    int changes = 0;
    std::vector<int> current = results.tracking.back();

    std::shuffle(order.begin(), order.end(), rng);

    for (int node : order) {
      int new_node_state = new_state(g, θ, node, current);

      if (new_node_state != current.at(node)) {
        current[node] = new_node_state;
        ++changes;
      }
    }

    results.tracking.push_back(current);
    ++results.sweep_count;

    state = std::move(current);

    if (changes == 0)
      break;
  }

  return results;
}

// Metrics
double count_flips(const SimResult &result) {
  const int l = result.sweep_count + 1;
  const int N = result.tracking[0].size();
  std::vector<int> flippers(N);
  std::vector<int> floppers(N);

  for (int i = 0; i < l - 1; ++i) {
    std::vector<int> prev = result.tracking[i];
    std::vector<int> next = result.tracking[i + 1];
    for (int m = 0; m < N; ++m) {
      int val = next[m] - prev[m];
      flippers[m] += std::abs(val);
    }
  }

  for (int i = 0; i < N; ++i) {
    floppers[i] = (flippers[i] > 1) ? 1 : 0;
  }

  double n_flips =
      std::reduce(floppers.begin(), floppers.end()) / static_cast<double>(N);
  return n_flips;
}

double count_ffrac(const SimResult &result) {
  int zeros = 0;
  for (int s : result.tracking.back()) {
    if (s == 0) {
      ++zeros;
    };
  };

  return static_cast<double>(zeros) / result.tracking[0].size();
}
