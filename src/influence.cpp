#include "influence.h"
#include "simulate.h"
#include <numeric>
#include <queue>
#include <vector>

InfluenceResult influence_sphere(const Graph &g, const std::vector<int> &state,
                                 int node, double θ) {
  const int N = g.nv();
  InfluenceResult result;
  result.dist.assign(N, -1);
  result.branching = 0;

  std::vector<int> flipped = state;
  flipped[node] = 1 - flipped[node];

  std::queue<int> Q;
  Q.push(node);
  result.dist[node] = 0;
  result.affected.push_back(node);

  while (!Q.empty()) {
    int v = Q.front();
    Q.pop();

    for (int nbr : g.neighbors(v)) {
      if (result.dist[nbr] != -1)
        continue;

      int orig_op = new_state(g, θ, nbr, state);
      int flip_op = new_state(g, θ, nbr, flipped);

      if (flip_op != state[nbr] && orig_op == state[nbr]) {
        flipped[nbr] = flip_op;
        result.dist[nbr] = result.dist[v] + 1;
        result.affected.push_back(nbr);
        if (v == node) {
          ++result.branching;
        }
        Q.push(nbr);
      }
    }
  }

  return result;
}

std::tuple<std::vector<double>, std::vector<double>>
total_instability(const Graph &g, const std::vector<int> &state, double θ) {
  const int N = g.nv();
  std::vector<double> I(N, 0.0);
  std::vector<double> b(N, 0.0);

  for (int node = 0; node < N; ++node) {
    auto sphere = influence_sphere(g, state, node, θ);
    b[node] = sphere.branching;

    for (int a : sphere.affected) {
      if (a == node)
        continue;
      I[a] += 1.0 / sphere.dist[a];
    }
  }

  double I_max = 1.0 / (N - 1.0);
  for (double &val : I)
    val *= I_max;

  return {I, b};
}

OpinionInstability opinion_instability(const std::vector<double> &I_total,
                                       const std::vector<int> &state) {
  int N = I_total.size();
  double n = static_cast<double>(N);
  OpinionInstability I01;
  I01.I0.assign(N, 0.0);
  I01.I1.assign(N, 0.0);

  double F0 = 0.0, F1 = 0.0;

  for (int i = 0; i < N; ++i) {
    if (state[i] == 0) {
      F0 += 1.0 / n;
    } else {
      F1 += 1.0 / n;
    }
  }

  for (int i = 0; i < N; ++i) {
    if (state[i] == 0 && F0 != 0.0) {
      I01.I0[i] = I_total[i] / F0;
    } else if (state[i] == 1 && F1 != 0.0) {
      I01.I1[i] = I_total[i] / F1;
    }
  }

  return I01;
}

double expectation(const std::vector<double> &I_vec) {
  int N = I_vec.size();
  double I_sum = std::reduce(I_vec.begin(), I_vec.end());

  return I_sum / static_cast<double>(N);
}

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
