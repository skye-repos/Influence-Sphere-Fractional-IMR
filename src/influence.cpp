#include "influence.h"
#include "simulate.h"
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

// Fast network branching factor.  Only depth-1 pure flips are counted, so the
// full influence-sphere BFS is unnecessary.  For each node we precompute its
// degree-weighted opinion sum; flipping a node only shifts the contribution of
// its own degree to its neighbours, applied in O(deg) per flip.  The reference
// applies each flipped neighbour before evaluating the next one, so the same
// sequential order is reproduced here.
double branching_factor(const Graph &g, const std::vector<int> &state,
                        double θ) {
  const int N = g.nv();

  std::vector<double> base_m1(N, 0.0);
  std::vector<double> tot(N, 0.0);
  for (int u = 0; u < N; ++u) {
    for (int nbr : g.neighbors(u)) {
      double w = g.deg[nbr];
      tot[u] += w;
      if (state[nbr] == 1)
        base_m1[u] += w;
    }
  }

  std::vector<int> orig_op(N, 0);
  for (int u = 0; u < N; ++u) {
    orig_op[u] = (tot[u] > 0.0 && base_m1[u] / tot[u] >= θ) ? 1 : 0;
  }

  std::vector<double> cur_m1 = base_m1;
  std::vector<int> touched;
  std::vector<char> mark(N, 0);
  touched.reserve(N);

  auto add_delta = [&](int x, double d) {
    if (!mark[x]) {
      mark[x] = 1;
      touched.push_back(x);
    }
    cur_m1[x] += d;
  };

  double sum = 0.0;
  for (int v = 0; v < N; ++v) {
    const double dv = static_cast<double>(g.deg[v]);

    // Flip the seed: neighbours of v see v's opinion inverted.
    double seed_delta = (state[v] == 0) ? dv : -dv;
    for (int x : g.neighbors(v))
      add_delta(x, seed_delta);

    int branch_v = 0;
    for (int u : g.neighbors(v)) {
      if (orig_op[u] != state[u])
        continue; // u is not stable in the original state
      int flip_op = (tot[u] > 0.0 && cur_m1[u] / tot[u] >= θ) ? 1 : 0;
      if (flip_op != state[u]) {
        ++branch_v;
        // Apply u's flip for the remaining neighbours of the seed.
        double du = static_cast<double>(g.deg[u]);
        double delta = (state[u] == 0) ? du : -du;
        for (int x : g.neighbors(u))
          add_delta(x, delta);
      }
    }
    sum += branch_v;

    for (int x : touched) {
      cur_m1[x] = base_m1[x];
      mark[x] = 0;
    }
    touched.clear();
  }

  return sum / static_cast<double>(N);
}

OpinionInstability opinion_instability(const std::vector<double> &I_total,
                                       std::vector<int> &state) {
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
