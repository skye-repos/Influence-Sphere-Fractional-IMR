#include "graph.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>

int Graph::ne() const {
  int total = 0;
  for (int d : deg)
    total += d;
  return total / 2;
}

int DiGraph::ne() const {
  int total = 0;
  for (int d_out : out_deg) {
    total += d_out;
  }
  for (int d_in : in_deg) {
    total += d_in;
  }
  return total;
}

void Graph::add_edge(int u, int v) {
  adj[u].push_back(v);
  adj[v].push_back(u);
  deg[u]++;
  deg[v]++;
}

void DiGraph::add_edge(int u, int v) {
  adj[u].push_back(v);
  out_deg[u]++;
  in_deg[v]++;
}

void Graph::rem_edge(int u, int v) {
  auto it = std::find(adj[u].begin(), adj[u].end(), v);

  if (it != adj[u].end()) {
    adj[u].erase(it);
    deg[u]--;
  }

  it = std::find(adj[v].begin(), adj[v].end(), u);

  if (it != adj[v].end()) {
    adj[v].erase(it);
    deg[v]--;
  }
}

void DiGraph::rem_edge(int u, int v) {
  auto it = std::find(adj[u].begin(), adj[u].end(), v);

  if (it != adj[u].end()) {
    adj[u].erase(it);
    out_deg[u]--;
    in_deg[v]--;
  }
}

bool Graph::has_edge(int u, int v) const {
  return std::find(adj[u].begin(), adj[u].end(), v) != adj[u].end();
}

bool DiGraph::has_edge(int u, int v) const {
  return std::find(adj[u].begin(), adj[u].end(), v) != adj[u].end();
}

Graph erdos_renyii(int N, double p, std::mt19937 &rng) {
  Graph g(N);
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  std::uniform_int_distribution<> dist_range(0, (N - 1));
  int n_edges = std::ceil(p * N * (N - 1) / 2);

  if (p <= 0.1) {
    int added = 0;
    while (added < n_edges) {
      int u = dist_range(rng);
      int v = dist_range(rng);
      if (u != v && !g.has_edge(u, v)) {
        g.add_edge(u, v);
        ++added;
      }
    }
  } else {
    for (int i = 0; i < N; ++i) {
      for (int j = i + 1; j < N; ++j) {
        if (dist(rng) < p)
          g.add_edge(i, j);
      }
    }
  }

  return g;
}
