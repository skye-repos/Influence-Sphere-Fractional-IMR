#pragma once
#include <cstddef>
#include <random>
#include <vector>

struct Graph {
  std::vector<std::vector<int>> adj;
  std::vector<int> deg;

  explicit Graph(int N) : adj(N), deg(N, 0) {}

  int nv() const { return static_cast<int>(adj.size()); }
  int ne() const;

  void add_edge(int u, int v);
  void rem_edge(int u, int v);
  bool has_edge(int u, int v) const;

  const std::vector<int> &neighbors(int u) const { return adj[u]; }
};

Graph erdos_renyii(int N, double p, std::mt19937 &rng);
