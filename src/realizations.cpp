#include "realizations.h"
#include "influence.h"
#include <iostream>
#include <omp.h>
#include <vector>

AvgSimResult avg_simulate_IMR(const int N, const double p,
                              const int num_realizations, const int &seed,
                              double θ_min, double θ_max, int Nθ, double F0_min,
                              double F0_max, int NF0) {

  double nr = static_cast<double>(num_realizations);

  std::vector<double> θ_list(Nθ);
  for (int i = 0; i < Nθ; ++i) {
    θ_list[i] = θ_min + i * (θ_max - θ_min) / Nθ;
  }

  std::vector<double> F0_list(NF0);
  for (int i = 0; i < NF0; ++i) {
    F0_list[i] = F0_min + i * (F0_max - F0_min) / NF0;
  }

  int num_threads = omp_get_max_threads();

  std::vector<std::vector<std::vector<double>>> t_flips(num_threads);
  std::vector<std::vector<std::vector<double>>> t_ffrac(num_threads);
  std::vector<std::vector<std::vector<double>>> t_ctime(num_threads);
  std::vector<std::vector<std::vector<double>>> t_brnch(num_threads);

  for (int t = 0; t < num_threads; ++t) {
    t_flips[t].assign(NF0, std::vector<double>(Nθ, 0.0));
    t_ffrac[t].assign(NF0, std::vector<double>(Nθ, 0.0));
    t_ctime[t].assign(NF0, std::vector<double>(Nθ, 0.0));
    t_brnch[t].assign(NF0, std::vector<double>(Nθ, 0.0));
  }

  std::cout << "N_realizations * Nθ * NF0 = " << num_realizations * Nθ * NF0
            << '\n';

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    std::mt19937 rng(seed + tid * 100000);

#pragma omp for schedule(dynamic, 1)
    for (int rel = 0; rel < num_realizations; ++rel) {
      Graph g = erdos_renyii(N, p, rng);
      std::cout << "Realization #" << (rel + 1) << "/" << num_realizations
                << '\n'
                << "N = " << N << ", edge count = " << g.ne() << '\n';

      for (int i = 0; i < Nθ; ++i) {
        double θ = θ_list[i];
        for (int j = 0; j < NF0; ++j) {
          double F0 = F0_list[j];
          SimResult result = simulate_IMR(g, rng, θ, F0, 100);

          auto n_flips = count_flips(result);
          auto n_ffrac = count_ffrac(result);
          auto Igraph = influence_graph(g, result.tracking[0], θ);

          t_flips[tid][j][i] += n_flips;
          t_ffrac[tid][j][i] += n_ffrac;
          t_ctime[tid][j][i] += static_cast<double>(result.sweep_count);
          t_brnch[tid][j][i] += Igraph.branching_factor;
          std::cout << "F0, θ: " << F0 << ", " << θ << '\n';
        }
      }
    }
  }
  std::vector<std::vector<double>> flips(NF0, std::vector<double>(Nθ, 0.0));
  std::vector<std::vector<double>> ffrac(NF0, std::vector<double>(Nθ, 0.0));
  std::vector<std::vector<double>> ctime(NF0, std::vector<double>(Nθ, 0.0));
  std::vector<std::vector<double>> brnch(NF0, std::vector<double>(Nθ, 0.0));

  for (int tid = 0; tid < num_threads; ++tid) {
    for (int fi = 0; fi < NF0; ++fi) {
      for (int ti = 0; ti < Nθ; ++ti) {
        flips[fi][ti] += t_flips[tid][fi][ti] / nr;
        ffrac[fi][ti] += t_ffrac[tid][fi][ti] / nr;
        ctime[fi][ti] += t_ctime[tid][fi][ti] / nr;
        brnch[fi][ti] += t_brnch[tid][fi][ti] / nr;
      }
    }
  }

  AvgSimResult avg_result;
  avg_result.flips_matrix = std::move(flips);
  avg_result.ffrac_matrix = std::move(ffrac);
  avg_result.ctime_matrix = std::move(ctime);
  avg_result.brnch_matrix = std::move(brnch);

  return avg_result;
}
