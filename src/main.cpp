#include "io.h"
#include "realizations.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

struct Config {
  int N = 20000;
  double kavg = 4.5;
  int nrel = 25;
  double tmin = 0.00;
  double tmax = 1.00;
  int Nt = 1000;
  double F0min = 0.10;
  double F0max = 0.55;
  int NF0 = 10;
  int seed = 12345;
  int max_sweeps = 100;
};

void print_usage(const char *prog) {
  std::cerr << "Usage: " << prog << " [options]\n"
            << "  -N <int>           Number of nodes (default 1000)\n"
            << "  -kavg <float>      Average degree (default 4.5)\n"
            << "  -nrel <int>        Number of realisations (default 25)\n"
            << "  -tmin <float>      Min theta (default 0.40)\n"
            << "  -tmax <float>      Max theta (default 0.80)\n"
            << "  -Nt <int>          # of Theta steps (default 100)\n"
            << "  -F0min <float>     Min F0 (default 0.05)\n"
            << "  -F0max <float>     Max F0 (default 0.50)\n"
            << "  -NF0 <int>         # of F0 steps (default 5)\n"
            << "  -seed <int>        Random seed (default 12345)\n"
            << "  -max-sweeps <int>  Max sweeps / branching-factor steps "
               "(default 100)\n"
            << "  -h                 Show this help\n";
};

Config parse_args(int argc, char *argv[]) {
  Config cfg;

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-N") == 0)
      cfg.N = std::atoi(argv[++i]);
    else if (std::strcmp(argv[i], "-kavg") == 0)
      cfg.kavg = std::atof(argv[++i]);
    else if (std::strcmp(argv[i], "-nrel") == 0)
      cfg.nrel = std::atoi(argv[++i]);
    else if (std::strcmp(argv[i], "-tmin") == 0)
      cfg.tmin = std::atof(argv[++i]);
    else if (std::strcmp(argv[i], "-tmax") == 0)
      cfg.tmax = std::atof(argv[++i]);
    else if (std::strcmp(argv[i], "-Nt") == 0)
      cfg.Nt = std::atoi(argv[++i]);
    else if (std::strcmp(argv[i], "-F0min") == 0)
      cfg.F0min = std::atof(argv[++i]);
    else if (std::strcmp(argv[i], "-F0max") == 0)
      cfg.F0max = std::atof(argv[++i]);
    else if (std::strcmp(argv[i], "-NF0") == 0)
      cfg.NF0 = std::atoi(argv[++i]);
    else if (std::strcmp(argv[i], "-seed") == 0)
      cfg.seed = std::atoi(argv[++i]);
    else if (std::strcmp(argv[i], "-max-sweeps") == 0)
      cfg.max_sweeps = std::atoi(argv[++i]);
    else if (std::strcmp(argv[i], "-h") == 0) {
      print_usage(argv[0]);
      std::exit(0);
    } else {
      std::cerr << "Unknown flag: " << argv[i] << '\n';
      print_usage(argv[0]);
      std::exit(1);
    }
  }

  return cfg;
}

int main(int argc, char *argv[]) {
  Config cfg = parse_args(argc, argv);

  const int N = std::move(cfg.N);
  const double p = cfg.kavg / (N - 1);
  const int num_realizations = std::move(cfg.nrel);

  const double θ_min = std::move(cfg.tmin);
  const double θ_max = std::move(cfg.tmax);
  const int Nθ = std::move(cfg.Nt);
  const double F0_min = std::move(cfg.F0min);
  const double F0_max = std::move(cfg.F0max);
  const int NF0 = std::move(cfg.NF0);

  AvgSimResult results =
      avg_simulate_IMR(N, p, num_realizations, cfg.seed, θ_min, θ_max, Nθ,
                       F0_min, F0_max, NF0, cfg.max_sweeps);

  std::vector<std::string> θ_labels(Nθ);
  for (int i = 0; i < Nθ; ++i) {
    double val = θ_min + i * (θ_max - θ_min) / Nθ;
    θ_labels[i] = std::to_string(val);
  }

  std::vector<std::string> F0_labels(NF0);
  for (int i = 0; i < NF0; ++i) {
    double val = F0_min + i * (F0_max - F0_min) / NF0;
    F0_labels[i] = std::to_string(val);
  }

  auto brnch_path = "results/N = " + std::to_string(N) + "/CSV/brnch.csv";
  write_csv(brnch_path, θ_labels, F0_labels, results.brnch_matrix);

  auto brnch_step_dir =
      "results/N = " + std::to_string(N) + "/CSV/brnch_step/";
  for (size_t s = 0; s < results.brnch_step_matrix.size(); ++s) {
    auto brnch_step_path =
        brnch_step_dir + "brnch_" + std::to_string(s) + ".csv";
    write_csv(brnch_step_path, θ_labels, F0_labels,
              results.brnch_step_matrix[s]);
  }
}
