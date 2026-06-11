#include "io.h"
#include "realizations.h"
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <utility>

struct Config {
  int N = 1000;
  double kavg = 4.5;
  int nrel = 5;
  double tmin = 0.40;
  double tmax = 0.80;
  int Nt = 100;
  double F0min = 0.05;
  double F0max = 0.50;
  int NF0 = 5;
  int seed = 12345;
};

void print_usage(const char *prog) {
  std::cerr << "Usage: " << prog << " [options]\n"
            << "  -N <int>           Number of nodes (default 1000)\n"
            << "  -kavg <float>      Average degree (default 4.5)\n"
            << "  -nrel <int>        Number of realisations (default 5)\n"
            << "  -tmin <float>      Min theta (default 0.40)\n"
            << "  -tmax <float>      Max theta (default 0.80)\n"
            << "  -Nt <int>          # of Theta steps (default 100)\n"
            << "  -F0min <float>     Min F0 (default 0.05)\n"
            << "  -F0max <float>     Max F0 (default 0.50)\n"
            << "  -NF0 <int>         # of F0 steps (default 5)\n"
            << "  -seed <int>        Random seed (default 12345)\n"
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

  AvgSimResult results = avg_simulate_IMR(N, p, num_realizations, cfg.seed, θ_min,
                                          θ_max, Nθ, F0_min, F0_max, NF0);

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

  auto flips_path = "results/N = " + std::to_string(N) + "/CSV/flips.csv";
  auto ffrac_path = "results/N = " + std::to_string(N) + "/CSV/ffrac.csv";
  auto ctime_path = "results/N = " + std::to_string(N) + "/CSV/ctime.csv";
  auto exin0_path = "results/N = " + std::to_string(N) + "/CSV/exin0.csv";
  auto exin1_path = "results/N = " + std::to_string(N) + "/CSV/exin1.csv";

  write_csv(flips_path, θ_labels, F0_labels, results.flips_matrix);
  write_csv(ffrac_path, θ_labels, F0_labels, results.ffrac_matrix);
  write_csv(ctime_path, θ_labels, F0_labels, results.ctime_matrix);
  write_csv(exin0_path, θ_labels, F0_labels, results.exin0_matrix);
  write_csv(exin1_path, θ_labels, F0_labels, results.exin1_matrix);
}
