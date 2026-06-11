#include "io.h"
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>

void write_csv(const std::string &path,
               const std::vector<std::string> &col_labels,
               const std::vector<std::string> &row_labels,
               const std::vector<std::vector<double>> &results) {

  ensure_dir(path);
  std::ofstream file(path);
  file << std::fixed << std::setprecision(6);

  file << "F0/θ";
  for (const auto &col : col_labels) {
    file << ", " << col;
  }
  file << '\n';

  for (size_t r = 0; r < results.size(); ++r) {
    file << row_labels[r];
    for (double val : results[r]) {
      file << ", " << val;
    }
    file << '\n';
  }
}

void ensure_dir(const std::string &path) {
  auto parent = std::filesystem::path(path).parent_path();
  std::filesystem::create_directories(parent);
}
