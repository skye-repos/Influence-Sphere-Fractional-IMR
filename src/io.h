#pragma once
#include <string>
#include <vector>

void write_csv(const std::string &path,
               const std::vector<std::string> &col_labels,
               const std::vector<std::string> &row_labels,
               const std::vector<std::vector<double>> &results);

void ensure_dir(const std::string &path);
