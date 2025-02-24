#pragma once

#include <string>

class Param {
public:
  int lower_bound_size;
  int gap_threshold;
  double overlap_threshold;
  double gap_ratio;
};

bool filter_by_seed(const std::string &s, int pos);