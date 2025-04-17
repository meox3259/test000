#pragma once

#include <cstdint>
#include <string>

class Param {
public:
  int lower_bound_size;
  int gap_threshold;
  double overlap_threshold;
  double gap_ratio;

  bool verbose;
  char *input_filename;
  char *result_filename;
};

bool filter_by_seed(const std::string &s, int pos);
char safeNumberToDnaChar(int num, char default_char = 'N');
uint8_t *alloc_uint8_t(const std::string &sequence);