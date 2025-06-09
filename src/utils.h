#pragma once

#include <cstdint>
#include <string>
#include <vector>

class Param {
public:
  int lower_bound_size;
  int gap_threshold;
  double overlap_threshold;
  double gap_ratio;

  bool verbose;
  char *input_filename;
  char *result_filename;
  std::vector<std::string> filter;
};

bool filter_by_seed(const std::string &s, int pos);
char safeNumberToDnaChar(int num, char default_char = 'N');
uint8_t *alloc_uint8_t(const std::string &sequence);
std::vector<std::string> split(const std::string &s, char delimiter);
bool check_substring(const char* text, const std::string& substring);