#pragma once
#include "utils.h"

#include <iostream>

#define overlap_threshold 0.5

const char *seeds[] = {
    "RRRRRRRYYRYY", "RRYRRYYYYYYY", "RRRRRRYYYRYY", "RRYRRRYYYYYY",
    "RRRRRYRYYRYY", "RRYRRYRYYYYY", "RRRRRYYRYRYY", "RRYRYRRYYYYY",
    "RRRRRYYYYRYY", "RRYRRRRYYYYY", "RRRRYRYYYRYY", "RRYRRRYRYYYY",
    "RRRRYYYRRYYY", "RRRYYRRRYYYY", "RRRRYYYRYRYY", "RRYRYRRRYYYY",
    "RRRYRRRYYRYY", "RRYRRYYYRYYY", "RRRYRYRYYRYY", "RRYRRYRYRYYY",
    "RRRYRYYRRYYY", "RRRYYRRYRYYY", "RRRYRYYRYRYY", "RRYRYRRYRYYY",
    "RRRYYRRYYRYY", "RRYRRYYRRYYY", "RRRYYRYRYRYY", "RRYRYRYRRYYY",
    "RRYRRRRYYRYY", "RRYRRYYYYRYY", "RRYRRYYRYRYY", "RRYRYRRYYRYY"};

const int seed_size = 32;

bool filter_by_seed(const std::string &s, int pos) {
  std::cerr << "s = " << s << std::endl;
  for (int i = 0; i < seed_size; ++i) {
    bool match = true;
    for (int j = 0; j < 12; ++j) {
      if ((seeds[i][j] == 'R' && (s[pos + j] != '0' || s[pos + j] != '2')) ||
          (seeds[i][j] == 'Y' && (s[pos + j] != '1' || s[pos + j] != '3'))) {
        match = false;
        break;
      }
    }
    if (match) {
      return true;
    }
  }
  return false;
}