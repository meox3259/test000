#include "utils.h"

#include <cstdint>
#include <iostream>
#include <sys/types.h>

#define overlap_threshold 0.5

enum class GENE { A, C, G, T };

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

char safeNumberToDnaChar(int num, char default_char) {
  if (num >= 0 && num <= 3) {
    static const char mapping[] = {'A', 'C', 'G', 'T'};
    return mapping[num];
  }
  return default_char;
}

uint8_t *alloc_uint8_t(const std::string &sequence) {
  int len = static_cast<int>(sequence.size());
  uint8_t *data = new uint8_t[len];
  return data;
}

bool N_filter(const std::string &sequence) {
  int N_count = 0;
  for (char c : sequence) {
    N_count += (c == 'N');
  }
  return N_count > 0.5 * sequence.size();
}

template <typename T>
bool N_filter(T *sequence, int len) {
  int N_count = 0;
  for (int i = 0; i < len; ++i) {
    N_count += (sequence[i] == 'N');
  }
  return N_count > 0.5 * len;
}