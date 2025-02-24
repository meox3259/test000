#pragma once

#include "ksw2_align.h"

int extend_to_left(int n_cigar, uint32_t *cigar, int left_remain_length) {
  int extend_to_left = 0;
  for (int i = n_cigar - 1; i >= 0; --i) {
    int op = cigar[i] & 0xf;
    int len = cigar[i] >> 4;
    switch (op) {
    case 0: // match/mismatch
      if (len >= left_remain_length) {
        extend_to_left += left_remain_length;
        return extend_to_left;
      } else {
        extend_to_left += len;
        left_remain_length -= len;
      }
      break;
    case 1: // insertion
      if (len >= left_remain_length) {
        return extend_to_left;
      } else {
        left_remain_length -= len;
      }
      break;
    case 2: // deletion
      extend_to_left += len;
      break;
    }
  }
}