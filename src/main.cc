#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "handle_file.h"
#include "sequence_handler.h"
#include "utils.h"

int main(int argc, char *argv[]) {
  std::ofstream ofs("result.txt", std::ios::app);
  Param opt{.lower_bound_size = 100,
            .gap_threshold = 4000,
            .overlap_threshold = 0.5,
            .gap_ratio = 0.8};
  char *filename = argv[1];
  printf("name = %s\n", filename);
  FILE *fp = init_file(filename);
  Read *Read = return_read(fp);
  std::string s = "";
  for (int i = 0; i < Read->len; ++i) {
    s += (char)(Read->Read[i] + '0');
  }
  printf("len = %d\n", Read->len);
  solve(ofs, Read->Read, Read->len, opt);
  //  solve(s, opt);
}