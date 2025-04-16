#include <cmath>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include "handle_file.h"
#include "logger.h"
#include "sequence_handler.h"
#include "utils.h"

bool LogStream::verbose = false;

int main(int argc, char *argv[]) {
  std::ofstream ofs("result.txt", std::ios::app);
  int opt;
  Param option{.lower_bound_size = 100,
               .gap_threshold = 4000,
               .overlap_threshold = 0.5,
               .gap_ratio = 0.8,
               .verbose = false};
  while ((opt = getopt(argc, argv, "vf:")) != -1) {
    switch (opt) {
    case 'v':
      option.verbose = true;
      LogStream::verbose = true;
      break;
    case 'f':
      option.filename = optarg;
      break;
    }
  }
  char *filename = argv[1];
  printf("name = %s\n", filename);
  FILE *fp = init_file(filename);
  Read *Read = return_read(fp);
  std::string s = "";
  for (int i = 0; i < Read->len; ++i) {
    s += (char)(Read->Read[i] + '0');
  }
  printf("len = %d\n", Read->len);
  solve(ofs, Read->Read, Read->len, option);
  delete Read;
}