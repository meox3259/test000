#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "handle_file.h"
#include "logger.h"
#include "sequence_handler.h"
#include "utils.h"

bool LogStream::verbose = false;

int main(int argc, char *argv[]) {
  static struct option long_options[] = {
      {"file", required_argument, 0, 'f'},   {"verbose", no_argument, 0, 'v'},
      {"number", required_argument, 0, 'n'}, {"help", no_argument, 0, 'h'},
      {"result", required_argument, 0, 'r'}, {0, 0, 0, 0}};

  int opt;
  Param param{.lower_bound_size = 100,
              .gap_threshold = 4000,
              .overlap_threshold = 0.5,
              .gap_ratio = 0.8,
              .verbose = false};
  while ((opt = getopt(argc, argv, "vf:r:")) != -1) {
    switch (opt) {
    case 'v':
      param.verbose = true;
      LogStream::verbose = true;
      break;
    case 'f':
      param.input_filename = optarg;
      break;
    case 'r':
      param.result_filename = optarg;
      break;
    }
  }

  std::ofstream ofs(param.result_filename);
  LOG << "input_filename = " << param.input_filename;

  FILE *fp = init_file(param.input_filename);
  Read *Read = return_read(fp);
  LOG << "len = " << Read->len;

  solve(ofs, Read->Read, Read->len, param);
  delete Read;
}