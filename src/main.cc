#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <unistd.h>

#include "handle_file.h"
#include "logger.h"
#include "sequence_handler.h"
#include "utils.h"
//#include "fasta_sequence.h"

bool LogStream::verbose = false;

class Timer {
  std::chrono::steady_clock::time_point start;

public:
  Timer() : start(std::chrono::steady_clock::now()) {}
  ~Timer() {
    auto end = std::chrono::steady_clock::now();
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "程序运行时间: " << duration.count() << " 毫秒" << std::endl;
  }
};

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

  Timer timer;

  std::ofstream ofs(param.result_filename);
  LOG << "input_filename = " << param.input_filename;

  FILE *fp = init_file(param.input_filename);
  while (true) {
    LOG << "start while";
    Read *Read = return_read(fp);
    LOG << "end return_read";
    if (Read == nullptr) {
      break;
    }
    if (Read->len == 0) {
      break;
    }
    LOG << "len = " << Read->len;
    LOG << "Read->ID = " << Read->ID;
    solve(ofs, Read->Read, Read->len, Read->is_N, param);
    delete Read;
    LOG << "end solve";
  }
}