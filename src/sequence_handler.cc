#pragma once

#include "sequence_handler.h"
#include "ksw2/ksw2.h"
#include "sam.h"
#include "spoa/include/spoa/spoa.hpp"
#include "utils.h"

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <cstring>
#include <memory>

int match = 1, mis = 2, gap_open = 2, gap_ext = 1;
int8_t mat[25] = {1,  -2, -2, -2, 0,  -2, 1, -2, -2, 0, -2, -2, 1,
                  -2, 0,  -2, -2, -2, 1,  0, 0,  0,  0, 0,  0};

namespace factor {

int lower_bound = 1000, upper_bound = 3000;

bool is_gap_valid(int gap) { return gap >= lower_bound && gap <= upper_bound; }

bool calculate_sparsity(const std::vector<int> &index, int gap_threshold,
                        double gap_ratio) {
  if (!std::is_sorted(index.begin(), index.end())) {
    std::cerr << "Error: index is not sorted" << std::endl;
    return false;
  }

  int tot = 0;
  for (int i = 0; i < static_cast<int>(index.size()) - 1; ++i) {
    int gap = index[i + 1] - index[i];
    if (!(gap >= lower_bound && gap <= upper_bound)) {
      return false;
    }
  }
}

std::vector<int> get_suspect_region(const std::vector<int> &index) {
  int len = index.size();
  for (int i = 0; i < len; ++i) {
    std::vector<int> gap_size;
    for (int j = 1; j < 10; ++j) {
      gap_size.push_back(index[i + j] - index[i + j - 1]);
    }
    for (int gap : gap_size) {
      if
    }
  }
}

bool get_overlap(const std::vector<int> &l, const std::vector<int> &r,
                 double overlap_threshold) {
  if (l.empty() || r.empty()) {
    std::cerr << "Error: l or r is empty" << std::endl;
    return false;
  }
  int l0 = l[0], l1 = l.back(), r0 = r[0], r1 = r.back();
  int l_cross =
      upper_bound(r.begin(), r.end(), l1) - lower_bound(r.begin(), r.end(), l0);
  int r_cross =
      upper_bound(l.begin(), l.end(), r1) - lower_bound(l.begin(), l.end(), r0);
  int l_size = l.size(), r_size = r.size();
  double overlap_l = (double)l_cross / (double)(l_size);
  double overlap_r = (double)r_cross / (double)(r_size);
  return overlap_l >= overlap_threshold && overlap_r >= overlap_threshold;
}

} // namespace factor

void parse_index_set(const std::vector<int> &index) {
  int len = index.size();
  for (int i = 0; i < len; ++i) {
    std::vector<int> gap_size;
    for (int j = 1; j < 10; ++j) {
      gap_size.push_back(index[i + j] - index[i + j - 1]);
    }
  }
}

void solve(uint8_t *s, int len, const Param &opt) {
  std::string sequence("");
  for (int i = 0; i < len; ++i) {
    sequence += safeNumberToDnaChar(s[i]);
  }
  std::shared_ptr<SuffixAutomation> sam =
      std::make_shared<SuffixAutomation>(sequence);
  sam->get_right_index(opt.lower_bound_size);

  std::cerr << "start sam calc" << std::endl;
  std::vector<std::vector<int>> index_set;
  for (int i = 0; i < sam->size(); ++i) {
    auto right_index_set = sam->get_right_index_by_node(i);
    if (right_index_set.size() < 4) {
      continue;
    }
  }
  for (int i = 0; i < sam->size(); ++i) {
    auto right_index_i = sam->get_right_index_by_node(i);
    if (right_index_i.size() < 4) {
      continue;
    }
    if (calculate_sparsity(right_index_i, opt.gap_threshold, opt.gap_ratio)) {
      index_set.push_back(std::move(right_index_i));
    }
  }

  sort(index_set.begin(), index_set.end(),
       [](const auto &lhs, const auto &rhs) {
         return *lhs.begin() < *rhs.begin();
       });

  auto merge = [](std::vector<int> &lhs, std::vector<int> &rhs) {
    if (lhs.size() > rhs.size()) {
      swap(lhs, rhs);
    }
    for (int i : lhs) {
      rhs.emplace_back(i);
    }
    std::sort(rhs.begin(), rhs.end());
    lhs.clear();
  };

  std::cerr << "set size = " << index_set.size() << std::endl;
  std::cerr << "start overlap" << std::endl;
  for (int i = 0; i < index_set.size(); ++i) {
  }
  while (true) {
    int sz = index_set.size();
    int merge_count = 0;
    for (int i = 0; i < sz - 1; ++i) {
      if (get_overlap(index_set[i], index_set[i + 1], opt.overlap_threshold)) {
        merge(index_set[i], index_set[i + 1]);
        merge_count++;
      }
    }
    if (merge_count == 0) {
      break;
    }
    std::vector<std::vector<int>> new_index_set;
    for (int i = 0; i < sz; ++i) {
      if (!index_set[i].empty()) {
        new_index_set.emplace_back(index_set[i]);
      }
    }
    swap(index_set, new_index_set);
  }

  std::cerr << "size = " << index_set.size() << std::endl;
  for (const auto &vec : index_set) {
    std::cerr << "l = " << vec[0] << " " << "r = " << vec.back() << " "
              << "size = " << vec.size() << std::endl;
  }
  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3); // linear gaps
  spoa::Graph graph{};

  for (const auto &vec : index_set) {
    std::cerr << "l = " << vec[0] << " " << "r = " << vec.back() << " "
              << "size = " << vec.size() << std::endl;
    int seg_count = vec.size();

    std::vector<std::string> strings(seg_count);
    for (int i = 0; i < seg_count - 1; ++i) {
      std::string sequence("");
      for (int j = vec[i]; j < vec[i + 1]; ++j) {
        sequence += safeNumberToDnaChar(s[j]);
      }
      strings.emplace_back(std::move(sequence));
      std::cerr << "vec[i] = " << vec[i] << " vec[i + 1] = " << vec[i + 1];
      std::cerr << " i = " << i << " sequence = " << sequence << std::endl;
    }

    for (auto &sequence : strings) {
      auto alignment = alignment_engine->Align(sequence, graph);
      graph.AddAlignment(alignment, sequence);
    }

    auto consensus = graph.GenerateConsensus();
    int consensus_len = static_cast<int>(consensus.size());

    int total_length = vec.back() - vec[0] + 1;
    int cur_seg = vec.size();
    int unit_length = (total_length + cur_seg - 1) / cur_seg;

    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

    uint8_t *target = alloc_uint8_t(consensus);
    ksw_extz2_sse(0, std::min(consensus_len, vec[0]),
                  s + std::max(vec[0] - consensus_len, 0), consensus_len,
                  target, 5, mat, gap_open, gap_ext, w, zdrop, end_bonus, flag,
                  &ez);

    int *n_cigar;
    *n_cigar = ez.n_cigar;
    uint32_t **cigar;
    *cigar = ez.cigar;
  }
}