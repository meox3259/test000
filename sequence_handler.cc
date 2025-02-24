#pragma once

#include "sequence_handler.h"
#include "ksw2/ksw2.h"
#include "ksw2_align.h"
#include "sam.h"

#include <algorithm>
#include <assert.h>
#include <memory>

bool calculate_sparsity(const std::vector<int> &index, int gap_threshold,
                        double gap_ratio) {
  if (!std::is_sorted(index.begin(), index.end())) {
    std::cerr << "Error: index is not sorted" << std::endl;
    return false;
  }

  int tot = 0;
  for (int i = 0; i < static_cast<int>(index.size()) - 1; ++i) {
    int gap = index[i + 1] - index[i];
    if (gap <= gap_threshold) {
      tot++;
    }
  }

  int tot_gap = static_cast<int>(index.size()) - 1;

  return static_cast<double>(tot) / static_cast<double>(tot_gap) >= gap_ratio;
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

void solve(int *s, int len, const Param &opt) {
  std::shared_ptr<SuffixAutomation> sam = std::make_shared<SuffixAutomation>(s);
  sam->get_right_index(opt.lower_bound_size);

  std::cerr << "start sam calc" << std::endl;
  std::vector<std::vector<int>> index_set;
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
  }

  std::cerr << "size = " << index_set.size() << std::endl;
  for (const auto &vec : index_set) {
    std::cerr << "l = " << vec[0] << " " << "r = " << vec.back() << " "
              << "size = " << vec.size() << std::endl;
    int seg_count = vec.size();
    for (int i = 0; i < seg_count - 1; ++i) {
      int *start = s + vec[i];
      int len = vec[i + 1] - vec[i] + 1;
    }
    int total_length = vec.back() - vec[0] + 1;
    int cur_seg = vec.size();
    int unit_length = (total_length + cur_seg - 1) / cur_seg;
    ksw2_global_with_cigar(s + vec[0], unit_length, bseq + e1 - k + 1,
                           e2 - e1 + k, &n_cigar, &cigar);

    int n_seqs;
    cons_len = abpoa_gen_cons(ab, abpt, bseq, seq_len, par_pos + i, j - i,
                              cons_bseq, cons_qual, mtp, &n_seqs);

    int ksw2_global_with_cigar(const uint8_t *query, int qlen,
                               const uint8_t *target, int tlen, int *n_cigar,
                               uint32_t **cigar);

    extend_to_left(n_cigar);
  }
}