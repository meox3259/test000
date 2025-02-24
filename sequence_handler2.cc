#pragma once

#include "sequence_handler.h"
#include "sam.h"
#include "utils.h"

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

void solve(const std::string &s, const Param &opt) {
  int len = s.size();
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

  sort(
      index_set.begin(), index_set.end(),
      [](const auto &lhs, const auto &rhs) { return lhs.back() < rhs.back(); });

  std::cerr << "set size = " << index_set.size() << std::endl;
  std::cerr << "start overlap" << std::endl;
  while (true) {
    int sz = index_set.size();
    std::cerr << "sz = " << sz << std::endl;
    std::vector<std::pair<int, int>> index_pair;
    std::vector<bool> mark(sz);
    for (int i = 0; i < sz; ++i) {
      if (mark[i]) {
        continue;
      }
      for (int j = i + 1; j < sz; ++j) {
        if (mark[j]) {
          continue;
        }
        if (get_overlap(index_set[i], index_set[j], opt.overlap_threshold)) {
          mark[i] = true;
          mark[j] = true;
          index_pair.emplace_back(i, j);
          break;
        }
      }
    }

    if (index_pair.empty()) {
      break;
    }

    std::vector<std::vector<int>> new_index_set;
    for (auto &[i, j] : index_pair) {
      std::vector<int> new_set;
      for (int item : index_set[i]) {
        new_set.push_back(item);
      }
      for (int item : index_set[j]) {
        new_set.push_back(item);
      }
      sort(new_set.begin(), new_set.end());
      new_index_set.push_back(std::move(new_set));
    }
    for (int i = 0; i < sz; ++i) {
      if (!mark[i]) {
        new_index_set.push_back(std::move(index_set[i]));
      }
    }

    std::swap(index_set, new_index_set);
  }

  std::cerr << "size = " << index_set.size() << std::endl;
  for (const auto &vec : index_set) {
    std::cerr << "l = " << vec[0] << " " << "r = " << vec.back() << " "
              << "size = " << vec.size() << std::endl;
  }
}