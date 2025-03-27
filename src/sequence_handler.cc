#pragma once

#include "sequence_handler.h"
#include "ksw2/ksw2.h"
#include "sam.h"
#include "spoa/include/spoa/spoa.hpp"
#include "utils.h"

#include <algorithm>
#include <assert.h>
#include <cinttypes>
#include <cmath>
#include <cstring>
#include <memory>
#include <numeric>
#include <string_view>
#include <unordered_map>

const int INF = std::numeric_limits<int>::max() / 2;

int match = 1, mis = 2, gap_open = 2, gap_ext = 1;
int8_t mat[25] = {1,  -2, -2, -2, 0,  -2, 1, -2, -2, 0, -2, -2, 1,
                  -2, 0,  -2, -2, -2, 1,  0, 0,  0,  0, 0,  0};

namespace factor {

const double error_rate = 0.01;
const int lower_bound = 1000, upper_bound = 5000;
const int gap_delta_threshold = 50;

const int min_anchor_size = 5, max_anchor_size = 20;

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

bool check_gap_set_is_valid(const std::vector<int> &index) {
  std::vector<int> gap_set;
  for (int i = 0; i < static_cast<int>(index.size()) - 1; ++i) {
    gap_set.push_back(index[i + 1] - index[i]);
    //  std::cerr << index[i + 1] - index[i] << ' ';
  }
  // std::cerr << std::endl;
  int len = gap_set.size();
  if (len < 1) {
    return false;
  }
  int avg = accumulate(gap_set.begin(), gap_set.end(), 0) / len;
  if (!is_gap_valid(avg)) {
    return false;
  }
  for (int gap : gap_set) {
    if (!is_gap_valid(gap) || std::abs(gap - avg) > gap_delta_threshold) {
      return false;
    }
  }
  return true;
}

std::vector<std::pair<int, std::vector<int>>>
collect_suspect_region(std::vector<int> &index) {
  std::sort(index.begin(), index.end());
  std::vector<std::pair<int, std::vector<int>>> suspect_region;
  int len = index.size();
  for (int i = 0; i < len; ++i) {
    for (int j = std::min(max_anchor_size, len);
         j >= min_anchor_size && i + j - 1 < index.size(); --j) {
      std::vector<int> anchor_index;
      for (int k = 0; k < j; ++k) {
        anchor_index.push_back(index[i + k]);
        //  std::cerr << index[i + k] << ' ';
      }
      // std::cerr << std::endl;
      if (check_gap_set_is_valid(anchor_index)) {
        suspect_region.emplace_back((index[i + j] - index[i]) / (j - 1),
                                    std::move(anchor_index));
        break;
      }
    }
  }
  return suspect_region;
}

std::vector<int> deal_with_gap_values(std::vector<int> &gap_values) {
  std::sort(gap_values.begin(), gap_values.end());
  int len = gap_values.size();
  std::vector<int> partition_index{-1};
  for (int i = 0; i < len - 1; ++i) {
    if (gap_values[i + 1] - gap_values[i] > gap_delta_threshold) {
      partition_index.push_back(i);
    }
  }
  partition_index.push_back(len);
  std::vector<int> result;
  int partition_index_size = partition_index.size();
  for (int i = 0; i < partition_index_size - 1; ++i) {
    int start = partition_index[i] + 1;
    int end = partition_index[i + 1];
    int avg =
        accumulate(gap_values.begin() + start, gap_values.begin() + end, 0) /
        (end - start);
    result.push_back(avg);
  }
  return result;
}

} // namespace factor

namespace alignment {

std::vector<int> ksw2_get_xid(uint32_t *cigar, int n_cigar,
                              const uint8_t *query, const uint8_t *target) {
  int qi = 0;
  int ti = 0;
  std::vector<int> xid(4);
  for (int i = 0; i < n_cigar; ++i) {
    int op = cigar[i] & 0xf;
    int len = cigar[i] >> 4;
    if (op == 0) {
      // match/mismatch
      for (int j = 0; j < len; ++j) {
        if (query[qi + j] == target[ti + j]) {
          xid[0]++;
        } else {
          xid[3]++;
        }
      }
      qi += len;
      ti += len;
    } else if (op == 1) {
      xid[1] += len;
      qi += len;
    } else if (op == 2) {
      xid[2] += len;
      ti += len;
    } else {
      return {};
    }
  }
  return xid;
}

std::pair<int, int> alignment(uint8_t *query, int qlen, uint8_t *target,
                              int tlen) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
  int *n_cigar;
  *n_cigar = ez.n_cigar;
  uint32_t **cigar;
  *cigar = ez.cigar;
  std::vector<int> xid = ksw2_get_xid(*cigar, *n_cigar, query, target);
  return {xid[0], accumulate(xid.begin(), xid.end(), 0)};
}

std::pair<int, int> alignment(const std::string &query_string,
                              const std::string &target_string) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  uint8_t *query = alloc_uint8_t(query_string);
  uint8_t *target = alloc_uint8_t(target_string);
  int qlen = static_cast<int>(query_string.size());
  int tlen = static_cast<int>(target_string.size());
  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
  int *n_cigar;
  *n_cigar = ez.n_cigar;
  uint32_t *cigar;
  cigar = ez.cigar;

  std::vector<int> xid = ksw2_get_xid(cigar, *n_cigar, query, target);
  free(query);
  free(target);
  return {xid[0], accumulate(xid.begin(), xid.end(), 0)};
}

int ksw2_backtrack_left_end(int n_cigar, uint32_t *cigar, int qlen, int tlen,
                            int q_left_ext) {
  int t_left_ext = 0, i, j, op, len;
  int q_remain_len = q_left_ext;
  for (i = n_cigar - 1; i >= 0; --i) {
    op = cigar[i] & 0xf;
    len = cigar[i] >> 4;
    if (op == 0) { // MATCH/MISMATCH
      if (len >= q_remain_len) {
        t_left_ext += q_remain_len;
        return t_left_ext;
      } else {
        t_left_ext += len;
        q_remain_len -= len;
      }
    } else if (op == 1) { // INSERTION
      if (len >= q_remain_len) {
        return t_left_ext;
      } else {
        q_remain_len -= len;
      }
    } else if (op == 2) { // DELETION
      t_left_ext += len;
    }
  }
  if (q_remain_len > 0) {
    std::cout << "Error: unmatched cigar and q_left_ext." << std::endl;
  }
  return t_left_ext;
}

int ksw2_backtrack_right_end(int n_cigar, uint32_t *cigar, int qlen, int tlen,
                             int q_right_ext) {
  int t_right_ext = 0, i, j, op, len;
  int q_remain_len = q_right_ext;
  for (i = 0; i < n_cigar; ++i) {
    op = cigar[i] & 0xf;
    len = cigar[i] >> 4;
    if (op == 0) { // MATCH/MISMATCH
      if (len >= q_remain_len) {
        t_right_ext += q_remain_len;
        return t_right_ext;
      } else {
        t_right_ext += len;
        q_remain_len -= len;
      }
    } else if (op == 1) { // INSERTION
      if (len >= q_remain_len) {
        return t_right_ext;
      } else {
        q_remain_len -= len;
      }
    } else if (op == 2) { // DELETION
      t_right_ext += len;
    }
  }
  if (q_remain_len > 0) {
    std::cout << "Error: unmatched cigar and q_right_ext." << std::endl;
  }
  return t_right_ext;
}

int extend_left_boundary(uint8_t *query, int qlen, uint8_t *target, int tlen) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
  int *n_cigar;
  *n_cigar = ez.n_cigar;
  uint32_t *cigar;
  cigar = ez.cigar;
  std::vector<int> xid = ksw2_get_xid(cigar, *n_cigar, query, target);
  return ksw2_backtrack_left_end(*n_cigar, cigar, qlen, tlen, xid[0]);
}

int extend_right_boundary(uint8_t *query, int qlen, uint8_t *target, int tlen) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
  int *n_cigar;
  *n_cigar = ez.n_cigar;
  uint32_t *cigar;
  cigar = ez.cigar;
  std::vector<int> xid = ksw2_get_xid(cigar, *n_cigar, query, target);
  return ksw2_backtrack_right_end(*n_cigar, cigar, qlen, tlen, xid[0]);
}

} // namespace alignment

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

  std::vector<std::pair<int, std::vector<int>>> suspect_region;
  std::vector<int> gap_values;
  for (int i = 0; i < sam->size(); ++i) {
    auto right_index_i = sam->get_right_index_by_node(i);
    if (right_index_i.size() < factor::min_anchor_size) {
      continue;
    }
    auto current_suspect_region = factor::collect_suspect_region(right_index_i);
    if (current_suspect_region.size() > 0) {
      std::cerr << "current_suspect_region.size() = "
                << current_suspect_region.size() << std::endl;
    }
    for (const auto &[gap, anchor_index] : current_suspect_region) {
      suspect_region.emplace_back(gap, std::move(anchor_index));
      gap_values.push_back(gap);
    }
  }

  std::cerr << "gap_values.size() = " << gap_values.size() << std::endl;

  std::cerr << "eeeeeeeee" << std::endl;

  std::vector<std::vector<std::vector<int>>> bucket_of_index_partition(
      gap_values.size());
  std::vector<std::vector<int>> bucket_of_index(gap_values.size());

  std::sort(
      suspect_region.begin(), suspect_region.end(),
      [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });

  int suspect_region_size = suspect_region.size();
  for (int i = 0; i < suspect_region_size; ++i) {
    const auto &[gap, anchor_index] = suspect_region[i];
    int gap_index =
        std::lower_bound(gap_values.begin(), gap_values.end(), gap) -
        gap_values.begin();
    int dis1 = gap_index > 0 ? std::abs(gap - gap_values[gap_index - 1]) : INF;
    int dis2 = gap_index < gap_values.size() - 1
                   ? std::abs(gap - gap_values[gap_index + 1])
                   : INF;
    if (dis1 < dis2) {
      bucket_of_index_partition[gap_index - 1].push_back(
          std::move(anchor_index));
      for (int index : anchor_index) {
        bucket_of_index[gap_index - 1].push_back(index);
      }
    } else {
      bucket_of_index_partition[gap_index].push_back(std::move(anchor_index));
      for (int index : anchor_index) {
        bucket_of_index[gap_index].push_back(index);
      }
    }
  }

  std::cerr << "bucket_of_index_partition.size() = "
            << bucket_of_index_partition.size() << std::endl;

  for (auto &vec : bucket_of_index) {
    std::sort(vec.begin(), vec.end(),
              [](const auto &lhs, const auto &rhs) { return lhs < rhs; });
  }

  // 暂时不需要用，用于确定真正的区间大小
  for (auto &vec : bucket_of_index_partition) {
    std::sort(vec.begin(), vec.end(),
              [](const auto &lhs, const auto &rhs) { return lhs[0] < rhs[0]; });
  }

  std::vector<std::pair<int, int>> raw_estimate_unit_region;

  // 排好序后，枚举每个端点求一个区间和
  for (int i = 0; i < gap_values.size(); ++i) {
    int unit_size = gap_values[i];
    int k = 0;
    int len = bucket_of_index[i].size();
    for (int j = 0; j < len; ++j) {
      int start_position = bucket_of_index[i][j];
      int end_position = start_position + unit_size;
      while (k < len && bucket_of_index[i][k] < end_position) {
        k++;
      }
      int num_of_index = k - j;

      // 如果区间内点数>= unit_size * error_rate
      if (num_of_index >= (int)(unit_size * factor::error_rate)) {
        raw_estimate_unit_region.emplace_back(start_position, unit_size);
      }
    }
  }

  std::sort(
      raw_estimate_unit_region.begin(), raw_estimate_unit_region.end(),
      [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });

  int max_r_boundary = -1;
  int max_l_boundary = -1;

  std::unordered_map<int, std::pair<int, int>> max_covered_region;

  auto already_covered = [&max_l_boundary, &max_r_boundary](int start_position,
                                                            int end_position,
                                                            int unit_size) {
    if (end_position <= max_r_boundary && start_position >= max_l_boundary) {
      if (end_position > max_r_boundary) {
        max_r_boundary = end_position;
        max_l_boundary = start_position;
      }
      return true;
    }
    if (end_position > max_r_boundary) {
      max_r_boundary = end_position;
      max_l_boundary = start_position;
    }
    return false;
  };

  std::cerr << "114514111111" << std::endl;

  std::vector<std::tuple<int, int, int>> raw_estimate_interval;
  for (const auto &[start_position, unit_size] : raw_estimate_unit_region) {
    int end_position = start_position + unit_size;
    if (already_covered(start_position, end_position, unit_size)) {
      continue;
    }
    // 向左延伸
    int raw_left_boundary = start_position;
    int raw_right_boundary = end_position;
    {
      for (int i = start_position; i >= unit_size; i -= unit_size) {
        // 做一下alignment，如果相似度太小就break
        auto [match, length] = alignment::alignment(
            s + i - unit_size, unit_size, s + start_position, unit_size);
        std::cerr << "length = " << length << " " << "match = " << match
                  << std::endl;
        if (1. - ((double)(match) / (double)(length)) <= factor::error_rate) {
          break;
        }
        raw_left_boundary = i - unit_size;
      }
    }
    // 向右延伸
    {
      for (int i = end_position; i + unit_size <= len; i += unit_size) {
        // 做一下alignment，如果相似度太小就break
        auto [match, length] = alignment::alignment(
            s + i, unit_size, s + start_position, unit_size);
        if (1. - ((double)match / (double)(length)) <= factor::error_rate) {
          break;
        }
        raw_right_boundary = i + unit_size;
      }
    }
    if (max_covered_region.find(start_position) == max_covered_region.end() ||
        max_covered_region[start_position].second < raw_right_boundary) {
      max_covered_region[start_position] = {raw_left_boundary,
                                            raw_right_boundary};
    }
    raw_estimate_interval.emplace_back(unit_size, raw_left_boundary,
                                       raw_right_boundary);
  }

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3); // linear gaps
  spoa::Graph graph{};

  for (auto [unit_size, left_boundary, right_boundary] :
       raw_estimate_interval) {
    std::cerr << "l = " << left_boundary << " " << "r = " << right_boundary
              << " " << "size = " << unit_size << std::endl;
    int seg_count =
        (right_boundary - left_boundary + unit_size - 1) / unit_size;

    std::vector<std::string> strings(seg_count);
    for (int i = left_boundary; i + unit_size - 1 <= right_boundary;
         i += unit_size) {
      std::string sequence("");
      for (int j = i; j < i + unit_size; ++j) {
        sequence += safeNumberToDnaChar(s[j]);
      }
      strings.emplace_back(std::move(sequence));
    }

    for (auto &sequence : strings) {
      auto alignment = alignment_engine->Align(sequence, graph);
      graph.AddAlignment(alignment, sequence);
    }

    auto consensus = graph.GenerateConsensus();
    int consensus_len = static_cast<int>(consensus.size());

    uint8_t *cons = alloc_uint8_t(consensus);
    // 向左拓展边界
    int left_ext = alignment::extend_left_boundary(
        s + std::max(0, left_boundary - unit_size),
        std::min(unit_size, left_boundary), cons, consensus_len);

    // 向右拓展边界
    int right_ext = alignment::extend_right_boundary(
        s + std::min(len, right_boundary + unit_size),
        std::min(unit_size, len - right_boundary), cons, consensus_len);
    left_boundary -= left_ext;
    right_boundary += right_ext;
    double all_match = 0;
    double all_total = 0;
    for (int i = left_boundary; i < right_boundary; i += unit_size) {
      auto [match, total] =
          alignment::alignment(s + i, unit_size, cons, consensus_len);
      all_match += match;
      all_total += total;
    }
    std::cout << "[" << left_boundary << "," << right_boundary << "]" << " "
              << all_match / all_total << ' ' << consensus << std::endl;
    free(cons);
  }
}