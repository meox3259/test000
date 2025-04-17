#include "sequence_handler.h"
#include "ksw2/ksw2.h"
#include "logger.h"
#include "radix_sort.h"
#include "spoa/include/spoa/spoa.hpp"
#include "suffix_array.h"
#include "utils.h"

#include <algorithm>
#include <array>
#include <assert.h>
#include <cinttypes>
#include <cmath>
#include <cstring>
#include <fstream>
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
const int min_anchor_size = 3, max_anchor_size = 20;
const int kmer_size = 30;

bool is_gap_valid(int gap) { return gap >= lower_bound && gap <= upper_bound; }

void collect_suspect_region(
    std::vector<int> index,
    std::vector<std::pair<int, std::vector<int>>> &suspect_region) {
  suspect_region.clear();
  if (index.size() < min_anchor_size) {
    return;
  }
  std::sort(index.begin(), index.end());
  assert(std::is_sorted(index.begin(), index.end()));
  int len = index.size();
  int sum_of_gap = 0;
  int window_size = min_anchor_size;
  std::vector<int> anchor_index;
  anchor_index.reserve(max_anchor_size);
  for (int i = 0; i < len;) {
    anchor_index.clear();
    int anchor_size;
    int sum_of_gap = 0;
    int min_gap = INF;
    int max_gap = -INF;
    int avg_gap = INF;
    anchor_index.push_back(index[i]);
    for (anchor_size = 1;
         anchor_size <= max_anchor_size && i + anchor_size < len;
         ++anchor_size) {
      int gap = index[i + anchor_size] - index[i + anchor_size - 1];
      sum_of_gap += gap;
      min_gap = std::min(min_gap, gap);
      max_gap = std::max(max_gap, gap);
      anchor_index.push_back(index[i + anchor_size]);
      avg_gap = sum_of_gap / anchor_size;
      if (std::abs(avg_gap - min_gap) > gap_delta_threshold ||
          std::abs(avg_gap - max_gap) > gap_delta_threshold) {
        break;
      }
    }

    if (avg_gap >= lower_bound && avg_gap <= upper_bound &&
        anchor_index.size() > min_anchor_size) {
      suspect_region.emplace_back(avg_gap, std::move(anchor_index));
    }

    i += anchor_size;
  }
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
  std::vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
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

  std::vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
  delete[] query;
  delete[] target;
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
  std::vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
  return ksw2_backtrack_left_end(ez.n_cigar, ez.cigar, qlen, tlen, xid[0]);
}

int extend_right_boundary(uint8_t *query, int qlen, uint8_t *target, int tlen) {
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  int w = -1, zdrop = -1, end_bonus = 0, flag = 0;

  ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w,
                zdrop, end_bonus, flag, &ez);
  std::vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
  return ksw2_backtrack_right_end(ez.n_cigar, ez.cigar, qlen, tlen, xid[0]);
}

} // namespace alignment

void solve(std::ofstream &ofs, uint8_t *s, int len, const Param &opt) {
  std::vector<int> sequence;
  for (int i = 0; i < len; ++i) {
    sequence.push_back(s[i]);
  }

  LOG << "Start to build suffix array: size = " << sequence.size();

  auto suffix_array = yosupo::suffix_array(sequence);
  auto lcp_array = yosupo::lcp_array(sequence, suffix_array);

  LOG << "End to build suffix array";

  std::vector<std::vector<int>> index_set;
  std::vector<std::pair<int, std::vector<int>>> suspect_region;
  std::vector<int> gap_values;
  std::vector<std::pair<int, std::vector<int>>> current_suspect_region;
  std::vector<int> right_index{suffix_array[0]};

  LOG << "Start to build suspect region";
  for (int i = 0; i < len - 1; ++i) {
    if (lcp_array[i] < factor::kmer_size) {
      if (right_index.size() <= factor::min_anchor_size) {
        continue;
      }
      factor::collect_suspect_region(right_index, current_suspect_region);
      for (const auto &[gap, anchor_index] : current_suspect_region) {
        suspect_region.emplace_back(gap, std::move(anchor_index));
        gap_values.push_back(gap);
      }
      right_index.clear();
    }
    right_index.push_back(suffix_array[i + 1]);
  }

  if (right_index.size() > factor::min_anchor_size) {
    factor::collect_suspect_region(right_index, current_suspect_region);
    for (const auto &[gap, anchor_index] : current_suspect_region) {
      suspect_region.emplace_back(gap, std::move(anchor_index));
      gap_values.push_back(gap);
    }
  }
  LOG << "End to build suspect region: suspect_region.size() = "
      << suspect_region.size();

  std::vector<std::vector<std::vector<int>>> bucket_of_index_partition(
      gap_values.size());
  std::vector<std::vector<int>> bucket_of_index(gap_values.size());

  std::function<int(const std::pair<int, std::vector<int>> &)> get_index =
      [](const std::pair<int, std::vector<int>> &item) -> int {
    return item.first;
  };

  LOG << "start radix_sort of suspect_region";

  radix_sort(suspect_region, std::forward<decltype(get_index)>(get_index), len);

  LOG << "end radix_sort of suspect_region";

  int suspect_region_size = suspect_region.size();
  for (int i = 0; i < suspect_region_size; ++i) {
    // remove the log
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

  LOG << "bucket_of_index_partition.size() = "
      << bucket_of_index_partition.size();

  for (auto &vec : bucket_of_index) {
    std::sort(vec.begin(), vec.end(),
              [](const auto &lhs, const auto &rhs) { return lhs < rhs; });
    assert(std::is_sorted(vec.begin(), vec.end()));
  }

  // 暂时不需要用，用于确定真正的区间大小
  // for (auto &vec : bucket_of_index_partition) {
  //   std::sort(vec.begin(), vec.end(),
  //             [](const auto &lhs, const auto &rhs) { return lhs[0] < rhs[0];
  //             });
  // }

  LOG << "gap_values.size() = " << gap_values.size();
  std::vector<std::pair<int, int>> raw_estimate_unit_region;

  // 排好序后，枚举每个端点求一个区间和
  for (int i = 0; i < gap_values.size(); ++i) {
    int unit_size = gap_values[i];
    int k = 0;
    int len_of_bucket = bucket_of_index[i].size();
    for (int j = 0; j < len_of_bucket; ++j) {
      int start_position = bucket_of_index[i][j];
      int end_position = start_position + unit_size;
      while (k < len_of_bucket && bucket_of_index[i][k] < end_position) {
        k++;
      }
      int num_of_index = k - j;

      // error_rate -> error model?
      // here calculate chains that fall in this interval[p, p + unit_size]
      // smarter algorithm
      if (num_of_index >= (int)(unit_size * factor::error_rate)) {
        raw_estimate_unit_region.emplace_back(start_position, unit_size);
      }
    }
  }

  LOG << "raw_estimate_unit_region.size() = "
      << raw_estimate_unit_region.size();

  radix_sort(
      raw_estimate_unit_region,
      [](const auto &item) -> int { return item.first; }, len);

  int max_r_boundary = -1;
  int max_l_boundary = -1;

  std::unordered_map<int, std::pair<int, int>> max_covered_region;

  auto already_covered = [&max_l_boundary, &max_r_boundary](int start_position,
                                                            int end_position,
                                                            int unit_size) {
    if (end_position <= max_r_boundary + unit_size &&
        start_position >= max_l_boundary) {
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

  LOG << "raw_estimate_unit_region.size() = "
      << raw_estimate_unit_region.size();

  int cnt = 0;

  std::vector<std::tuple<int, int, int>> raw_estimate_interval;
  for (const auto &[start_position, unit_size] : raw_estimate_unit_region) {
    int end_position = start_position + unit_size;
    if (already_covered(start_position, end_position, unit_size)) {
      continue;
    }
    cnt++;
    // 向左延伸
    int raw_left_boundary = start_position;
    int raw_right_boundary = end_position;
    {
      for (int i = start_position; i >= unit_size; i -= unit_size) {
        // 做一下alignment，如果相似度太小就break
        auto [match, length] = alignment::alignment(
            s + i - unit_size, unit_size, s + i, unit_size);
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
            s + i, unit_size, s + i - unit_size, unit_size);
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

  LOG << "cnt = " << cnt;
  LOG << "raw_estimate_interval.size() = " << raw_estimate_interval.size();

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3); // linear gaps
  spoa::Graph graph{};

  ofs << "fasta:" << std::endl;

  for (auto [unit_size, left_boundary, right_boundary] :
       raw_estimate_interval) {
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

    LOG << "start spoa alignment";

    for (auto &sequence : strings) {
      auto alignment = alignment_engine->Align(sequence, graph);
      graph.AddAlignment(alignment, sequence);
    }

    LOG << "end spoa alignment";

    auto consensus = graph.GenerateConsensus();
    int consensus_len = static_cast<int>(consensus.size());
    uint8_t *cons = alloc_uint8_t(consensus);

    // 向左拓展边界，这里需要把序列倒过来然后匹配
    std::string reverse_sequence("");
    for (int i = 0; i < std::min(unit_size, left_boundary); ++i) {
      reverse_sequence += safeNumberToDnaChar(s[left_boundary - i]);
    }

    std::reverse(consensus.begin(), consensus.end());
    std::string reverse_consensus = consensus;
    std::reverse(reverse_consensus.begin(), reverse_consensus.end());

    uint8_t *reverse_cons = alloc_uint8_t(reverse_consensus);
    uint8_t *reverse_seq = alloc_uint8_t(reverse_sequence);

    int qlen = std::min(unit_size, left_boundary);
    int left_ext = alignment::extend_left_boundary(reverse_seq, qlen,
                                                   reverse_cons, consensus_len);

    delete[] reverse_cons;
    delete[] reverse_seq;

    // 向右拓展边界
    int right_ext = alignment::extend_right_boundary(
        s + std::min(len, right_boundary + unit_size),
        std::min(unit_size, len - right_boundary), cons, consensus_len);
    left_boundary -= left_ext;
    right_boundary += right_ext;
    // tmp
    left_boundary = std::max(0, left_boundary);

    double all_match = 0;
    double all_total = 0;
    for (int i = left_boundary; i < right_boundary; i += unit_size) {
      auto [match, total] = alignment::alignment(
          s + i, std::min(unit_size, right_boundary - i), cons, consensus_len);
      all_match += match;
      all_total += total;
    }
    ofs << "[" << left_boundary << "," << right_boundary << "]" << " "
        << all_match / all_total << ' ' << consensus << std::endl;
    delete[] cons;
  }
}