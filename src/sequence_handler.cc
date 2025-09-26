#include "sequence_handler.h"
#include "edlib_align.h"
#include "fenwick_tree.h"
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

int prime_number[] = {
    -100000000, 1009, 1061, 1109, 1163, 1213, 1277, 1327, 1381,    1429,
    1481,       1531, 1579, 1627, 1693, 1741, 1789, 1847, 1901,    1949,
    1997,       2053, 2111, 2161, 2213, 2267, 2333, 2381, 2437,    2503,
    2551,       2609, 2657, 2707, 2767, 2819, 2879, 2927, 2999,    3049,
    3109,       3163, 3217, 3271, 3319, 3371, 3433, 3491, 3539,    3593,
    3643,       3691, 3739, 3793, 3847, 3907, 3967, 4000, 4050,    4100,
    4150,       4200, 4350, 4300, 4350, 4400, 4450, 4500, 4550,    4600,
    4650,       4700, 4750, 4800, 4850, 4900, 4950, 5000, +1000000};

namespace factor {

const double error_rate = 0.01;
const double kmer_rate = 0.05; // change
const double alignment_error_rate = 0.2;
const int lower_bound = 1000, upper_bound = 5000;
const int gap_delta_threshold = 50;
const int min_anchor_size = 3, max_anchor_size = 20;
const int kmer_size = 30;
const double N_threshold = 0.5;

bool is_gap_valid(int gap) { return gap >= lower_bound && gap <= upper_bound; }

void collect_suspect_region(
    std::vector<int> index,
    std::vector<std::pair<int, std::vector<int>>> &suspect_region) {
  suspect_region.clear();
  if (index.size() < min_anchor_size) {
    return;
  }
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
  for (int i = 0; i < ez.n_cigar; ++i) // print CIGAR
    printf("%d%c", ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
  putchar('\n');
  fflush(stdout);
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

int ksw2_global_with_cigar(const uint8_t *query, int qlen, const uint8_t *target, int tlen, int *n_cigar, uint32_t **cigar) {
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int w=-1, zdrop=-1, end_bonus=0, flag = 0;
    LOG << "ksw2_global_with_cigar start";
    LOG << "qlen = " << qlen << ", tlen = " << tlen;
    ksw_extz2_sse(0, qlen, query, tlen, target, 5, mat, gap_open, gap_ext, w, zdrop, end_bonus, flag, &ez);
#ifdef __DEBUG__
    print_cigar(ez.n_cigar, ez.cigar);
#endif
    LOG << "ksw2_global_with_cigar end";
    vector<int> xid = ksw2_get_xid(ez.cigar, ez.n_cigar, query, target);
    int iden_n = 0;
    if (!xid.empty()) {
        iden_n = xid[0];
    }
    *n_cigar = ez.n_cigar;
    *cigar = ez.cigar;
    // if (ez.cigar) free(ez.cigar);
    return iden_n;
}

} // namespace alignment

void solve(std::ofstream &ofs, uint8_t *seq, int len, int *is_N, const Param &opt) {
  LOG << "Start to build suffix array: size = " << len;
  std::vector<int> sequence;
  for (int i = 0; i < len; ++i) {
    sequence.push_back(seq[i]);
  }

  auto suffix_array = yosupo::suffix_array(sequence);
  auto lcp_array = yosupo::lcp_array(sequence, suffix_array);

  LOG << "End to build suffix array";

  std::vector<std::vector<int>> index_set;
  std::vector<int> gap_values;
  std::vector<std::pair<int, std::vector<int>>> current_suspect_region;
  std::vector<int> right_index{suffix_array[0]};

  LOG << "Start to build suspect region";
  std::vector<int> prefix_sum_N(len + 1);
  for (int i = 0; i < len; ++i) {
    prefix_sum_N[i + 1] = prefix_sum_N[i] + is_N[i];
  }
  LOG << "prefix_sum_N_last = " << prefix_sum_N[len - 1];
  // 构建 suspect region

  std::vector<std::pair<int, int>> gaps;

  // 这里可以进行整体的基数排序 - M
  // 目前改进的方法为，每个k-mer匹配向后5个的所有k-mer，作为gap
  int id = 0;
  std::vector<std::pair<int, int>> right_index_pair;
  for (int i = 0; i < len - 1; ++i) {
    if (lcp_array[i] < factor::kmer_size) {
      if (right_index.size() <= factor::min_anchor_size) {
        continue;
      }
      // 这里非常慢
      //  int start_position = *min_element(right_index.begin(),
      //  right_index.end()); int end_position =
      //  *max_element(right_index.begin(), right_index.end());
      int start_position = 0;
      int end_position = 0;
      int N_count = prefix_sum_N[end_position] - prefix_sum_N[start_position];
      if (N_count >
          /*factor::N_threshold * (end_position - start_position + 1)*/ 100) {
        continue;
      }
      for (int j : right_index) {
        right_index_pair.emplace_back(j, id);
      }
      id++;
      right_index.clear();
    }
    right_index.push_back(suffix_array[i + 1]);
  }

  if (right_index.size() > factor::min_anchor_size) {
    int start_position = right_index[0];
    int end_position = right_index.back();
    int N_count = prefix_sum_N[end_position + 1] - prefix_sum_N[start_position];
    if (!(N_count >
          factor::N_threshold * (end_position - start_position + 1))) {
      for (int j : right_index) {
        right_index_pair.emplace_back(j, id);
      }
      id++;
    }
  }
  LOG << "first id = " << id;
  radix_sort(
      right_index_pair, [](const auto &item) -> int { return item.first; },
      len);
  LOG << "out of radix sort";
  std::vector<std::vector<int>> right_index_2d(id);
  for (const auto &[index, id] : right_index_pair) {
    right_index_2d[id].push_back(index);
  }

  LOG << "right_index_2d.size() = " << right_index_2d.size();

  right_index = {suffix_array[0]};
  id = 0;
  int start_position = suffix_array[0];
  int end_position = suffix_array[0];
  for (int i = 0; i < len - 1; ++i) {
    if (lcp_array[i] < factor::kmer_size) {
      if (right_index.size() <= factor::min_anchor_size) {
        continue;
      }
      // 这里非常慢
      assert(start_position <= end_position);
      int N_count = prefix_sum_N[end_position] - prefix_sum_N[start_position];
      if (N_count > /*factor::N_threshold * (end_position - start_position + 1)*/ 100) {
        right_index.clear();
        start_position = len - 1;
        end_position = 0;
        continue;
      }
      const auto &sorted_right_index = right_index_2d[id];
      //   LOG << "sorted_right_index.size() = " << sorted_right_index.size();
      for (int j = 0; j < sorted_right_index.size(); ++j) {
        for (int k = 1; k <= 10 && j + k < sorted_right_index.size(); ++k) {
          int gap = sorted_right_index[j + k] - sorted_right_index[j];
          if (gap >= factor::lower_bound && gap <= factor::upper_bound) {
            gap_values.push_back(gap);
            gaps.emplace_back(gap, sorted_right_index[j]);
          }
        }
      }
      id++;
      right_index.clear();
      start_position = len - 1;
      end_position = 0;
    }
    right_index.push_back(suffix_array[i + 1]);
    start_position = min(start_position, suffix_array[i + 1]);
    end_position = max(end_position, suffix_array[i + 1]);
  }

  if (right_index.size() > factor::min_anchor_size) {
    int start_position = *min_element(right_index.begin(), right_index.end());
    int end_position = *max_element(right_index.begin(), right_index.end());
    int N_count = prefix_sum_N[end_position + 1] - prefix_sum_N[start_position];
    if (!(N_count > factor::N_threshold * (end_position - start_position + 1))) {
      auto sorted_right_index = right_index_2d[id];
      assert(is_sorted(
          sorted_right_index.begin(), sorted_right_index.end(),
          [](const auto &a, const auto &b) -> bool { return a < b; }));
      for (int j = 0; j < sorted_right_index.size(); ++j) {
        for (int k = 1; k <= 10 && j + k < sorted_right_index.size(); ++k) {
          int gap = sorted_right_index[j + k] - sorted_right_index[j];
          //   LOG << "gap = " << gap;
          if (gap >= factor::lower_bound && gap <= factor::upper_bound) {
            gap_values.push_back(gap);
            gaps.emplace_back(gap, sorted_right_index[j]);
          }
        }
      }
    }
  }

  LOG << "gaps.size() = " << gaps.size();

  std::function<int(const std::pair<int, std::vector<int>> &)> get_index =
      [](const std::pair<int, std::vector<int>> &item) -> int {
    return item.first;
  };

  radix_sort(gaps, [](const auto &item) -> int { return item.second; }, len);
  radix_sort(gap_values, [](const auto &item) -> int { return item; }, len);
  gap_values.erase(std::unique(gap_values.begin(), gap_values.end()), gap_values.end());

  std::vector<std::vector<pair<int, int>>> bucket_of_index(gap_values.size());

  // 对gaps处理一下，挑出一些质数，然后处理
  // 这里对gap排序了，所以之后bucket_of_index内部不要排序
  for (int i = 0; i < gaps.size(); ++i) {
    int start_position = gaps[i].second;
    int end_position = start_position + gaps[i].first;
    int gap = gaps[i].first;
    int gap_index =
        std::lower_bound(prime_number, prime_number + 78, gap) - prime_number;
    int dis1 = std::abs(gap - prime_number[gap_index]);
    int dis2 = std::abs(gap - prime_number[gap_index - 1]);
    if (dis1 <= 50) {
      bucket_of_index[gap_index].emplace_back(start_position, gap);
    }
    if (dis2 <= 50) {
      bucket_of_index[gap_index - 1].emplace_back(start_position, gap);
    }
  }

  LOG << "bucket_of_index.size() = "
      << bucket_of_index.size();
  // 校验一下
  for (auto &vec : bucket_of_index) {
    assert(std::is_sorted(vec.begin(), vec.end(), [](const auto &a, const auto &b) -> bool { return a.first < b.first; }));
  }

  LOG << "gap_values.size() = " << gap_values.size();
  std::vector<std::tuple<int, int, int>> raw_estimate_unit_region;

  // 排好序后，计算原始的
  for (int i = 0; i < bucket_of_index.size(); ++i) {
    int len_of_bucket = bucket_of_index[i].size();
    int k = 0;
    int num_of_index = 1;
    for (int j = 0; j < len_of_bucket; ++j) {
      int start_position = bucket_of_index[i][j].first;
      int end_position = start_position + bucket_of_index[i][j].second;
      int gap = bucket_of_index[i][j].second;
      num_of_index--;
      while (k < len_of_bucket && bucket_of_index[i][k].first < end_position) {
        k++;
        num_of_index++;
      }
      if (num_of_index >= gap * factor::kmer_rate) {
        //  LOG << "start_position = " << start_position << ", gap = " << gap
        //      << ", num_of_index = " << num_of_index;
        raw_estimate_unit_region.emplace_back(start_position, gap, i);
      }
    }
  }

  LOG << "raw_estimate_unit_region.size() = "
      << raw_estimate_unit_region.size();

  // 按照gap从小到大排序去做
  radix_sort(
      raw_estimate_unit_region,
      [](const auto &item) -> int { return std::get<1>(item); }, len);

  LOG << "radix_sort done";

  int max_r_boundary = -1;
  int max_l_boundary = -1;

  std::unordered_map<int, std::pair<int, int>> max_covered_region;

  // 树状数组
  MaxFenwickTree ft(len);
  auto already_covered = [&max_l_boundary, &max_r_boundary, &ft,
                          &len](int start_position, int end_position,
                                int unit_size) {
    int l_max = ft.queryPrefixMax(min(start_position + unit_size, len));
    //  LOG << "l_max = " << l_max << ", start_position = " << start_position
    //      << ", end_position = " << end_position << ", unit_size = " <<
    //      unit_size;
    return l_max >= end_position;
  };

  char *query = new char[len];
  for (int i = 0; i < len; ++i) {
    query[i] = safeNumberToDnaChar(seq[i]);
  }

  int cnt = 0;
  std::vector<std::tuple<int, int, int, std::vector<int>>> raw_estimate_interval;
  // 现在是按照unit_size排序
  // 找到当前gap对应的bucket，假设上次区间分割点为[s, e]， 下次查找两个匹配的k-mer (s1, e1)，其中(s < s1 < e < e1)，且s1离s和e1离e都尽量近
  // 再做alignment，向右延伸
  for (const auto &[start_position, unit_size, index] : raw_estimate_unit_region) {
    int end_position = start_position + unit_size;
    if (already_covered(start_position, end_position, unit_size)) {
      continue;
    }
    LOG << "try_to_extend: start_position = " << start_position
        << ", unit_size = " << unit_size << ", index = " << index;
    cnt++;
    // bucket是一个vector，存储了gap=gap_index[index]的所有下标
    const auto &bucket = bucket_of_index[index];
    // 向左延伸
    // 参考了tidehunter的做法
    int raw_left_boundary = start_position;
    int raw_right_boundary = end_position;
    std::vector<int> anchors;

    int S = start_position;
    int E = end_position;
    int n_cigar;
    uint32_t *cigar;
    while (E < len) {
      auto iter = lower_bound(bucket.begin(), bucket.end(), std::make_pair(E, 0));
      if (iter == bucket.end()) {
        break;
      }
      LOG << "E = " << E;
      if (iter == bucket.begin()) {
        LOG << "Error: iter == bucket.begin()";
        break;
      }
      auto iter_prev = prev(iter);

      int s1 = iter_prev->first;
      int e1 = s1 + iter_prev->second;
      int s2 = iter->first;
      int e2 = s2 + iter->second;
      int qlen = e2 - e1 + factor::kmer_size;
      int tlen = s2 - s1 + factor::kmer_size;
      if (qlen > 10000 || tlen > 10000) {
        break;
      }
      int iden_n = alignment::ksw2_global_with_cigar(
          seq + e1 - factor::kmer_size + 1, qlen,
          seq + s1 - factor::kmer_size + 1, tlen, &n_cigar, &cigar);
      LOG << "iden_n = " << iden_n << ' ' << s2 - s1 + factor::kmer_size << ' '
          << e2 - e1 + factor::kmer_size;
      if (iden_n >=
          min(s2 - s1 + factor::kmer_size, e2 - e1 + factor::kmer_size) * 0.9) {
        S = E;
        // 这里待修改
        E = e2 + alignment::extend_right_boundary(
                     seq + e1, e2 - e1 + factor::kmer_size, seq + s1,
                     s2 - s1 + factor::kmer_size);
        anchors.push_back(S);
        LOG << "S = " << S << ' ' << "E = " << E;
      } else {
        // 插入分割位置
        break;
      }
    }
    raw_right_boundary = E;

    S = start_position;
    E = end_position;
    while (S > 0) {
      auto iter = lower_bound(bucket.begin(), bucket.end(), std::make_pair(S, 0));
      if (iter == bucket.end()) {
        break;
      }
      if (iter == bucket.begin()) {
        LOG << "Error: iter == bucket.begin()";
        break;
      }
      auto iter_prev = prev(iter);
      int s1 = iter_prev->first;
      int e1 = s1 + iter_prev->second;
      int s2 = iter->first;
      int e2 = s2 + iter->second;
      int tlen = s2 - s1 + factor::kmer_size;
      int qlen = e2 - e1 + factor::kmer_size;
      if (qlen > 10000 || tlen > 10000) {
        break;
      }
      int iden_n = alignment::ksw2_global_with_cigar(
          seq + e1 - factor::kmer_size + 1, qlen,
          seq + s1 - factor::kmer_size + 1, tlen, &n_cigar, &cigar);
      if (iden_n >= min(s2 - s1 + factor::kmer_size, e2 - e1 + factor::kmer_size) * 0.9) {
        E = S;
        // 这里待修改
        S = s1 - alignment::extend_left_boundary(seq + s1,  s2 - s1 + factor::kmer_size, seq + e1, e2 - e1 + factor::kmer_size);
        anchors.push_back(E);
      } else {
        // 插入分割位置
        break;
      }
    }
    raw_left_boundary = S;
    LOG << "raw_left_boundary = " << raw_left_boundary
        << ", raw_right_boundary = " << raw_right_boundary;
    anchors.push_back(raw_left_boundary);
    anchors.push_back(raw_right_boundary);
    sort(anchors.begin(), anchors.end());
    anchors.erase(unique(anchors.begin(), anchors.end()), anchors.end());
    if (anchors.size() < 5) {
      continue;
    }

    // 树状数组更新
    ft.updateMax(raw_left_boundary, raw_right_boundary);
    raw_estimate_interval.emplace_back(unit_size, raw_left_boundary,
                                       raw_right_boundary, std::move(anchors));
  }

  LOG << "cnt = " << cnt;
  LOG << "raw_estimate_interval.size() = " << raw_estimate_interval.size();

  auto alignment_engine = spoa::AlignmentEngine::Create(
      spoa::AlignmentType::kNW, 3, -5, -3); // linear gaps
  spoa::Graph graph{};

  LOG << "start spoa alignment111";

  ofs << "fasta:" << std::endl;

  for (auto [unit_size, left_boundary, right_boundary, anchors] :
       raw_estimate_interval) {
    int seg_count =
        (right_boundary - left_boundary + unit_size - 1) / unit_size;

    int anchor_size = static_cast<int>(anchors.size());
    
    std::vector<std::string> sequences;
    // 1804 5099 8429 11753 15077 18408 21777 25101 28432
    for (int i = 0; i < anchor_size - 1; ++i) {
      std::string cur_seq = "";
      for (int j = anchors[i]; j < anchors[i + 1]; ++j) {
        cur_seq += safeNumberToDnaChar(seq[j]);
      }
      sequences.push_back(cur_seq);
    }

    LOG << "start spoa alignment";

    for (auto &sequence : sequences) {
      auto alignment = alignment_engine->Align(sequence, graph);
      graph.AddAlignment(alignment, sequence);
    }

    LOG << "end spoa alignment";

    auto consensus = graph.GenerateConsensus();
    int consensus_len = static_cast<int>(consensus.size());
    uint8_t *cons = alloc_uint8_t(consensus);

    double all_match = 0;
    double all_total = 0;
    int tot_seq = sequences.size();
    LOG << "tot_seq = " << tot_seq;
    uint8_t *seq0 = alloc_uint8_t(sequences[tot_seq / 2]);
    int seq0_len = static_cast<int>(sequences[tot_seq / 2].size());
    for (int i = 0; i < tot_seq; ++i) {
      auto &sequence = sequences[i];
      uint8_t *seq = alloc_uint8_t(sequence);
      int seq_len = static_cast<int>(sequence.size());
      auto [match, total] =
          alignment::alignment(cons, consensus_len, seq, seq_len);
      //  auto [match, total] =
      //      alignment::alignment(cons, consensus_len, seq, seq_len);
      all_match += match;
      all_total += total;
      LOG << "match = " << match << ", total = " << total;
      delete[] seq;
    }
    delete[] seq0;

    ofs << "[" << left_boundary << "," << right_boundary << "]" << " " << unit_size << " "
        << all_match / all_total << ' ' << consensus << std::endl;

    for (int i = 0; i < anchors.size(); ++i) {
      ofs << anchors[i] << ' ';
    }
    ofs << endl;
  }
}