#pragma once

#include <cassert>
#include <functional>
#include <vector>

#include "logger.h"

template <typename T, typename F>
void radix_sort(std::vector<T> &arr, F &&get_index, int max_digit) {
  LOG << "max_digit = " << max_digit << " " << "arr.size() = " << arr.size();
  std::vector<int> cnt(max_digit + 1, 0);
  for (const auto &item : arr) {
    cnt[get_index(item)]++;
  }
  for (int i = 1; i <= max_digit; ++i) {
    cnt[i] += cnt[i - 1];
  }
  std::vector<T> output(arr.size());
  for (auto &item : arr) {
    output[--cnt[get_index(item)]] = std::move(item);
  }
  std::swap(arr, output);
  assert(std::is_sorted(arr.begin(), arr.end(),
                        [&get_index](const T &a, const T &b) {
                          return get_index(a) < get_index(b);
                        }));
}