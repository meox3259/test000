#include <functional>
#include <vector>

template <typename T>
void radix_sort(std::vector<T> &arr,
                const std::function<int(const T &)> &get_index, int max_digit) {
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
  assert(std::is_sorted(arr.begin(), arr.end()));
}

template <typename T, typename F>
void radix_sort(std::vector<T> &arr, F &&get_index, int max_digit) {
  radix_sort<T>(arr, std::forward<F>(get_index), max_digit);
}