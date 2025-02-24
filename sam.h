#pragma once

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

class SuffixAutomation {
public:
  SuffixAutomation(const std::string &s);
  void new_node();
  void add(char c, int pos);
  std::vector<int> toposort();
  int size();
  void solve();
  void get_right_index(int bound);
  bool judge_len_range(int x, int bound);
  std::vector<int> get_right_index_by_node(int x);

  std::vector<std::vector<int>> right_index;

private:
  std::vector<int> len = {0};
  std::vector<int> link = {-1};
  std::vector<std::map<char, int>> nxt = {{}};
  std::map<int, int> node_to_pos;
  int sz = 1;
  int last = 0;
};