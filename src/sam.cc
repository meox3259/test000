#include "sam.h"
#include "utils.h"

#include <cassert>
#include <queue>

SuffixAutomation::SuffixAutomation(const std::string &s) {
  for (int i = 0; i < s.size(); ++i) {
    add(s[i], i);
  }
}

void SuffixAutomation::new_node() {
  sz++;
  len.resize(sz);
  link.resize(sz);
  nxt.resize(sz);
}

void SuffixAutomation::add(char c, int pos) {
  int now = sz;
  new_node();
  len[now] = len[last] + 1;
  int p;
  for (p = last; p >= 0 && nxt[p].find(c) == nxt[p].end(); p = link[p]) {
    nxt[p][c] = now;
  }
  if (p < 0) {
    link[now] = 0;
  } else {
    int q = nxt[p][c];
    if (len[q] == len[p] + 1) {
      link[now] = q;
    } else {
      int clone = sz;
      new_node();
      len[clone] = len[p] + 1;
      link[clone] = link[q];
      nxt[clone] = nxt[q];
      for (; p >= 0 && nxt[p][c] == q; p = link[p]) {
        nxt[p][c] = clone;
      }
      link[q] = link[now] = clone;
    }
  }
  node_to_pos[now] = pos;
  last = now;
}

std::vector<int> SuffixAutomation::toposort() {
  std::vector<int> c(sz);
  for (int i = 0; i < sz; ++i) {
    c[len[i]]++;
  }
  for (int i = 0; i + 1 < sz; ++i) {
    c[i + 1] += c[i];
  }
  std::vector<int> topo(sz);
  for (int i = 0; i < sz; ++i) {
    topo[--c[len[i]]] = i;
  }
  return topo;
}

void SuffixAutomation::get_right_index(int bound) {
  right_index.resize(sz);
  auto topo = toposort();

  std::vector<std::vector<int>> adj(sz);
  for (int i = 0; i < sz; ++i) {
    if (link[i] >= 0) {
      adj[link[i]].push_back(i);
    }
  }

  auto dfs = [&adj, bound, this](auto &&dfs, int x) -> void {
    if (judge_len_range(x, 30)) {
      std::queue<int> q;
      q.push(x);
      while (!q.empty()) {
        int now = q.front();
        if (node_to_pos.find(now) != node_to_pos.end()) {
          right_index[x].push_back(node_to_pos[now]);
        }
        q.pop();
        for (int nxt : adj[now]) {
          q.push(nxt);
        }
      }
      return;
    }
    for (int y : adj[x]) {
      dfs(dfs, y);
    }
  };

  dfs(dfs, 0);

  for (int i = 0; i < sz; ++i) {
    sort(right_index[i].begin(), right_index[i].end());
  }
}

int SuffixAutomation::size() { return sz; }

bool SuffixAutomation::judge_len_range(int x, int bound) {
  if (link[x] < 0) {
    return 0 >= bound;
  }
  return len[link[x]] + 1 >= bound;
}

std::vector<int> SuffixAutomation::get_right_index_by_node(int x) {
  if (x >= sz || x < 0) {
    std::cerr << "Error: node " << x << " not in range" << std::endl;
    return {};
  }
  return right_index[x];
}