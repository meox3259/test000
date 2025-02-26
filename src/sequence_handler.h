#pragma once

#include <memory>
#include <string>
#include <vector>

#include "utils.h"

bool calculate_sparsity(const std::vector<int> &index);
std::pair<int, int> get_overlap(const std::vector<int> &l,
                                const std::vector<int> &r);
void solve(uint8_t *s, int len, const Param &param);