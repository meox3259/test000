#pragma once

#include <memory>
#include <string>
#include <vector>

#include "utils.h"

bool calculate_sparsity(const std::vector<int> &index);
void solve(uint8_t *s, int len, const Param &param);