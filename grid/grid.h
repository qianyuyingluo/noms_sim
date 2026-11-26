#pragma once

#include <vector>

// cuts 是网格线必须经过的坐标，count 是要分出的格子数。
// 返回从小到大排列的 count + 1 条网格线，单位和 cuts 相同。
std::vector<double> divideAxis(std::vector<double> cuts, int count);
