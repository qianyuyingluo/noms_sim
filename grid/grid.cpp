#include "grid.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// 在必须经过的边界之间插入网格线，一共分成 count 格。
std::vector<double> divideAxis(std::vector<double> cuts, int count) {
    for (double cut : cuts) {
        if (!std::isfinite(cut)) {
            throw std::invalid_argument("Mesh coordinates must be finite.");
        }
    }
    if (cuts.size() < 2) {
        throw std::invalid_argument("At least two boundary coordinates are required.");
    }
    std::sort(cuts.begin(), cuts.end());
    const double span = cuts.back() - cuts.front();
    if (!std::isfinite(span) || span <= 0.0) {
        throw std::invalid_argument("Mesh boundaries must span a finite positive length.");
    }
    const double tolerance = 16.0 * std::numeric_limits<double>::epsilon() * span;
    std::vector<double> boundaries{cuts.front()};
    for (double cut : cuts) {
        if (cut - boundaries.back() > tolerance) {
            boundaries.push_back(cut);
        }
    }

    const int sections = static_cast<int>(boundaries.size()) - 1;
    if (count < sections || count <= 0) {
        throw std::invalid_argument("Mesh count is too small for the required boundaries.");
    }

    // 每个区段至少分到一格，其余格子按区段长度分配。
    const int extra = count - sections;
    int assignedExtra = 0;
    std::vector<double> lines{boundaries.front()};
    for (int i = 0; i < sections; ++i) {
        const double left = boundaries[i];
        const double right = boundaries[i + 1];
        // 用累计长度分配格子，避免逐段四舍五入后总数不对。
        const double fraction = (right - boundaries.front()) / span;
        int totalExtra = static_cast<int>(std::lround(extra * fraction));
        if (i == sections - 1) {
            totalExtra = extra;
        }
        const int cells = 1 + totalExtra - assignedExtra;
        assignedExtra = totalExtra;
        const double step = (right - left) / cells;
        for (int k = 1; k < cells; ++k) {
            lines.push_back(left + k * step);
        }
        lines.push_back(right);
    }
    return lines;
}
