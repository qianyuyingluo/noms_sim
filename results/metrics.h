#pragma once

#include "../solver/sweep.h"

struct MetricSettings {
    double highVds;         // DIBL 的高漏压，也用于提取 gm、ro 和 Cgg
    double differenceStep; // 中心差分两侧的电压变化量
};

// 输出指标表、高漏压转移数据和 gm/Id 曲线数据。
void saveMetrics(const DeviceMesh& mesh, const MOSFET& device,
                 const DcResult& equilibrium, const DcResult& off, const DcResult& on,
                 const std::vector<TransferPoint>& transfer, const MetricSettings& settings);
