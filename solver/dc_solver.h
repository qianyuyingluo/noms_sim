#pragma once

#include "../device/device.h"

struct DcResult {
    double vgs = 0.0;
    double vds = 0.0;
    std::vector<double> potentialV;
    std::vector<double> electronsPerM3;
    std::vector<double> holesPerM3;
    std::vector<double> fieldX;
    std::vector<double> fieldY;
    std::vector<double> currentX;
    std::vector<double> currentY;
    double drainCurrentA = 0.0;
    double sourceCurrentA = 0.0;
    double bodyCurrentA = 0.0;
    double gateChargeC = 0.0;

    // 下一个偏压点可沿用上个点的解作为初值。
    std::vector<double> solution;
};

DcResult solveDc(const DeviceMesh& mesh, const MOSFET& device, double vgs, double vds,
                 const DcResult* previous = nullptr);
