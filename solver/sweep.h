#pragma once

#include "dc_solver.h"

struct TransferPoint {
    double vgs;
    double drainCurrentA;
};

void saveOutputSweep(const DeviceMesh& mesh, const MOSFET& device, const DcResult& equilibrium);
std::vector<TransferPoint> saveTransferSweep(const DeviceMesh& mesh, const MOSFET& device,
                                              const DcResult& equilibrium, double vds,
                                              const char* path);
