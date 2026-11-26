#pragma once

#include "../device/device.h"

#include <vector>

// 四个接触边的电势边界值，单位 V。
struct PoissonBoundary {
    double source = 0.0;
    double drain = 0.0;
    double gate = 0.0;
    double body = 0.0;
};

// 求解 -div(epsilon * grad(phi)) = rho。
// rho 是总电荷密度，单位 C/m^3，与 mesh.cells 一一对应。
// epsSi 和 epsOx 是相对介电常数，可取自 NMOS.h。
// 接触边使用 boundary 中的电势，其余外边界使用零法向电通量。
// 返回每个格子中心的电势（V），域外格子的结果用 NaN 标记。
std::vector<double> solvePoisson(const DeviceMesh& mesh,
                                 const std::vector<double>& rho,
                                 double epsSi,
                                 double epsOx,
                                 const PoissonBoundary& boundary);
