#include "device/device.h"
#include "results/fields.h"
#include "results/metrics.h"
#include "solver/poisson.h"
#include "solver/sweep.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

MOSFET exampleDevice() {
    MOSFET device{};

    // 示例剖面尺寸。栅氧厚度用 nm，其余长度用 um。
    device.SiL = 4.0;
    device.SiD = 1.0;
    device.Lg = 2.0;
    device.Ls = 1.0;
    device.Ld = 1.0;
    device.xjs = 0.2;
    device.xjd = 0.2;
    device.tox = 10.0;

    // 掺杂浓度，单位为 cm^-3。
    device.NA = 1e16;
    device.NDs = 1e18;
    device.NDd = 1e18;

    device.epsSi = 11.7;
    device.epsOx = 3.9;
    device.W = 1.0;
    device.ni = 1e10;
    device.mu_n = 1000;
    device.mu_p = 450;
    device.T = 300;
    device.useVFB = true;
    device.VFB = -0.25;

    device.Vgs_off = 0.3;
    device.Vds_off = 0.1;
    device.Vgs_on = 1.4;
    device.Vds_on = 0.1;
    device.Vgs_list = {0.6, 0.9, 1.2, 1.4};
    device.VdsFrom = 0.0;
    device.VdsTo = 1.0;
    device.VdsStep = 0.1;
    device.VgsFrom = 0.3;
    device.VgsTo = 1.4;
    device.VgsStep = 0.05;
    device.Vds_transfer = 0.1;
    return device;
}

int saveElectrostaticExample(const DeviceMesh& mesh, const MOSFET& device) {
    // 静电示例：栅极为 1 V，其余接触为 0 V，体电荷密度为零。
    PoissonBoundary boundary;
    boundary.gate = 1.0;
    const std::vector<double> rho(mesh.cells.size(), 0.0);
    const auto potential = solvePoisson(mesh, rho, device.epsSi, device.epsOx, boundary);

    std::ofstream csv("output/potential.csv");
    csv << "row,col,x_left_um,x_right_um,y_top_um,y_bottom_um,x_um,y_um,potential_V\n";
    csv << std::setprecision(16);

    // 每行保存一个格子的边界、中心和电势。写 CSV 时坐标换回微米。
    const std::size_t columns = mesh.x.size() - 1;
    for (std::size_t row = 0; row + 1 < mesh.y.size(); ++row) {
        for (std::size_t col = 0; col < columns; ++col) {
            const double left = mesh.x[col] * 1e6;
            const double right = mesh.x[col + 1] * 1e6;
            const double top = mesh.y[row] * 1e6;
            const double bottom = mesh.y[row + 1] * 1e6;
            csv << row << ',' << col << ',' << left << ',' << right << ','
                << top << ',' << bottom << ',' << (left + right) / 2.0 << ','
                << (top + bottom) / 2.0 << ',' << potential[row * columns + col] << '\n';
        }
    }
    csv.close();
    if (!csv) {
        std::cerr << "Could not write output/potential.csv\n";
        return 1;
    }
    std::cout << "Saved output/potential.csv (" << mesh.cells.size() << " cells)\n";
    return 0;
}

int main() {
    try {
        MOSFET device = exampleDevice();

        // ===== 网格密度：做网格加密对比时，在这里修改 =====
        device.Nx = 40;                   // 从源极到漏极方向的总格数
        device.NySi = 40;                 // 硅区深度方向的总格数
        device.NyOx = 5;                  // 栅氧层厚度方向的总格数
        device.surfaceMeshDepthUm = 0.012; // 在距硅表面 12 nm 以内加密，设为 0 关闭

        const DeviceMesh mesh = buildMesh(device);
        std::filesystem::create_directories("output");
        if (saveElectrostaticExample(mesh, device) != 0) {
            return 1;
        }

        const DcResult equilibrium = solveDc(mesh, device, 0, 0);
        const DcResult off = solveDc(mesh, device, device.Vgs_off, device.Vds_off, &equilibrium);
        const DcResult on = solveDc(mesh, device, device.Vgs_on, device.Vds_on, &off);
        saveFields("output/fields_off.csv", mesh, off);
        saveFields("output/fields_on.csv", mesh, on);
        std::cout << "Saved fields_off.csv and fields_on.csv\n";

        saveOutputSweep(mesh, device, equilibrium);
        const auto transfer = saveTransferSweep(mesh, device, equilibrium, device.Vds_transfer,
                                                 "output/transfer_iv.csv");
        std::cout << "Saved output_iv.csv and transfer_iv.csv\n";
        const MetricSettings metrics{device.VdsTo, 0.005};
        saveMetrics(mesh, device, equilibrium, off, on, transfer, metrics);
        std::cout << "Saved metrics.csv, gm_id.csv and transfer_high_vds.csv\n";
        std::cout << "Drain current: off " << off.drainCurrentA << " A, on " << on.drainCurrentA << " A\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
