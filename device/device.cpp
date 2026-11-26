#include "device.h"
#include "../grid/grid.h"

#include <cmath>
#include <stdexcept>

namespace {

bool positiveFinite(double value) {
    return std::isfinite(value) && value > 0.0;
}

void checkDevice(const MOSFET& device) {
    const double lengths[] = {device.SiL, device.SiD, device.Lg, device.Ls,
                              device.Ld, device.xjs, device.xjd, device.tox};
    for (double length : lengths) {
        if (!positiveFinite(length)) {
            throw std::invalid_argument("Device dimensions must be finite and positive.");
        }
    }
    if (device.Lg >= device.SiL || device.Ls + device.Ld >= device.SiL ||
        device.xjs >= device.SiD || device.xjd >= device.SiD) {
        throw std::invalid_argument("Gate, source, or drain does not fit in the silicon domain.");
    }
    if (!positiveFinite(device.NA) || !positiveFinite(device.NDs) ||
        !positiveFinite(device.NDd) || device.NDs <= device.NA || device.NDd <= device.NA) {
        throw std::invalid_argument("Source and drain donor doping must exceed the positive body acceptor doping.");
    }
}

// 根据一个格子的中心坐标判定材料、净掺杂和电极接触。
MeshCell describeCell(const MOSFET& device, double x, double y,
                      std::size_t row, std::size_t lastRow) {
    const auto surfaceRow = static_cast<std::size_t>(device.NyOx);
    const bool underGate = x >= device.gateStartUm() && x < device.gateEndUm();
    MeshCell cell{Material::Outside, 0.0, Contact::None, Contact::None};
    if (row < surfaceRow) {
        if (underGate) {
            cell.material = Material::Oxide;
            if (row == 0) {
                cell.top = Contact::Gate;
            }
        }
        return cell;
    }

    cell.material = Material::Silicon;
    double netDoping = -device.NA;
    const bool exposedSurface = row == surfaceRow && !underGate;
    if (x < device.Ls && y < device.xjs) {
        netDoping += device.NDs;
        if (exposedSurface) {
            cell.top = Contact::Source;
        }
    } else if (x >= device.SiL - device.Ld && y < device.xjd) {
        netDoping += device.NDd;
        if (exposedSurface) {
            cell.top = Contact::Drain;
        }
    }
    cell.netDopingPerM3 = netDoping * 1e6; // cm⁻³ 转 m⁻³
    if (row == lastRow) {
        cell.bottom = Contact::Body;
    }
    return cell;
}

} // namespace

DeviceMesh buildMesh(const MOSFET& device) {
    checkDevice(device);
    constexpr double umToM = 1e-6;
    constexpr double nmToM = 1e-9;

    // 1. 先用微米画网格线，栅氧和硅区共用 y = 0。
    DeviceMesh mesh;
    mesh.x = divideAxis({0.0, device.Ls, device.SiL - device.Ld,
                         device.gateStartUm(), device.gateEndUm(), device.SiL}, device.Nx);
    mesh.y = divideAxis({-device.tox * nmToM / umToM, 0.0}, device.NyOx);
    std::vector<double> siliconCuts{0.0, device.xjs, device.xjd, device.SiD};
    if (device.surfaceMeshDepthUm > 0) {
        const double depth = device.surfaceMeshDepthUm;
        if (depth >= device.SiD) {
            throw std::invalid_argument("Surface refinement depth must be inside silicon.");
        }
        // 反型沟道靠近硅表面，先在浅处增加几条网格线。
        siliconCuts.push_back(depth / 8);
        siliconCuts.push_back(depth / 3);
        siliconCuts.push_back(depth);
    }
    const auto siliconY = divideAxis(siliconCuts, device.NySi);
    for (std::size_t j = 1; j < siliconY.size(); ++j) {
        mesh.y.push_back(siliconY[j]); // 跳过已经存过的 y = 0
    }

    // 2. 两条相邻网格线围出一个格子，用中心坐标给它贴标签。
    const std::size_t columns = mesh.x.size() - 1;
    const std::size_t rows = mesh.y.size() - 1;
    for (std::size_t j = 0; j < rows; ++j) {
        const double y = (mesh.y[j] + mesh.y[j + 1]) / 2.0;
        for (std::size_t i = 0; i < columns; ++i) {
            const double x = (mesh.x[i] + mesh.x[i + 1]) / 2.0;
            mesh.cells.push_back(describeCell(device, x, y, j, rows - 1));
        }
    }
    // 3. 把坐标换成米，供方程求解器使用。
    for (double& x : mesh.x) {
        x *= umToM;
    }
    for (double& y : mesh.y) {
        y *= umToM;
    }
    return mesh;
}

const MeshCell& DeviceMesh::cell(std::size_t i, std::size_t j) const {
    if (x.empty() || y.empty() || i >= x.size() - 1 || j >= y.size() - 1) {
        throw std::out_of_range("Mesh cell index is out of range.");
    }
    return cells.at(j * (x.size() - 1) + i);
}
