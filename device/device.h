#pragma once

#include "../NMOS.h"

#include <cstddef>
#include <vector>

enum class Material { Outside, Silicon, Oxide };
enum class Contact { None, Source, Drain, Gate, Body };

struct MeshCell {
    Material material;
    double netDopingPerM3;
    Contact top;
    Contact bottom;
};

struct DeviceMesh {
    // 网格线从小到大排列，坐标单位为 m。纵坐标向硅内部增大。
    std::vector<double> x;
    std::vector<double> y;
    // cells[j * (x.size() - 1) + i] 对应第 i 列、第 j 行。
    std::vector<MeshCell> cells;

    const MeshCell& cell(std::size_t i, std::size_t j) const;
};

DeviceMesh buildMesh(const MOSFET& device);
