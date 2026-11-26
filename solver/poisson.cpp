#include "poisson.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

struct Neighbor {
    std::size_t index = 0;
    double weight = 0.0;
};

struct CellEquation {
    std::array<Neighbor, 4> neighbors{};
    double totalWeight = 0.0;
    double rhs = 0.0;
};

// 相邻两格的半格厚度相当于串联的两层介质。
void addNeighbor(CellEquation& equation, int side, std::size_t index,
                 double faceLength, double ownHalf, double otherHalf,
                 double ownEps, double otherEps) {
    if (otherEps == 0.0) {
        return; // 域外边界的法向电通量为零
    }
    const double weight = faceLength / (ownHalf / ownEps + otherHalf / otherEps);
    equation.neighbors[side] = {index, weight};
    equation.totalWeight += weight;
}

void addContact(CellEquation& equation, Contact contact, double weight,
                const PoissonBoundary& boundary, bool& anchored) {
    double voltage = 0.0;
    switch (contact) {
        case Contact::None: return;
        case Contact::Source: voltage = boundary.source; break;
        case Contact::Drain: voltage = boundary.drain; break;
        case Contact::Gate: voltage = boundary.gate; break;
        case Contact::Body: voltage = boundary.body; break;
    }
    equation.totalWeight += weight;
    equation.rhs += weight * voltage;
    anchored = true;
}

std::vector<double> cellWidths(const std::vector<double>& lines) {
    if (lines.size() < 2) {
        throw std::invalid_argument("Each axis needs at least two mesh lines.");
    }
    std::vector<double> widths;
    for (std::size_t i = 1; i < lines.size(); ++i) {
        const double width = lines[i] - lines[i - 1];
        if (!std::isfinite(width) || width <= 0.0) {
            throw std::invalid_argument("Mesh lines must be finite and strictly increasing.");
        }
        widths.push_back(width);
    }
    return widths;
}

} // namespace

std::vector<double> solvePoisson(const DeviceMesh& mesh,
                                 const std::vector<double>& rho,
                                 double epsSi, double epsOx,
                                 const PoissonBoundary& boundary) {
    const auto dx = cellWidths(mesh.x);
    const auto dy = cellWidths(mesh.y);
    const std::size_t columns = dx.size();
    const std::size_t rows = dy.size();
    const std::size_t count = mesh.cells.size();
    if (count != columns * rows || rho.size() != count) {
        throw std::invalid_argument("Mesh and charge density sizes must match.");
    }
    if (!std::isfinite(epsSi) || !std::isfinite(epsOx) || epsSi <= 0.0 || epsOx <= 0.0 ||
        !std::isfinite(boundary.source) || !std::isfinite(boundary.drain) ||
        !std::isfinite(boundary.gate) || !std::isfinite(boundary.body)) {
        throw std::invalid_argument("Permittivities must be positive and boundary potentials must be finite.");
    }

    constexpr double eps0 = 8.8541878128e-12; // F/m
    std::vector<double> epsilon(count, 0.0);
    for (std::size_t k = 0; k < count; ++k) {
        if (mesh.cells[k].material == Material::Silicon) {
            epsilon[k] = eps0 * epsSi;
        } else if (mesh.cells[k].material == Material::Oxide) {
            epsilon[k] = eps0 * epsOx;
        }
    }

    // 先给每个格子写下通量平衡式，再迭代更新电势。
    std::vector<CellEquation> equations(count);
    std::vector<double> phi(count, std::numeric_limits<double>::quiet_NaN());
    std::vector<bool> reached(count, false);
    std::vector<std::size_t> connectedToContact;
    std::size_t activeCount = 0;
    for (std::size_t j = 0; j < rows; ++j) {
        for (std::size_t i = 0; i < columns; ++i) {
            const std::size_t k = j * columns + i;
            const double eps = epsilon[k];
            if (eps == 0.0) {
                continue;
            }
            if (!std::isfinite(rho[k])) {
                throw std::invalid_argument("Charge density must be finite inside the domain.");
            }
            ++activeCount;
            phi[k] = 0.0;
            auto& eq = equations[k];
            eq.rhs = rho[k] * dx[i] * dy[j];
            if (i > 0) {
                addNeighbor(eq, 0, k - 1, dy[j], dx[i] / 2, dx[i - 1] / 2, eps, epsilon[k - 1]);
            }
            if (i + 1 < columns) {
                addNeighbor(eq, 1, k + 1, dy[j], dx[i] / 2, dx[i + 1] / 2, eps, epsilon[k + 1]);
            }
            bool anchored = false;
            const auto& cell = mesh.cells[k];
            if (j > 0 && epsilon[k - columns] > 0.0) {
                addNeighbor(eq, 2, k - columns, dx[i], dy[j] / 2, dy[j - 1] / 2, eps, epsilon[k - columns]);
            } else {
                addContact(eq, cell.top, dx[i] * eps / (dy[j] / 2), boundary, anchored);
            }
            if (j + 1 < rows && epsilon[k + columns] > 0.0) {
                addNeighbor(eq, 3, k + columns, dx[i], dy[j] / 2, dy[j + 1] / 2, eps, epsilon[k + columns]);
            } else {
                addContact(eq, cell.bottom, dx[i] * eps / (dy[j] / 2), boundary, anchored);
            }
            if (anchored) {
                reached[k] = true;
                connectedToContact.push_back(k);
            }
        }
    }

    // 每片连通区域都需要接到至少一个已知电势的边界。
    for (std::size_t i = 0; i < connectedToContact.size(); ++i) {
        for (const auto& neighbor : equations[connectedToContact[i]].neighbors) {
            if (neighbor.weight > 0.0 && !reached[neighbor.index]) {
                reached[neighbor.index] = true;
                connectedToContact.push_back(neighbor.index);
            }
        }
    }
    if (activeCount == 0 || connectedToContact.size() != activeCount) {
        throw std::invalid_argument("Each active region needs a fixed-potential contact.");
    }

    constexpr int maxSweeps = 100000;
    constexpr double tolerance = 1e-10; // V
    for (int sweep = 0; sweep < maxSweeps; ++sweep) {
        double largestChange = 0.0;
        for (std::size_t k = 0; k < count; ++k) {
            if (epsilon[k] == 0.0) {
                continue;
            }
            const auto& eq = equations[k];
            double numerator = eq.rhs;
            for (const auto& neighbor : eq.neighbors) {
                if (neighbor.weight > 0.0) {
                    numerator += neighbor.weight * phi[neighbor.index];
                }
            }
            // Gauss–Seidel：算出新电势后立即供后面的格子使用。
            const double updated = numerator / eq.totalWeight;
            if (!std::isfinite(updated)) {
                throw std::runtime_error("Poisson iteration produced a non-finite potential.");
            }
            largestChange = std::max(largestChange, std::abs(updated - phi[k]));
            phi[k] = updated;
        }
        if (largestChange < tolerance) {
            return phi;
        }
    }
    throw std::runtime_error("Poisson solver did not converge.");
}
