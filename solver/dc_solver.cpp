#include "dc_solver.h"
#include "transport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

struct Face {
    int neighbor = -1;
    Contact contact = Contact::None;
    double half = 0;
    double distance = 0;
    double length = 0;
    double electrostaticWeight = 0;
};

struct Cell {
    std::array<Face, 4> faces{}; // 左、右、上、下
    double doping = 0;
    double epsilon = 0;
    double area = 0;
    int color = 0;
    bool silicon = false;
};

struct Model {
    std::vector<Cell> cells;
    int columns = 0;
    double ni = 0;
    double ut = 0;
    double dn = 0;
    double dp = 0;
    double bodyPhi = 0;
    double flatBand = 0;
};

struct Sample {
    double phi;
    double electronLevel;
    double holeLevel;
    Carriers carriers;
};

Model makeModel(const DeviceMesh& mesh, const MOSFET& device) {
    Model model;
    model.columns = static_cast<int>(mesh.x.size()) - 1;
    const int rows = static_cast<int>(mesh.y.size()) - 1;
    if (model.columns <= 0 || rows <= 0 ||
        mesh.cells.size() != static_cast<std::size_t>(model.columns) * rows ||
        !(device.ni > 0 && device.T > 0 && device.mu_n > 0 && device.mu_p > 0 &&
          device.W > 0 && device.epsSi > 0 && device.epsOx > 0)) {
        throw std::invalid_argument("DC solver needs a valid mesh and positive material parameters.");
    }
    model.ni = device.ni * 1e6;
    model.ut = thermalVoltage(device.T);
    model.dn = device.mu_n * 1e-4 * model.ut;
    model.dp = device.mu_p * 1e-4 * model.ut;
    model.bodyPhi = std::asinh(-device.NA * 1e6 / (2 * model.ni));
    // 以硅的本征能级近似位于禁带中央来估算硅功函数。
    model.flatBand = device.useVFB ? device.VFB :
        device.phiM - (4.05 + 1.12 / 2 - model.ut * model.bodyPhi);

    model.cells.resize(mesh.cells.size());
    for (std::size_t k = 0; k < mesh.cells.size(); ++k) {
        Cell& cell = model.cells[k];
        const Material material = mesh.cells[k].material;
        cell.silicon = material == Material::Silicon;
        if (material != Material::Outside) {
            cell.epsilon = vacuumPermittivity * (cell.silicon ? device.epsSi : device.epsOx);
        }
        cell.doping = mesh.cells[k].netDopingPerM3;
    }

    for (int j = 0; j < rows; ++j) {
        const double dy = mesh.y[j + 1] - mesh.y[j];
        for (int i = 0; i < model.columns; ++i) {
            const double dx = mesh.x[i + 1] - mesh.x[i];
            const int k = j * model.columns + i;
            Cell& cell = model.cells[k];
            cell.area = dx * dy;
            cell.color = (i + 2 * j) % 5;
            if (cell.epsilon == 0) {
                continue;
            }
            for (int direction = 0; direction < 4; ++direction) {
                Face& face = cell.faces[direction];
                const int ii = i + (direction == 1) - (direction == 0);
                const int jj = j + (direction == 3) - (direction == 2);
                const bool horizontal = direction < 2;
                face.half = (horizontal ? dx : dy) / 2;
                face.length = horizontal ? dy : dx;
                if (ii >= 0 && ii < model.columns && jj >= 0 && jj < rows) {
                    const int other = jj * model.columns + ii;
                    const double otherEps = model.cells[other].epsilon;
                    if (otherEps > 0) {
                        const double otherHalf = horizontal ?
                            (mesh.x[ii + 1] - mesh.x[ii]) / 2 : (mesh.y[jj + 1] - mesh.y[jj]) / 2;
                        face.neighbor = other;
                        face.distance = face.half + otherHalf;
                        face.electrostaticWeight = face.length / (face.half / cell.epsilon + otherHalf / otherEps);
                        continue;
                    }
                }
                if (direction == 2) {
                    face.contact = mesh.cells[k].top;
                } else if (direction == 3) {
                    face.contact = mesh.cells[k].bottom;
                }
                if (face.contact != Contact::None) {
                    face.distance = face.half;
                    face.electrostaticWeight = face.length * cell.epsilon / face.half;
                }
            }
        }
    }
    return model;
}

Sample readSample(const Model& model, const std::vector<double>& values, int k) {
    const double phi = values[3 * k];
    const double electronLevel = values[3 * k + 1];
    const double holeLevel = values[3 * k + 2];
    Carriers carriers{};
    if (model.cells[k].silicon) {
        carriers = carrierDensity(model.ni, phi, electronLevel, holeLevel);
    }
    return {phi, electronLevel, holeLevel, carriers};
}

Sample contactSample(const Model& model, const Cell& cell, Contact contact, double vgs, double vds) {
    if (contact == Contact::Gate) {
        return {model.bodyPhi + (vgs - model.flatBand) / model.ut, 0, 0, {}};
    }
    const double applied = contact == Contact::Drain ? vds / model.ut : 0.0;
    const double phi = std::asinh(cell.doping / (2 * model.ni)) + applied;
    return {phi, applied, applied, carrierDensity(model.ni, phi, applied, applied)};
}

bool carriesCurrent(const Model& model, const Cell& cell, const Face& face) {
    if (!cell.silicon) {
        return false;
    }
    if (face.neighbor >= 0) {
        return model.cells[face.neighbor].silicon;
    }
    return face.contact == Contact::Source || face.contact == Contact::Drain ||
           face.contact == Contact::Body;
}

CarrierFlux fluxAcross(const Model& model, const Face& face, Sample a, Sample b) {
    return carrierFlux(b.phi - a.phi, a.carriers, b.carriers,
                       b.electronLevel - a.electronLevel, b.holeLevel - a.holeLevel,
                       model.dn * face.length / face.distance, model.dp * face.length / face.distance);
}

std::vector<double> residual(const Model& model, const std::vector<double>& values, double vgs, double vds) {
    std::vector<Sample> samples;
    samples.reserve(model.cells.size());
    for (int k = 0; k < static_cast<int>(model.cells.size()); ++k) {
        samples.push_back(readSample(model, values, k));
    }
    std::vector<double> result = values;
    for (std::size_t k = 0; k < model.cells.size(); ++k) {
        const Cell& cell = model.cells[k];
        if (cell.epsilon == 0) {
            continue;
        }
        const Sample a = samples[k];
        const double chargeScale = elementaryCharge * cell.area / model.ut;
        double poisson = -chargeScale * (a.carriers.holes - a.carriers.electrons + cell.doping);
        double poissonScale = chargeScale * (2 * model.ni + std::abs(cell.doping));
        double electron = 0, hole = 0, electronScale = 0, holeScale = 0;
        for (const Face& face : cell.faces) {
            if (face.electrostaticWeight == 0) {
                continue;
            }
            const Sample b = face.neighbor >= 0 ? samples[face.neighbor] :
                contactSample(model, cell, face.contact, vgs, vds);
            poisson += face.electrostaticWeight * (a.phi - b.phi);
            poissonScale += face.electrostaticWeight;
            if (carriesCurrent(model, cell, face)) {
                const auto flux = fluxAcross(model, face, a, b);
                electron += flux.electrons;
                hole += flux.holes;
                electronScale += flux.electronScale;
                holeScale += flux.holeScale;
            }
        }
        result[3 * k] = poisson / poissonScale;
        if (cell.silicon) {
            result[3 * k + 1] = electron / electronScale;
            result[3 * k + 2] = hole / holeScale;
        }
    }
    return result;
}

// 按格子编号后，每条方程只依赖本格和四邻格，矩阵只需存对角线附近。
class BandMatrix {
public:
    BandMatrix(int size, int halfWidth)
        : size(size), halfWidth(halfWidth),
          values(static_cast<std::size_t>(size) * (2 * halfWidth + 1), 0.0) {}

    double& at(int row, int col) {
        return values[static_cast<std::size_t>(row) * (2 * halfWidth + 1) + col - row + halfWidth];
    }

    bool solve(std::vector<double>& rhs) {
        for (int k = 0; k < size; ++k) {
            const double pivot = at(k, k);
            if (!std::isfinite(pivot) || std::abs(pivot) < 1e-14) {
                return false;
            }
            const int end = std::min(size, k + halfWidth + 1);
            for (int i = k + 1; i < end; ++i) {
                const double factor = at(i, k) / pivot;
                if (factor == 0) {
                    continue;
                }
                for (int j = k + 1; j < end; ++j) {
                    at(i, j) -= factor * at(k, j);
                }
                rhs[i] -= factor * rhs[k];
            }
        }
        for (int i = size - 1; i >= 0; --i) {
            const int end = std::min(size, i + halfWidth + 1);
            for (int j = i + 1; j < end; ++j) {
                rhs[i] -= at(i, j) * rhs[j];
            }
            rhs[i] /= at(i, i);
        }
        return true;
    }

private:
    int size;
    int halfWidth;
    std::vector<double> values;
};

double maxResidual(const std::vector<double>& values) {
    double largest = 0;
    for (double value : values) {
        if (!std::isfinite(value)) {
            return std::numeric_limits<double>::infinity();
        }
        largest = std::max(largest, std::abs(value));
    }
    return largest;
}

BandMatrix makeJacobian(const Model& model, const std::vector<double>& values, double vgs, double vds) {
    BandMatrix matrix(static_cast<int>(values.size()), 3 * model.columns + 2);
    constexpr double step = 1e-4;
    // 五种颜色使任意一个格子的四个邻居与本格颜色都不同。
    // 一次可以扰动同色格子，再用差商得到每条方程的导数。
    for (int color = 0; color < 5; ++color) {
        for (int variable = 0; variable < 3; ++variable) {
            auto plus = values;
            auto minus = values;
            for (std::size_t k = 0; k < model.cells.size(); ++k) {
                if (model.cells[k].color == color) {
                    plus[3 * k + variable] += step;
                    minus[3 * k + variable] -= step;
                }
            }
            const auto rp = residual(model, plus, vgs, vds);
            const auto rm = residual(model, minus, vgs, vds);
            for (int k = 0; k < static_cast<int>(model.cells.size()); ++k) {
                int columnCell = k;
                if (model.cells[k].color != color) {
                    columnCell = -1;
                    for (const Face& face : model.cells[k].faces) {
                        if (face.neighbor >= 0 && model.cells[face.neighbor].color == color) {
                            columnCell = face.neighbor;
                            break;
                        }
                    }
                }
                if (columnCell >= 0) {
                    for (int equation = 0; equation < 3; ++equation) {
                        const int row = 3 * k + equation;
                        matrix.at(row, 3 * columnCell + variable) = (rp[row] - rm[row]) / (2 * step);
                    }
                }
            }
        }
    }
    return matrix;
}

bool solvePoint(const Model& model, std::vector<double>& values, double vgs, double vds) {
    auto error = residual(model, values, vgs, vds);
    for (int iteration = 0; iteration < 65; ++iteration) {
        const double norm = maxResidual(error);
        if (norm < 2e-14) {
            return true;
        }
        auto correction = error;
        for (double& value : correction) {
            value = -value;
        }
        auto matrix = makeJacobian(model, values, vgs, vds);
        if (!matrix.solve(correction)) {
            return false;
        }
        bool accepted = false;
        for (double scale = 1.0; scale >= 1.0 / 4096; scale /= 2) {
            auto trial = values;
            for (std::size_t i = 0; i < trial.size(); ++i) {
                trial[i] += scale * correction[i];
            }
            auto trialError = residual(model, trial, vgs, vds);
            if (maxResidual(trialError) < norm) {
                values = std::move(trial);
                error = std::move(trialError);
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            return false;
        }
    }
    return false;
}

DcResult readResults(const Model& model, const MOSFET& device, std::vector<double> values,
                     double vgs, double vds) {
    DcResult result;
    result.vgs = vgs;
    result.vds = vds;
    const std::size_t count = model.cells.size();
    const double blank = std::numeric_limits<double>::quiet_NaN();
    for (auto* field : {&result.potentialV, &result.electronsPerM3, &result.holesPerM3,
                        &result.fieldX, &result.fieldY, &result.currentX, &result.currentY}) {
        field->resize(count, blank);
    }
    for (std::size_t k = 0; k < count; ++k) {
        const Cell& cell = model.cells[k];
        if (cell.epsilon == 0) {
            continue;
        }
        const Sample a = readSample(model, values, static_cast<int>(k));
        result.potentialV[k] = a.phi * model.ut;
        std::array<double, 4> electric{};
        std::array<double, 4> current{};
        for (int direction = 0; direction < 4; ++direction) {
            const Face& face = cell.faces[direction];
            if (face.electrostaticWeight == 0) {
                continue;
            }
            const Sample b = face.neighbor >= 0 ? readSample(model, values, face.neighbor) :
                contactSample(model, cell, face.contact, vgs, vds);
            electric[direction] = face.electrostaticWeight * (a.phi - b.phi) * model.ut /
                                  (cell.epsilon * face.length);
            if (face.contact == Contact::Gate) {
                // 从栅极指向氧化层的电位移通量，积分后得到栅极电荷。
                result.gateChargeC += face.electrostaticWeight * (b.phi - a.phi) *
                                      model.ut * device.W * 1e-6;
            }
            if (carriesCurrent(model, cell, face)) {
                const auto flux = fluxAcross(model, face, a, b);
                const double outward = elementaryCharge * (flux.electrons + flux.holes);
                current[direction] = outward / face.length;
                const double inward = -outward * device.W * 1e-6;
                if (face.contact == Contact::Drain) result.drainCurrentA += inward;
                if (face.contact == Contact::Source) result.sourceCurrentA += inward;
                if (face.contact == Contact::Body) result.bodyCurrentA += inward;
            }
        }
        result.fieldX[k] = (electric[1] - electric[0]) / 2;
        result.fieldY[k] = (electric[3] - electric[2]) / 2;
        if (cell.silicon) {
            result.electronsPerM3[k] = a.carriers.electrons;
            result.holesPerM3[k] = a.carriers.holes;
            result.currentX[k] = (current[1] - current[0]) / 2;
            result.currentY[k] = (current[3] - current[2]) / 2;
        }
    }
    result.solution = std::move(values);
    return result;
}

} // namespace

DcResult solveDc(const DeviceMesh& mesh, const MOSFET& device, double vgs, double vds,
                 const DcResult* previous) {
    if (!std::isfinite(vgs) || !std::isfinite(vds)) {
        throw std::invalid_argument("Gate and drain voltages must be finite.");
    }
    const Model model = makeModel(mesh, device);
    std::vector<double> values(model.cells.size() * 3, 0.0);
    double fromGate = 0, fromDrain = 0;
    if (previous && previous->solution.size() == values.size()) {
        values = previous->solution;
        fromGate = previous->vgs;
        fromDrain = previous->vds;
    } else {
        for (std::size_t k = 0; k < model.cells.size(); ++k) {
            if (model.cells[k].epsilon > 0) {
                values[3 * k] = std::asinh(model.cells[k].doping / (2 * model.ni));
            }
        }
        if (!solvePoint(model, values, 0, 0)) {
            throw std::runtime_error("Could not solve NMOS equilibrium.");
        }
    }

    // 逐步调高偏压。一次求解失败时把步长减半，沿用上一点的解重试。
    const double change = std::max(std::abs(vgs - fromGate), std::abs(vds - fromDrain));
    if (change > 0) {
        double step = std::min(1.0, 0.05 / change);
        double fraction = 0.0;
        while (fraction < 1.0) {
            const double next = std::min(1.0, fraction + step);
            auto trial = values;
            if (solvePoint(model, trial, fromGate + next * (vgs - fromGate),
                           fromDrain + next * (vds - fromDrain))) {
                values = std::move(trial);
                fraction = next;
            } else {
                step /= 2;
                if (step * change < 1e-5) {
                    throw std::runtime_error("DC point did not converge near Vgs=" +
                        std::to_string(fromGate + next * (vgs - fromGate)) + ", Vds=" +
                        std::to_string(fromDrain + next * (vds - fromDrain)));
                }
            }
        }
    }
    return readResults(model, device, std::move(values), vgs, vds);
}
