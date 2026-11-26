#include "sweep.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {

std::vector<double> voltageSteps(double start, double end, double step) {
    if (!std::isfinite(start) || !std::isfinite(end) ||
        !std::isfinite(step) || step <= 0 || start > end) {
        throw std::invalid_argument("Voltage scan needs a positive step and start <= end.");
    }
    std::vector<double> values;
    for (int i = 0; i < 100000; ++i) {
        const double voltage = start + i * step;
        if (voltage >= end - 1e-12) {
            values.push_back(end);
            return values;
        }
        values.push_back(voltage);
    }
    throw std::invalid_argument("Too many voltage points.");
}

void writePoint(std::ofstream& csv, const DcResult& point) {
    csv << point.vgs << ',' << point.vds << ',' << point.drainCurrentA << ','
        << point.sourceCurrentA << ',' << point.bodyCurrentA << '\n';
}

} // namespace

void saveOutputSweep(const DeviceMesh& mesh, const MOSFET& device, const DcResult& equilibrium) {
    std::ofstream csv("output/output_iv.csv");
    csv << "Vgs_V,Vds_V,Id_A,Is_A,Ib_A\n" << std::setprecision(16);
    const auto voltages = voltageSteps(device.VdsFrom, device.VdsTo, device.VdsStep);
    for (double gate : device.Vgs_list) {
        auto point = equilibrium;
        for (double drain : voltages) {
            point = solveDc(mesh, device, gate, drain, &point);
            writePoint(csv, point);
        }
    }
    csv.close();
    if (!csv) {
        throw std::runtime_error("Could not write output/output_iv.csv.");
    }
}

std::vector<TransferPoint> saveTransferSweep(const DeviceMesh& mesh, const MOSFET& device,
                                              const DcResult& equilibrium, double vds,
                                              const char* path) {
    std::ofstream csv(path);
    if (!csv) {
        throw std::runtime_error(std::string("Could not open ") + path);
    }
    csv << "Vgs_V,Vds_V,Id_A,Is_A,Ib_A\n" << std::setprecision(16);
    auto point = equilibrium;
    std::vector<TransferPoint> samples;
    for (double gate : voltageSteps(device.VgsFrom, device.VgsTo, device.VgsStep)) {
        point = solveDc(mesh, device, gate, vds, &point);
        writePoint(csv, point);
        samples.push_back({point.vgs, point.drainCurrentA});
    }
    csv.close();
    if (!csv) {
        throw std::runtime_error(std::string("Could not write ") + path);
    }
    return samples;
}
