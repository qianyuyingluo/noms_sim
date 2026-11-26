#include "metrics.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace {

const double missing = std::numeric_limits<double>::quiet_NaN();

// 在相邻两个转移数据点间按 log10(Id) 插值，求目标电流对应的栅压。
double gateAtCurrent(const std::vector<TransferPoint>& curve, double target) {
    for (std::size_t i = 1; i < curve.size(); ++i) {
        const auto& left = curve[i - 1];
        const auto& right = curve[i];
        if (left.drainCurrentA > 0 && left.drainCurrentA <= target &&
            right.drainCurrentA >= target && right.drainCurrentA > left.drainCurrentA) {
            const double fraction = std::log10(target / left.drainCurrentA) /
                                   std::log10(right.drainCurrentA / left.drainCurrentA);
            return left.vgs + fraction * (right.vgs - left.vgs);
        }
    }
    return missing;
}

void writeMetric(std::ofstream& csv, const char* name, double value, const char* unit,
                 double vgs, double vds, double secondVds, const char* method) {
    csv << name << ',' << value << ',' << unit << ',' << vgs << ',' << vds << ','
        << secondVds << ',' << method << ',' << (std::isfinite(value) ? "ok" : "unavailable") << '\n';
}

// d(ln Id)/dVgs 就是 gm/Id。用两侧斜率估计，端点使用单侧斜率。
double logSlope(const std::vector<TransferPoint>& curve, std::size_t i) {
    const auto& point = curve[i];
    if (point.drainCurrentA <= 0) {
        return missing;
    }
    double leftStep = 0, rightStep = 0, leftSlope = 0, rightSlope = 0;
    if (i > 0 && curve[i - 1].drainCurrentA > 0) {
        leftStep = point.vgs - curve[i - 1].vgs;
        leftSlope = std::log(point.drainCurrentA / curve[i - 1].drainCurrentA) / leftStep;
    }
    if (i + 1 < curve.size() && curve[i + 1].drainCurrentA > 0) {
        rightStep = curve[i + 1].vgs - point.vgs;
        rightSlope = std::log(curve[i + 1].drainCurrentA / point.drainCurrentA) / rightStep;
    }
    if (leftStep > 0 && rightStep > 0) {
        return (rightStep * leftSlope + leftStep * rightSlope) / (leftStep + rightStep);
    }
    if (leftStep > 0) return leftSlope;
    if (rightStep > 0) return rightSlope;
    return missing;
}

} // namespace

void saveMetrics(const DeviceMesh& mesh, const MOSFET& device,
                 const DcResult& equilibrium, const DcResult& off, const DcResult& on,
                 const std::vector<TransferPoint>& transfer, const MetricSettings& settings) {
    const double lowVds = device.Vds_transfer;
    const double highVds = settings.highVds;
    const double delta = settings.differenceStep;
    if (!std::isfinite(highVds) || !std::isfinite(delta) || delta <= 0 ||
        lowVds <= 0 || highVds <= lowVds || highVds <= delta) {
        throw std::invalid_argument("Metrics need 0 < low Vds < high Vds and a positive difference step.");
    }

    const auto highTransfer = saveTransferSweep(mesh, device, equilibrium, highVds,
                                                 "output/transfer_high_vds.csv");
    const double criterion = 1e-7 * device.W / device.Lg;
    const double vthLow = gateAtCurrent(transfer, criterion);
    const double vthHigh = gateAtCurrent(highTransfer, criterion);
    const double ssLow = gateAtCurrent(transfer, criterion / 1000);
    const double ssHigh = gateAtCurrent(transfer, criterion / 10);
    const double ss = (ssHigh - ssLow) * 1000 / 2; // 两个电流截面相差两个数量级
    const double dibl = (vthLow - vthHigh) * 1000 / (highVds - lowVds);

    const double gate = on.vgs;
    const auto bias = solveDc(mesh, device, gate, highVds, &on);
    const auto gateMinus = solveDc(mesh, device, gate - delta, highVds, &bias);
    const auto gatePlus = solveDc(mesh, device, gate + delta, highVds, &bias);
    const auto drainMinus = solveDc(mesh, device, gate, highVds - delta, &bias);
    const auto drainPlus = solveDc(mesh, device, gate, highVds + delta, &bias);
    const double gm = (gatePlus.drainCurrentA - gateMinus.drainCurrentA) / (2 * delta);
    const double gds = (drainPlus.drainCurrentA - drainMinus.drainCurrentA) / (2 * delta);
    const double ro = gds > 0 ? 1 / gds : missing;
    const double cgg = (gatePlus.gateChargeC - gateMinus.gateChargeC) / (2 * delta);
    const double ratio = off.drainCurrentA > 0 ? on.drainCurrentA / off.drainCurrentA : missing;

    std::ofstream csv("output/metrics.csv");
    csv << "metric,value,unit,Vgs_V,Vds_V,second_Vds_V,method,status\n" << std::setprecision(16);
    writeMetric(csv, "Vth", vthLow, "V", vthLow, lowVds, missing, "Id=100nA*(W/L); log interpolation");
    writeMetric(csv, "Vth_high_vds", vthHigh, "V", vthHigh, highVds, missing,
                "Id=100nA*(W/L); log interpolation");
    writeMetric(csv, "Ion", on.drainCurrentA, "A", on.vgs, on.vds, missing, "Configured on bias");
    writeMetric(csv, "Ioff", off.drainCurrentA, "A", off.vgs, off.vds, missing, "Configured off bias");
    writeMetric(csv, "Ion_over_Ioff", ratio, "1", missing, on.vds, off.vds, "Ion and Ioff rows define biases");
    writeMetric(csv, "gm", gm, "S", gate, highVds, missing, "Centered gate-bias difference");
    writeMetric(csv, "ro", ro, "ohm", gate, highVds, missing, "Inverse centered drain-current slope");
    writeMetric(csv, "gm_ro", gm * ro, "V/V", gate, highVds, missing, "gm and ro at the same bias");
    writeMetric(csv, "SS", ss, "mV/dec", missing, lowVds, missing, "Two decades from 0.001 to 0.1 times threshold current");
    writeMetric(csv, "DIBL", dibl, "mV/V", missing, lowVds, highVds, "Low Vds Vth minus high Vds Vth");
    writeMetric(csv, "Cgg", cgg, "F", gate, highVds, missing, "Quasi-static centered derivative of gate charge");
    writeMetric(csv, "threshold_current", criterion, "A", missing, missing, missing, "100nA times W/L");
    writeMetric(csv, "difference_step", delta, "V", missing, missing, missing, "Half-step for centered differences");
    csv.close();
    if (!csv) {
        throw std::runtime_error("Could not write output/metrics.csv.");
    }

    std::ofstream curve("output/gm_id.csv");
    curve << "Vgs_V,Vds_V,Id_A,gm_S,gm_over_Id_per_V\n" << std::setprecision(16);
    for (std::size_t i = 0; i < transfer.size(); ++i) {
        const double slope = logSlope(transfer, i);
        const double current = transfer[i].drainCurrentA;
        curve << transfer[i].vgs << ',' << lowVds << ',' << current << ','
              << current * slope << ',' << slope << '\n';
    }
    curve.close();
    if (!curve) {
        throw std::runtime_error("Could not write output/gm_id.csv.");
    }
}
