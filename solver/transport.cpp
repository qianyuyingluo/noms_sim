#include "transport.h"

#include <cmath>

namespace {

double bernoulli(double x) {
    if (std::abs(x) < 1e-4) {
        return 1.0 - x / 2.0 + x * x / 12.0;
    }
    if (x > 50.0) {
        return x * std::exp(-x);
    }
    return x / std::expm1(x);
}

} // namespace

double thermalVoltage(double temperature) {
    return 1.380649e-23 * temperature / elementaryCharge;
}

Carriers carrierDensity(double intrinsic, double phi, double electronLevel, double holeLevel) {
    return {intrinsic * std::exp(phi - electronLevel),
            intrinsic * std::exp(holeLevel - phi)};
}

CarrierFlux carrierFlux(double deltaPhi, Carriers a, Carriers b,
                        double deltaElectronLevel, double deltaHoleLevel,
                        double electronFactor, double holeFactor) {
    const double forward = bernoulli(deltaPhi);
    const double backward = bernoulli(-deltaPhi);

    // 从较大的一项取系数。expm1 在两项很接近时仍能保留微小电流。
    const double electronScale = electronFactor *
        (deltaElectronLevel >= 0.0 ? a.electrons * backward : b.electrons * forward);
    const double holeScale = holeFactor *
        (deltaHoleLevel >= 0.0 ? b.holes * backward : a.holes * forward);
    double electronFlux = electronScale * std::expm1(-std::abs(deltaElectronLevel));
    double holeFlux = holeScale * std::expm1(-std::abs(deltaHoleLevel));
    if (deltaElectronLevel < 0.0) {
        electronFlux = -electronFlux;
    }
    if (deltaHoleLevel < 0.0) {
        holeFlux = -holeFlux;
    }
    return {electronFlux, holeFlux, electronScale, holeScale};
}
