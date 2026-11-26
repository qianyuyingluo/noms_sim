#pragma once

// 由局部电势和两个准费米电势计算载流子浓度。
// 三个电势都已经除以热电压 kT/q。
struct Carriers {
    double electrons;
    double holes;
};

struct CarrierFlux {
    double electrons;
    double holes;
    double electronScale;
    double holeScale;
};

constexpr double elementaryCharge = 1.602176634e-19;
constexpr double vacuumPermittivity = 8.8541878128e-12;

double thermalVoltage(double temperature);
Carriers carrierDensity(double intrinsic, double phi, double electronLevel, double holeLevel);

// Scharfetter–Gummel 面通量除以元电荷，正方向为从 A 到 B。
CarrierFlux carrierFlux(double deltaPhi, Carriers a, Carriers b,
                        double deltaElectronLevel, double deltaHoleLevel,
                        double electronFactor, double holeFactor);
