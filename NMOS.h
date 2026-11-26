#pragma once

#include <string>
#include <vector>

// 二维硅 NMOS 器件参数
class MOSFET {
    public:
        // 1. 器件剖面尺寸。x 从硅区左边界向右，深度从硅表面向下。
        // 矩形源区的 x 范围为 [0, Ls]，矩形漏区为 [SiL - Ld, SiL]。
        double Lg;          // 栅长，um（沿源->漏）
        double W;           // 沟道宽度，um（垂直于剖面）
        double tox;         // 栅氧厚度，nm

        double SiL;         // 仿真硅区总长，um
        double SiD;         // 仿真硅区总深，um

        double Ls;          // 源区横向长度，um
        double Ld;          // 漏区横向长度，um
        double xjs;         // 源结结深，um（矩形源区的深度）
        double xjd;         // 漏结结深，um（矩形漏区的深度）

        // 栅极以硅区为中心，坐标单位为 um。
        double gateStartUm() 
        const noexcept { 
            return (SiL - Lg) / 2.0; }
        double gateEndUm() 
        const noexcept { 
            return gateStartUm() + Lg; }

        // 2. 掺杂浓度，cm^-3
        double NA;          // P 型衬底受主浓度
        double NDs;         // N+ 源区施主浓度
        double NDd;         // N+ 漏区施主浓度

        // 3. 栅极
        std::string gate_mat;   // 栅极材料，如 "n+ polysilicon"
        double phiM;            // 功函数，eV
        bool   useVFB;          // 题目直接给 VFB 时置 true
        double VFB;             // 平带电压，V

        // 4. 材料与温度
        double epsSi;       // 硅相对介电常数
        double epsOx;       // 二氧化硅相对介电常数
        double ni;          // 本征载流子浓度，cm^-3
        double mu_n;        // 电子迁移率，cm^2/(V·s)
        double mu_p;        // 空穴迁移率，cm^2/(V·s)
        double T;           // 温度，K

        // 网格划分数。硅区和栅氧层在纵向分别划分，共用交界面。
        int Nx;             // 横向网格数
        int NySi;           // 硅区纵向网格数
        int NyOx;           // 栅氧层纵向网格数
        double surfaceMeshDepthUm = 0.0; // 硅表面加密深度，0 表示关闭

        // 衬底接源极，所有工作点和扫压共用 VBS = 0 V。
        static constexpr double Vbs = 0.0;

        // 5. 六张分布图的两个工作点
        double Vgs_off, Vds_off;   // 关断
        double Vgs_on,  Vds_on;    // 导通

        // 6. 输出特性：几条固定 VGS，扫 VDS
        std::vector<double> Vgs_list;       // 固定的几个 VGS，V
        double VdsFrom, VdsTo, VdsStep;     // VDS 扫描范围与步长，V

        // 7. 转移特性：固定 VDS，扫 VGS
        double VgsFrom, VgsTo, VgsStep;     // VGS 扫描范围与步长，V
        double Vds_transfer;                // 转移特性的固定 VDS，V
};
