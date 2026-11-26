#include "fields.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

void saveFields(const std::string& path, const DeviceMesh& mesh, const DcResult& point) {
    std::ofstream csv(path);
    csv << std::setprecision(16);
    csv << "row,col,x_left_um,x_right_um,y_top_um,y_bottom_um,x_um,y_um,vgs_V,vds_V,"
           "potential_V,Ex_V_m,Ey_V_m,field_V_m,electrons_m3,holes_m3,"
           "Jx_A_m2,Jy_A_m2,current_A_m2,net_doping_m3\n";
    const std::size_t columns = mesh.x.size() - 1;
    for (std::size_t row = 0; row + 1 < mesh.y.size(); ++row) {
        for (std::size_t col = 0; col < columns; ++col) {
            const std::size_t k = row * columns + col;
            const double left = mesh.x[col] * 1e6;
            const double right = mesh.x[col + 1] * 1e6;
            const double top = mesh.y[row] * 1e6;
            const double bottom = mesh.y[row + 1] * 1e6;
            const double doping = mesh.cells[k].material == Material::Silicon ?
                mesh.cells[k].netDopingPerM3 : std::numeric_limits<double>::quiet_NaN();
            csv << row << ',' << col << ',' << left << ',' << right << ',' << top << ',' << bottom << ','
                << (left + right) / 2 << ',' << (top + bottom) / 2 << ','
                << point.vgs << ',' << point.vds << ',' << point.potentialV[k] << ','
                << point.fieldX[k] << ',' << point.fieldY[k] << ','
                << std::hypot(point.fieldX[k], point.fieldY[k]) << ','
                << point.electronsPerM3[k] << ',' << point.holesPerM3[k] << ','
                << point.currentX[k] << ',' << point.currentY[k] << ','
                << std::hypot(point.currentX[k], point.currentY[k]) << ',' << doping << '\n';
        }
    }
    csv.close();
    if (!csv) {
        throw std::runtime_error("Could not write " + path);
    }
}
