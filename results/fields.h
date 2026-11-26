#pragma once

#include "../solver/dc_solver.h"

#include <string>

void saveFields(const std::string& path, const DeviceMesh& mesh, const DcResult& point);
