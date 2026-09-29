#pragma once

#include <string>

#include "core/model.h"

namespace pramana {

struct MpsReadOptions {
  // Upper bound given to integer columns declared between MARKER INTORG/INTEND
  // that receive no explicit bound. Industry convention (CPLEX, Gurobi, HiGHS):
  // +infinity. Some legacy readers use 1.
  double defaultIntegerUpper = kInf;
};

// Reads free or fixed MPS / QPS (optionally gzip-compressed). Supports
// OBJSENSE, RANGES (sign-dependent semantics on E rows), all BOUNDS types
// (UP/LO/FX/FR/MI/PL/BV/LI/UI), MARKER INTORG/INTEND, QUADOBJ, QMATRIX, QSECTION.
Model readMps(const std::string& path, const MpsReadOptions& opts = {});
Model readMpsFromString(const std::string& text, const MpsReadOptions& opts = {});

// Writes free MPS (QUADOBJ for Q). Round-trips through readMps exactly
// (values printed with 17 significant digits).
std::string writeMpsString(const Model& model);
void writeMps(const Model& model, const std::string& path);

}  // namespace pramana
