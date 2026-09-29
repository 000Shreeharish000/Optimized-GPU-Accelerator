// PRAMANA — common definitions shared by every module.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace pramana {

constexpr double kInf = std::numeric_limits<double>::infinity();
// MPS files use 1e30 (and some writers 1e20) to encode "infinite" bounds.
constexpr double kInfBoundThreshold = 1e20;

inline bool isInf(double v) { return std::fabs(v) >= kInfBoundThreshold; }
inline bool isFiniteBound(double v) { return std::fabs(v) < kInfBoundThreshold; }
inline double sqr(double v) { return v * v; }

class PramanaError : public std::runtime_error {
 public:
  explicit PramanaError(const std::string& msg) : std::runtime_error(msg) {}
};

#define PRAMANA_CHECK(cond, msg)                                              \
  do {                                                                        \
    if (!(cond)) throw ::pramana::PramanaError(std::string("check failed: ") + \
                                               (msg) + " [" #cond "]");       \
  } while (0)

// Final status of a solve. Every status except Optimal/Infeasible/Unbounded is
// a "no proof" status; Optimal/Infeasible/Unbounded are only reported after the
// certifier has accepted the claim (see cert/certifier.h).
enum class Status {
  NotSolved,
  Optimal,
  Infeasible,
  Unbounded,
  InfeasibleOrUnbounded,
  TimeLimit,
  IterationLimit,
  NodeLimit,
  NumericalFailure,
  Nonconvex,
  Error,
};

const char* statusName(Status s);

// Basis status of a column or row logical.
enum class BasisStatus : int8_t {
  Lower = 0,   // nonbasic at lower bound
  Upper = 1,   // nonbasic at upper bound
  Basic = 2,
  Zero = 3,    // nonbasic free variable at zero
  Fixed = 4,   // nonbasic, lower == upper
};

const char* basisStatusName(BasisStatus s);

}  // namespace pramana
