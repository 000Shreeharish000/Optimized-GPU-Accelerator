#include "util/log.h"

#include <atomic>
#include <cstdio>
#include <mutex>
#include <vector>

#include "util/common.h"

namespace pramana {

namespace {
std::atomic<int> gLevel{static_cast<int>(LogLevel::Info)};
std::mutex gMutex;
}  // namespace

void Log::setLevel(LogLevel level) { gLevel = static_cast<int>(level); }
LogLevel Log::level() { return static_cast<LogLevel>(gLevel.load()); }
bool Log::enabled(LogLevel level) { return static_cast<int>(level) <= gLevel.load(); }

void Log::vprintf(LogLevel level, const char* fmt, va_list args) {
  if (!enabled(level)) return;
  std::lock_guard<std::mutex> lock(gMutex);
  std::vfprintf(stderr, fmt, args);
  std::fputc('\n', stderr);
  std::fflush(stderr);
}

void Log::printf(LogLevel level, const char* fmt, ...) {
  if (!enabled(level)) return;
  va_list args;
  va_start(args, fmt);
  vprintf(level, fmt, args);
  va_end(args);
}

std::string formatString(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  va_list copy;
  va_copy(copy, args);
  int n = std::vsnprintf(nullptr, 0, fmt, copy);
  va_end(copy);
  std::vector<char> buf(static_cast<size_t>(n > 0 ? n : 0) + 1);
  std::vsnprintf(buf.data(), buf.size(), fmt, args);
  va_end(args);
  return std::string(buf.data(), static_cast<size_t>(n > 0 ? n : 0));
}

const char* statusName(Status s) {
  switch (s) {
    case Status::NotSolved: return "NOT_SOLVED";
    case Status::Optimal: return "OPTIMAL";
    case Status::Infeasible: return "INFEASIBLE";
    case Status::Unbounded: return "UNBOUNDED";
    case Status::InfeasibleOrUnbounded: return "INFEASIBLE_OR_UNBOUNDED";
    case Status::TimeLimit: return "TIME_LIMIT";
    case Status::IterationLimit: return "ITERATION_LIMIT";
    case Status::NodeLimit: return "NODE_LIMIT";
    case Status::NumericalFailure: return "NUMERICAL_FAILURE";
    case Status::Nonconvex: return "NONCONVEX";
    case Status::Error: return "ERROR";
  }
  return "UNKNOWN";
}

const char* basisStatusName(BasisStatus s) {
  switch (s) {
    case BasisStatus::Lower: return "LOWER";
    case BasisStatus::Upper: return "UPPER";
    case BasisStatus::Basic: return "BASIC";
    case BasisStatus::Zero: return "ZERO";
    case BasisStatus::Fixed: return "FIXED";
  }
  return "?";
}

}  // namespace pramana
