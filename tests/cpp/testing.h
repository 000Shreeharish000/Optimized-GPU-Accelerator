// Minimal unit-test framework (no third-party dependency).
#pragma once

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace ptest {

struct Case {
  const char* name;
  std::function<void()> fn;
};
inline std::vector<Case>& registry() {
  static std::vector<Case> r;
  return r;
}
struct Registrar {
  Registrar(const char* n, std::function<void()> f) { registry().push_back({n, std::move(f)}); }
};
struct Failure {
  std::string msg;
};
inline int& checks() {
  static int c = 0;
  return c;
}

}  // namespace ptest

#define PTEST(name)                                           \
  static void name();                                         \
  static ptest::Registrar reg_##name(#name, name);            \
  static void name()

#define EXPECT(cond)                                                                            \
  do {                                                                                          \
    ++ptest::checks();                                                                          \
    if (!(cond)) throw ptest::Failure{std::string(__FILE__) + ":" + std::to_string(__LINE__) + \
                                      ": expected " #cond};                                     \
  } while (0)

#define EXPECT_NEAR(a, b, tol)                                                                          \
  do {                                                                                                  \
    ++ptest::checks();                                                                                  \
    double _a = (a), _b = (b);                                                                          \
    if (!(std::fabs(_a - _b) <= (tol) * (1.0 + std::fabs(_b))))                                         \
      throw ptest::Failure{std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " #a " = " +     \
                           std::to_string(_a) + " vs " #b " = " + std::to_string(_b)};                  \
  } while (0)

#ifndef PRAMANA_SOURCE_DIR
#define PRAMANA_SOURCE_DIR "."
#endif
inline std::string dataPath(const std::string& rel) { return std::string(PRAMANA_SOURCE_DIR) + "/" + rel; }
