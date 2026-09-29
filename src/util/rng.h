#pragma once

#include <cstdint>

namespace pramana {

// Deterministic SplitMix64 generator: identical streams on every platform, so
// perturbations and tie-breaking are reproducible run-to-run (determinism is a
// validation requirement, see docs/VALIDATION.md #20).
class Rng {
 public:
  explicit Rng(uint64_t seed = 0x9E3779B97F4A7C15ull) : state_(seed) {}
  uint64_t next() {
    uint64_t z = (state_ += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  // Uniform in [0, 1).
  double uniform() { return (next() >> 11) * (1.0 / 9007199254740992.0); }
  int below(int n) { return static_cast<int>(next() % static_cast<uint64_t>(n)); }

 private:
  uint64_t state_;
};

}  // namespace pramana
