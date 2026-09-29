#pragma once

#include <atomic>
#include <chrono>

namespace pramana {

class Timer {
 public:
  Timer() { reset(); }
  void reset() { start_ = Clock::now(); }
  double seconds() const {
    return std::chrono::duration<double>(Clock::now() - start_).count();
  }

 private:
  using Clock = std::chrono::steady_clock;
  Clock::time_point start_;
};

// Deadline shared by nested algorithms (B&B -> node LP -> simplex).
class Deadline {
 public:
  Deadline() = default;
  explicit Deadline(double limitSeconds) : limit_(limitSeconds) {}
  // A shared cancellation flag (used by "race" mode to stop the losing engine).
  void setCancelFlag(const std::atomic<bool>* flag) { cancel_ = flag; }
  bool expired() const {
    if (cancel_ && cancel_->load(std::memory_order_relaxed)) return true;
    return limit_ > 0 && timer_.seconds() >= limit_;
  }
  bool cancelled() const { return cancel_ && cancel_->load(std::memory_order_relaxed); }
  double remaining() const { return limit_ > 0 ? limit_ - timer_.seconds() : 1e300; }
  double elapsed() const { return timer_.seconds(); }
  double limit() const { return limit_; }

 private:
  Timer timer_;
  double limit_ = 0;  // <= 0 means unlimited
  const std::atomic<bool>* cancel_ = nullptr;
};

}  // namespace pramana
