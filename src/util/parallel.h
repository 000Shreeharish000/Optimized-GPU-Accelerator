// Tiny deterministic parallel-for over a persistent pool of std::threads.
// Work is split into contiguous static chunks, so floating-point results of a
// chunked reduction are identical run-to-run for a fixed thread count.
#pragma once

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace pramana {

class ThreadPool {
 public:
  static ThreadPool& instance();
  static void setThreads(int n);  // call before first use; 0 = hardware
  int threads() const { return static_cast<int>(workers_.size()) + 1; }

  // Runs f(chunkIndex, begin, end) over [0,n) split into threads() chunks.
  void parallelFor(int n, const std::function<void(int, int, int)>& f, int minPerChunk = 2048);

  ~ThreadPool();

 private:
  explicit ThreadPool(int n);
  void workerLoop(int id);

  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable cv_, done_;
  const std::function<void(int, int, int)>* job_ = nullptr;
  int jobN_ = 0, jobChunks_ = 0;
  uint64_t generation_ = 0;
  int pending_ = 0;
  bool stop_ = false;
};

}  // namespace pramana
