#include "util/parallel.h"

#include <algorithm>

namespace pramana {

namespace {
int gRequestedThreads = 0;
}

void ThreadPool::setThreads(int n) { gRequestedThreads = n; }

ThreadPool& ThreadPool::instance() {
  static ThreadPool pool(gRequestedThreads > 0
                             ? gRequestedThreads
                             : std::max(1, static_cast<int>(std::thread::hardware_concurrency())));
  return pool;
}

ThreadPool::ThreadPool(int n) {
  for (int i = 1; i < n; ++i) workers_.emplace_back([this, i] { workerLoop(i); });
}

ThreadPool::~ThreadPool() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    stop_ = true;
  }
  cv_.notify_all();
  for (auto& t : workers_) t.join();
}

static void runChunk(const std::function<void(int, int, int)>& f, int n, int chunks, int c) {
  int begin = static_cast<int>(static_cast<long long>(n) * c / chunks);
  int end = static_cast<int>(static_cast<long long>(n) * (c + 1) / chunks);
  if (begin < end) f(c, begin, end);
}

void ThreadPool::workerLoop(int id) {
  uint64_t seen = 0;
  for (;;) {
    const std::function<void(int, int, int)>* job;
    int n, chunks;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_.wait(lock, [&] { return stop_ || generation_ != seen; });
      if (stop_) return;
      seen = generation_;
      job = job_;
      n = jobN_;
      chunks = jobChunks_;
    }
    if (id < chunks) runChunk(*job, n, chunks, id);
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (--pending_ == 0) done_.notify_one();
    }
  }
}

void ThreadPool::parallelFor(int n, const std::function<void(int, int, int)>& f, int minPerChunk) {
  int chunks = std::min(threads(), std::max(1, n / std::max(1, minPerChunk)));
  if (chunks <= 1 || workers_.empty()) {
    if (n > 0) f(0, 0, n);
    return;
  }
  {
    std::lock_guard<std::mutex> lock(mutex_);
    job_ = &f;
    jobN_ = n;
    jobChunks_ = chunks;
    pending_ = static_cast<int>(workers_.size());
    ++generation_;
  }
  cv_.notify_all();
  runChunk(f, n, chunks, 0);
  std::unique_lock<std::mutex> lock(mutex_);
  done_.wait(lock, [&] { return pending_ == 0; });
}

}  // namespace pramana
