#include <cstdio>
#include <cstring>

#include "testing.h"
#include "util/log.h"
#include "util/timer.h"

int main(int argc, char** argv) {
  pramana::Log::setLevel(pramana::LogLevel::Error);
  const char* filter = argc > 1 ? argv[1] : nullptr;
  int passed = 0, failed = 0;
  for (auto& c : ptest::registry()) {
    if (filter && !std::strstr(c.name, filter)) continue;
    pramana::Timer t;
    try {
      c.fn();
      std::printf("[PASS] %-40s %7.3fs\n", c.name, t.seconds());
      ++passed;
    } catch (ptest::Failure& f) {
      std::printf("[FAIL] %-40s %s\n", c.name, f.msg.c_str());
      ++failed;
    } catch (std::exception& e) {
      std::printf("[FAIL] %-40s exception: %s\n", c.name, e.what());
      ++failed;
    }
  }
  std::printf("\n%d passed, %d failed, %d checks\n", passed, failed, ptest::checks());
  return failed == 0 ? 0 : 1;
}
