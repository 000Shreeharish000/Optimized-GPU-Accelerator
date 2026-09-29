#pragma once

#include <cstdarg>
#include <string>

namespace pramana {

enum class LogLevel { Quiet = 0, Error = 1, Info = 2, Detail = 3, Debug = 4 };

// Process-wide logger. Solver components log through this so the CLI can
// route everything to stderr while JSON results go to files/stdout.
class Log {
 public:
  static void setLevel(LogLevel level);
  static LogLevel level();
  static bool enabled(LogLevel level);
  static void printf(LogLevel level, const char* fmt, ...);
  static void vprintf(LogLevel level, const char* fmt, va_list args);
};

#define PLOG_INFO(...) ::pramana::Log::printf(::pramana::LogLevel::Info, __VA_ARGS__)
#define PLOG_DETAIL(...) ::pramana::Log::printf(::pramana::LogLevel::Detail, __VA_ARGS__)
#define PLOG_DEBUG(...) ::pramana::Log::printf(::pramana::LogLevel::Debug, __VA_ARGS__)
#define PLOG_ERROR(...) ::pramana::Log::printf(::pramana::LogLevel::Error, __VA_ARGS__)

std::string formatString(const char* fmt, ...);

}  // namespace pramana
