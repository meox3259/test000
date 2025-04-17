#pragma once

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

class LogStream {
public:
  static bool verbose;
  LogStream() {
    auto now = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t(now);

    std::tm tm_buf;
#ifdef _WIN32
    localtime_s(&tm_buf, &now_time_t);
#else
    localtime_r(&now_time_t, &tm_buf);
#endif

    stream_ << "[" << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S") << "] ";
  }

  ~LogStream() {
    if (verbose) {
      std::cout << stream_.str() << std::endl;
    }
  }

  template <typename T> LogStream &operator<<(const T &value) {
    stream_ << value;
    return *this;
  }

private:
  std::ostringstream stream_;
};

#define LOG LogStream()