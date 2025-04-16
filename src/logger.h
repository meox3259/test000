#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>

// 定义日志类，用于收集流式输出
class LogStream {
public:
  static bool verbose;
  LogStream() {
    // 获取当前时间
    auto now = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t(now);

    // 格式化时间（精确到秒）
    std::tm tm_buf;
#ifdef _WIN32
    localtime_s(&tm_buf, &now_time_t);
#else
    localtime_r(&now_time_t, &tm_buf);
#endif

    // 将时间添加到流的开头
    stream_ << "[" << std::put_time(&tm_buf, "%Y-%m-%d %H:%M:%S") << "] ";
  }

  ~LogStream() {
    // 在对象销毁时输出收集到的所有内容
    std::cout << stream_.str() << std::endl;
  }

  // 流式操作符重载，支持各种类型的输入
  template <typename T> LogStream &operator<<(const T &value) {
    stream_ << value;
    return *this;
  }

private:
  std::ostringstream stream_; // 用于收集输出内容的字符串流
};

// 定义宏，创建临时LogStream对象以支持流式输出
#define LOG LogStream()