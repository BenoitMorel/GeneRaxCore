
#include <chrono>
#include <thread>
#include <cstdio>

using Time = std::chrono::time_point<std::chrono::system_clock>;

class Timer {
public:
  Timer() {
    reset();
  }

  long getElapsedMs() const {
    auto end = getTime();
    return getElapsedMs(_start, end);
  }

  void reset() {
    _start = std::chrono::system_clock::now();
  }
  
  static Time getTime() {
    return std::chrono::system_clock::now();
  }
  
  static long getElapsedMs(Time begin, Time end) {
    return std::chrono::duration_cast<std::chrono::milliseconds>
      (end-begin).count();
  }
private:
  Time _start;
};
