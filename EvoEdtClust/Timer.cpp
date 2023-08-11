#include "Timer.h"

std::chrono::high_resolution_clock::time_point Timer::set_start() {
    return std::chrono::high_resolution_clock::now();
}

void Timer::time_profile(const std::string & info, const std::chrono::high_resolution_clock::time_point & start) {
    std::cerr << info << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count() << "ms" << std::endl;
}