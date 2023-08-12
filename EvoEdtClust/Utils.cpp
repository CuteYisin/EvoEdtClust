#include "Utils.h"

int Config::n_threads = 1;

std::chrono::high_resolution_clock::time_point Timer::set_start() {
    return std::chrono::high_resolution_clock::now();
}

void Timer::time_profile(const std::string & info, const std::chrono::high_resolution_clock::time_point & start) {
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
    if(time_elapsed >= 10) {
        std::cerr << info << " " << time_elapsed << "ms" << std::endl;
    }
}