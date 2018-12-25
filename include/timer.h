#pragma once

#include <string>
#include <chrono>

struct timer {
    std::chrono::system_clock::time_point nt;

    timer() : nt(std::chrono::high_resolution_clock::now()) {}

    double get_seconds() {
        auto t = std::chrono::high_resolution_clock::now();
        auto millisecond = std::chrono::duration<double, std::milli>(t - nt).count();
        auto second = millisecond / 1000;
        nt = t;
        return second;
    }
};
