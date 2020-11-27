#pragma once

#include <chrono>

/// Simple Stopwatch class
class StopWatch
{
private:
    /// First time point.
    std::chrono::high_resolution_clock::time_point mStart;

public:
    /// Start is automatically called.
    StopWatch();

    /// Start stop watch.
    void start();

    /// Stop and return time in milliseconds.
    double stop();

    /// Restarts timer and returns time until now.
    double restart();
};
