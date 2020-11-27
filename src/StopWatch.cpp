#include "StopWatch.hpp"


StopWatch::StopWatch()
{
    start();
}

void StopWatch::start()
{
    mStart = std::chrono::high_resolution_clock::now();
}

double StopWatch::restart()
{
    double milliseconds = 0;

    auto now = std::chrono::high_resolution_clock::now();
    milliseconds
        = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(now - mStart).count()) / 1000;
    mStart = now;

    return milliseconds;
}

double StopWatch::stop()
{
    auto start = mStart;

    auto milliseconds = restart();

    mStart = start;

    return milliseconds;
}
