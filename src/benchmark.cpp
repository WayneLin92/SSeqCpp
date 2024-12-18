#include "benchmark.h"
#include "myio.h"
#include "utility.h"
#include <fmt/format.h>
#include <fmt/color.h>

namespace bench {

std::array<int, 4> Counter::counts_ = {0, 0, 0, 0};
std::array<double, 4> AccTimer::counts_ = {0, 0, 0, 0};
std::array<int, 4> MaxGetter::max_ = {0, 0, 0, 0};

//const Counter a; /* Print when destructs. */

Timer::Timer() : bPrinted_(false)
{
    start_ = std::chrono::system_clock::now();
}

Timer::~Timer()
{
    if (!bPrinted_) {
        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start_;
        fmt::print(fmt::fg(fmt::color::green), "Time: {}s ({})\n", elapsed.count(), ut::get_time());
        if (elapsed.count() > 30.0)
            fmt::print("\a");
    }
    fflush(stdout);
}

void Timer::print(const char* msg)
{
    bPrinted_ = true;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start_;
    fmt::print(fmt::fg(fmt::color::green), "{}: {}s ({})\n", msg, elapsed.count(), ut::get_time());
    Reset();
}

std::string Timer::print2str()
{
    bPrinted_ = true;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start_;
    std::string result = fmt::format("Time: {}s ({})\n", elapsed.count(), ut::get_time());
    Reset();
    return result;
}

void AccTimer::print()
{
    fmt::print("AccTimer:\n");
    for (size_t i = 0; i < counts_.size(); ++i)
        fmt::print("{}: {}s\n", i, counts_[i]);
}

void Counter::print()
{
    fmt::print("Counter:\n");
    for (size_t i = 0; i < counts_.size(); ++i)
        fmt::print("{}: {}\n", i, counts_[i]);
}

void MaxGetter::print()
{
    fmt::print("MaxGetter:\n");
    for (size_t i = 0; i < max_.size(); ++i)
        fmt::print("{}: {}\n", i, max_[i]);
}

}  // namespace bench
