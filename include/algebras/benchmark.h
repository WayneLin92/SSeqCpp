#ifndef BENCHMARK_INCLUDED
#define BENCHMARK_INCLUDED

#include "utility.h"
#include <iostream>

namespace bench {
/**
 * The class `Timer` is used to record the time.
 */
class Timer
{
private:
    std::chrono::time_point<std::chrono::system_clock> start_;
    bool bPrinted_;
    int n_;

public:
    static std::vector<double> counts_;

public:
    Timer() : bPrinted_(false), n_(-1)
    {
        start_ = std::chrono::system_clock::now();
    }
    Timer(int n) : bPrinted_(false), n_(n)
    {
        start_ = std::chrono::system_clock::now();
    }
    ~Timer()
    {
        if (!bPrinted_) {
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed = end - start_;
            if (n_ >= 0)
                counts_[n_] += elapsed.count();
            std::cout << "\033[0;32m" << elapsed.count() << "s (" << ut::get_time() << ")\033[0m\n";
        }
    }
    void print(const char* msg = "")
    {
        bPrinted_ = true;
        auto end = std::chrono::system_clock::now();
        auto elapsed = end - start_;
        if (n_ >= 0)
            counts_[n_] += elapsed.count();
        std::cout << "\033[0;32m" << msg << elapsed.count() << "s (" << ut::get_time() << ")\033[0m\n";
        Reset();
    }
    double Elapsed() const
    {
        std::chrono::duration<double> result = std::chrono::system_clock::now() - start_;
        return result.count();
    }
    void Reset()
    {
        start_ = std::chrono::system_clock::now();
    }
    void SuppressPrint()
    {
        bPrinted_ = true;
    }

    static void print_counts_()
    {
        for (size_t i = 0; i < counts_.size(); ++i)
            std::cout << i << ": " << counts_[i] << "s\n";
    }
};

/**
* The class `Counter` records counts in a static variable.
*/
class Counter
{
public:
    static std::vector<int> counts_;

public:

    Counter(int n)
    {
        ++counts_[n];
    }

    static void print_counts_()
    {
        for (size_t i = 0; i < counts_.size(); ++i)
            std::cout << i << ": " << counts_[i] << "\n";
    }
};

}
#endif