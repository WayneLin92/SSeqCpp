#ifndef BENCHMARK_INCLUDED
#define BENCHMARK_INCLUDED
#include <vector>
#include <chrono>
#include <string>

namespace bench {
/**
 * The class `Timer` is used to record the time.
 */
class Timer
{
private:
    std::chrono::time_point<std::chrono::system_clock> start_;
    bool bPrinted_;

public:
    Timer();
    ~Timer();
    void print(const char* msg = "");
    std::string print2str();
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
};

/**
 * The class `Timer` is used to accumulate the time.
 */
class AccTimer
{
private:
    std::chrono::time_point<std::chrono::system_clock> start_;
    bool ended_;
    int n_;

public:
    static std::vector<double> counts_;

public:
    AccTimer(int n) : ended_(false), n_(n)
    {
        start_ = std::chrono::system_clock::now();
    }
    ~AccTimer()
    {
        if (!ended_) {
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed = end - start_;
            counts_[n_] += elapsed.count();
        }
    }
    void end()
    {
        ended_ = true;
    }
    double Elapsed() const
    {
        std::chrono::duration<double> result = std::chrono::system_clock::now() - start_;
        return result.count();
    }
    void Reset(int n = -1)
    {
        start_ = std::chrono::system_clock::now();
        ended_ = false;
        n_ = n == -1 ? n : n_;
    }
    static void print();
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

    static void print();
};

/**
 * The class `Counter` records counts in a static variable.
 */
class MaxGetter
{
public:
    static std::vector<int> max_;

public:
    MaxGetter(int n, int value)
    {
        if (value > max_[n])
            max_[n] = value;
    }

    static void print();
};

}  // namespace bench
#endif