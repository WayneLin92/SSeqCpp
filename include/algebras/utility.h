#ifndef UTILITY_INCLUDED
#define UTILITY_INCLUDED

#include <algorithm>
#include <chrono>
#include <iterator>
#include <string>
#include <vector>

#ifdef _MSC_VER
#include <execution>
#else
#include <future>
#endif

#ifndef __unix__
#ifndef _MSC_VER
#include <mutex>
#endif
#endif

/**
 * The namespace `ut` provides basic utility functions
 */
namespace ut {
     
inline constexpr int FUTURE_NUM_THREADS = 256;

/*
 * A "random access" iterator of [begin, end) generated on the fly.
 */
class Range
{
public:
    class iterator
    {
        friend class Range;

    private:
        int i_;

    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = int;
        using difference_type = int;
        using pointer = int*;
        using reference = int;

    public:
        iterator() : i_(0) {}

        int operator*() const
        {
            return i_;
        }
        const iterator& operator++()
        {
            ++i_;
            return *this;
        }
        iterator operator++(int)
        {
            iterator copy(*this);
            ++i_;
            return copy;
        }
        bool operator==(const iterator& other) const
        {
            return i_ == other.i_;
        }
        bool operator!=(const iterator& other) const
        {
            return i_ != other.i_;
        }
        bool operator<(const iterator& other) const
        {
            return i_ < other.i_;
        }
        iterator operator+(int other) const
        {
            return iterator(i_ + other);
        }
        int operator[](int other) const
        {
            return i_ + other;
        }
        int operator-(const iterator& other) const
        {
            return i_ - other.i_;
        }

    protected:
        iterator(int start) : i_(start) {}
    };

private:
    iterator begin_;
    iterator end_;

public:
    Range(int begin, int end) : begin_(begin), end_(end) {}

    iterator begin() const
    {
        return begin_;
    }
    iterator end() const
    {
        return end_;
    }
};

/**
 * Create the int array 1, 2, ..., n
 */
inline std::vector<int> range(int n)
{
    std::vector<int> result;
    for (int i = 0; i < n; ++i)
        result.push_back(i);
    return result;
}

/**
 * Remove elements of `cont` for which the Predicate is true
 */
template <typename Container1d, typename FnPred>
inline void RemoveIf(Container1d& cont, FnPred pred)
{
    cont.erase(std::remove_if(cont.begin(), cont.end(), pred), cont.end());
}

/**
 * Remove elements of `cont` which are empty containers
 */
template <typename Container1d>
inline void RemoveEmptyElements(Container1d& cont)
{
    RemoveIf(cont, [](const typename Container1d::value_type& g) { return g.empty(); });
}

/**
 * Remove elements of `cont` which are zero
 */
template <typename Container>
inline void RemoveZeroElements(Container& cont)
{
    RemoveIf(cont, [](const typename Container::value_type& g) { return !g; });
}

namespace detail {
    /* A safe way to convert time_t to std::tm */
    inline std::tm localtime_xp(std::time_t timer)
    {
        std::tm bt{};
#if defined(__unix__)
        localtime_r(&timer, &bt);
#elif defined(_MSC_VER)
        localtime_s(&bt, &timer);
#else
        static std::mutex mtx;
        std::lock_guard<std::mutex> lock(mtx);
        bt = *std::localtime(&timer);
#endif
        return bt;
    }
}  // namespace detail

/* default = "YYYY-MM-DD HH:MM:SS" */
inline std::string get_time(const std::string& fmt = "%F %T")
{
    auto bt = detail::localtime_xp(std::time(0));
    char buf[64];
    return std::string{buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt)};
}

inline uint64_t BindPair(uint32_t i, uint32_t j)
{
    return (uint64_t(i) << 32) + uint64_t(j);
}

inline void GetPair(uint64_t ij, int& i, int& j)
{
    i = int(ij >> 32);
    j = int(ij & 0x7fffffff);
}

/* Reverse the bits */
inline uint64_t Reverse(uint64_t b)
{
    b = (b & 0xFFFFFFFF00000000) >> 32 | (b & 0x00000000FFFFFFFF) << 32;
    b = (b & 0xFFFF0000FFFF0000) >> 16 | (b & 0x0000FFFF0000FFFF) << 16;
    b = (b & 0xFF00FF00FF00FF00) >> 8 | (b & 0x00FF00FF00FF00FF) << 8;
    b = (b & 0xF0F0F0F0F0F0F0F0) >> 4 | (b & 0x0F0F0F0F0F0F0F0F) << 4;
    b = (b & 0xCCCCCCCCCCCCCCCC) >> 2 | (b & 0x3333333333333333) << 2;
    b = (b & 0xAAAAAAAAAAAAAAAA) >> 1 | (b & 0x5555555555555555) << 1;
    return b;
}

inline int popcount(unsigned int i)
{
#if defined(_MSC_VER)
    return std::_Popcount(i);
#elif defined(__GNUC__)
    return __builtin_popcount(i);
#elif defined(__clang__)
    return std::__popcount(i);
#else
    i = i - ((i >> 1) & 0x55555555);                 // add pairs of bits
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);  // quads
    i = (i + (i >> 4)) & 0x0F0F0F0F;                 // groups of 8
    return (i * 0x01010101) >> 24;                   // horizontal sum of bytes
#endif
}

inline int popcount(uint64_t i)
{
#if defined(_MSC_VER)
    return std::_Popcount(i);
#elif defined(__GNUC__)
    return __builtin_popcount(i);
#elif defined(__clang__)
    return std::__popcount(i);
#else
    i = i - ((i >> 1) & 0x55555555);                 // add pairs of bits
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);  // quads
    i = (i + (i >> 4)) & 0x0F0F0F0F;                 // groups of 8
    return (i * 0x01010101) >> 24;                   // horizontal sum of bytes
#endif
}

/**
 * For i=0,...,n-1, execute f(i) in sequence.
 */
template <typename FnOp>
void for_each_seq(int n, FnOp f)
{
    for (int i = 0; i < n; ++i)
        f(i);
}

/**
 * For i=0,...,n-1, execute f(i) in parallel.
 */
template <typename FnOp>
void for_each_par(int n, FnOp f)
{
#ifdef _MSC_VER
    ut::Range r(0, n);
    std::for_each(std::execution::par_unseq, r.begin(), r.end(), f);
#else
    std::vector<std::future<void>> futures;
    for (int i = 0; i < FUTURE_NUM_THREADS; ++i)
        futures.push_back(std::async(std::launch::async, [i, n, f]() {
            for (int j = i; j < n; j += FUTURE_NUM_THREADS)
                f(i);
        }));
    for (int i = 0; i < FUTURE_NUM_THREADS; ++i)
        futures[i].wait();
#endif
}

}  // namespace ut

#endif