#ifndef UTILITY_INCLUDED
#define UTILITY_INCLUDED

#include <algorithm>
#include <chrono>
#include <iterator>
#include <string>
#include <vector>
#include <future>

#ifndef __unix__
#ifndef _MSC_VER
#include <mutex>
#endif
#endif

/**
 * The namespace `ut` provides basic utility functions
 */
namespace ut {

inline constexpr size_t FUTURE_NUM_THREADS = 256;

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
        size_t i_;

    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = size_t;
        using difference_type = size_t;
        using pointer = size_t*;
        using reference = size_t;

    public:
        iterator() : i_(0) {}

        size_t operator*() const
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
        iterator operator+(size_t other) const
        {
            return iterator(i_ + other);
        }
        size_t operator[](size_t other) const
        {
            return i_ + other;
        }
        size_t operator-(const iterator& other) const
        {
            return i_ - other.i_;
        }

    protected:
        iterator(size_t start) : i_(start) {}
    };

private:
    iterator begin_;
    iterator end_;

public:
    Range(size_t begin, size_t end) : begin_(begin), end_(end) {}

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
inline std::vector<int> int_range(int n)
{
    std::vector<int> result;
    result.reserve(n);
    for (int i = 0; i < n; ++i)
        result.push_back(i);
    return result;
}

/**
 * Create the int array 1, 2, ..., n
 */
inline std::vector<size_t> size_t_range(size_t n)
{
    std::vector<size_t> result;
    result.reserve(n);
    for (size_t i = 0; i < n; ++i)
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

inline uint64_t Bind(uint64_t i, uint64_t j)
{
    return (i << 32) + j;
}

inline void UnBind(uint64_t ij, uint64_t& i, uint64_t& j)
{
    i = ij >> 32;
    j = ij & ((uint64_t(1) << 32) - 1);
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
    return std::__popcount(i);
#endif
}

inline int popcount(uint64_t i)
{
#if defined(_MSC_VER)
    return std::_Popcount(i);
#elif defined(__GNUC__)
    return __builtin_popcountll(i);
#elif defined(__clang__)
    return std::__popcountll(i);
#else
    return std::__popcountll(i);
#endif
}

/**
 * For i=0,...,n-1, execute f(i) in sequence.
 */
template <typename FnOp>
void for_each_seq(size_t n, FnOp f)
{
    for (size_t i = 0; i < n; ++i)
        f(i);
}

/**
 * For i=0,...,n-1, execute f(i) in parallel.
 */
template <typename FnOp>
void for_each_par(size_t n, FnOp f)
{
    std::vector<std::future<void>> futures;
    size_t nThreads = std::min(FUTURE_NUM_THREADS, n);
    for (int i = 0; i < nThreads; ++i)
        futures.push_back(std::async(std::launch::async, [i, n, f]() {
            for (int j = i; j < n; j += FUTURE_NUM_THREADS)
                f(i);
        }));
    for (int i = 0; i < nThreads; ++i)
        futures[i].wait();
}

inline void hash_combine(uint64_t& seed, uint64_t value)
{
    seed ^= value + 0x9e3779b97f4a7c15 + (seed << 12) + (seed >> 4);
}

}  // namespace ut

#endif