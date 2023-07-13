#ifndef UTILITY_INCLUDED
#define UTILITY_INCLUDED

#include <algorithm>
#include <chrono>
#include <future>
#include <iterator>
#include <string>
#include <vector>

#ifndef __unix__
#ifndef _MSC_VER
#include <mutex>
#endif
#endif

/**
 * The namespace `ut` provides basic utility functions
 */
namespace ut {

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
 * emulate std::map<int, T> by a vector
 * T should be a type easy to copy
 */
template <typename T, T d>
struct map_seq
{
    std::vector<T> data;
    T operator[](size_t index) const
    {
        if (index < data.size())
            return data[index];
        else
            return d;
    }
    T& operator[](size_t index)
    {
        if (index >= data.size())
            data.resize(index + 1, d);
        return data[index];
    }
    T at(size_t index) const
    {
        return (*this)[index];
    }
};

/**
 * emulate std::map<(int, int), T> by a 2d vector
 * T should be a type easy to copy
 */
template <typename T, T d>
struct map_seq2d
{
    std::vector<std::vector<T>> data;
    T operator()(size_t i, size_t j) const
    {
        if (i < data.size() && j < data[i].size())
            return data[i][j];
        else
            return d;
    }
    T& operator()(size_t i, size_t j)
    {
        if (i >= data.size())
            data.resize(i + 1);
        if (j >= data[i].size())
            data[i].resize(j + 1, d);
        return data[i][j];
    }
    T at(size_t i, size_t j) const
    {
        return (*this)(i, j);
    }
};

/**
 * Remove elements of `cont` for which the Predicate is true
 */
template <typename Container1d, typename FnPred>
void RemoveIf(Container1d& cont, FnPred pred)
{
    cont.erase(std::remove_if(cont.begin(), cont.end(), pred), cont.end());
}

/**
 * Remove elements of `cont` which are empty containers
 */
template <typename Container1d>
void RemoveEmptyElements(Container1d& cont)
{
    RemoveIf(cont, [](const typename Container1d::value_type& g) { return g.empty(); });
}

/**
 * Remove elements of `cont` whose bool() evaluates to false
 */
template <typename Container>
void RemoveZeroElements(Container& cont)
{
    RemoveIf(cont, [](const typename Container::value_type& g) { return !g; });
}

/* This function tries to avoid dest capacity way too bigger than size */
template <typename T>
void copy(std::vector<T>& tmp, std::vector<T>& dest)
{
    if (tmp.capacity() <= tmp.size() + tmp.size() / 2)
        std::swap(dest, tmp);
    else {
        dest.clear();
        dest.insert(dest.end(), tmp.cbegin(), tmp.cend());
    }
}

/* The container `map` maps a key to a collection of type T */
template <typename T>
T& get(std::vector<T>& map, size_t k)
{
    if (map.size() <= k)
        map.resize(k + 1);
    return map[k];
}

template <typename T>
const T get(std::vector<T>& map, size_t k, T default_)
{
    if (map.size() <= k)
        return default_;
    return map[k];
}

/* Obtain a vector of keys */
template <typename T>
auto get_keys(const T& map)
{
    std::vector<std::remove_cv_t<decltype(map.begin()->first)>> result;
    for (auto p = map.begin(); p != map.end(); ++p)
        result.push_back(p->first);
    return result;
}

template <typename T>
bool has(const std::vector<T>& sorted, const T& key)
{
    return std::binary_search(sorted.begin(), sorted.end(), key);
}

template <typename T, size_t N>
bool has(const std::array<T, N>& sorted, const T& key)
{
    return std::binary_search(sorted.begin(), sorted.end(), key);
}

template <typename T, typename K>
bool has(const T& map, const K& key)
{
    return map.find(key) != map.end();
}

template <typename T>
int IndexOf(const std::vector<T>& vec, const T& key)
{
    auto it = find(vec.begin(), vec.end(), key);
    return it != vec.end() ? int(it - vec.begin()) : -1;
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

/**
 * Combine two 32 bit integers into a 64 bit integer
 */
inline uint64_t Bind(uint64_t i, uint64_t j)
{
    return (i << 32) + j;
}

/**
 * Extract two 32 bit integers from a 64 bit integer
 */
inline void UnBind(uint64_t ij, uint64_t& i, uint64_t& j)
{
    i = ij >> 32;
    j = ij & ((uint64_t(1) << 32) - 1);
}

/**
 * Compute the number of 1 in the binary data
 */
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
template <typename Fn>
void for_each_seq(size_t n, Fn f)
{
    for (size_t i = 0; i < n; ++i)
        f(i);
}

/**
 * For i=0,...,n-1, execute f(i) in parallel.
 */
template <typename Fn>
void for_each_par128(size_t n, Fn f)
{
    if (n == 0)
        return;

    static constexpr size_t THREADS_MAX = 128;
    std::vector<std::future<void>> futures;
    size_t nThreads = std::min(THREADS_MAX, n);
    for (size_t t = 0; t < nThreads; ++t)
        futures.push_back(std::async(std::launch::async, [&f, t, n]() {
            for (size_t i = t; i < n; i += THREADS_MAX)
                f(i);
        }));
    for (auto& f : futures)
        f.wait();
}

/**
 * For 0<=i<j<n, execute f(i,j) in parallel.
 */
template <typename Fn>
void for_each_pair_par(size_t n, Fn f)
{
    if (n < 2)
        return;

    static constexpr size_t THREADS_MAX = 128;
    std::vector<std::future<void>> futures;
    const size_t m = (n + 1) / 2;
    const size_t m1 = (n % 2 == 0) ? m : m - 1;
    size_t nThreads = std::min(THREADS_MAX, m1);
    for (size_t k = 0; k < 2 * m - 1; ++k) {
        for (size_t t = 0; t < nThreads; ++t) {
            futures.push_back(std::async(std::launch::async, [&f, t, m1, k, m]() {
                for (size_t i = k + 1 + t; i <= k + m1; i += THREADS_MAX) {
                    if (i == k + m)
                        f(k, 2 * m - 1);
                    else {
                        size_t ii = i % (2 * m - 1);
                        size_t jj = (2 * k + (2 * m - 1) - ii) % (2 * m - 1);
                        if (ii > jj)
                            std::swap(ii, jj);
                        f(ii, jj);
                    }
                }
            }));
        }
        for (auto& f : futures)
            f.wait();
        futures.clear();
    }
}

}  // namespace ut

#endif
