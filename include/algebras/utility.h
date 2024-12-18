#ifndef UTILITY_INCLUDED
#define UTILITY_INCLUDED

#include <algorithm>
#include <chrono>
#include <future>
#include <iterator>
#include <string>
#include <vector>
#include <array>

#ifndef __unix__
#ifndef _MSC_VER
#include <mutex>
#endif
#endif

/**
 * The namespace `ut` provides basic utility functions
 */
namespace ut {

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

template <typename T, size_t N>
class vector
{
protected:
    std::array<T, N> data_;
    size_t size_ = 0;

public:
    constexpr auto size() const
    {
        return size_;
    }
    auto begin() const
    {
        return data_.begin();
    }
    auto end() const
    {
        return data_.begin() + size_;
    }
    auto rbegin() const
    {
        return data_.rbegin() + N - size_;
    }
    auto rend() const
    {
        return data_.rend();
    }
    auto& operator[](size_t i) const
    {
        return data_[i];
    }
    auto& front() const
    {
        return data_[0];
    }
    auto& front()
    {
        return data_[0];
    }
    auto& back() const
    {
        return data_[size_ - 1];
    }
    auto& back()
    {
        return data_[size_ - 1];
    }
    auto& operator[](size_t i)
    {
#ifdef MYDEBUG
        if (i >= size_) {
            printf("ut::vector out of scope\n");
            throw "DEBUG";
        }
#endif
        return data_[i];
    }

    void push_back(T p)
    {
        if (size_ >= N) {
            printf("ut::vector overflow\n");
            // throw "DEBUG";
            std::exit(1974285);
        }
        data_[size_++] = std::move(p);
    }
    void pop_back()
    {
        --size_;
    }
    void clear()
    {
        size_ = 0;
    }
};

// TODO: Deprecated. Use ut::get instead
/*
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

/** TODO: Deprecated. Use ut::get instead
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

/* An iterator adapter for 2d data */
using IJ = std::pair<size_t, size_t>;
template <typename Cont2d>
class Iter2d
{
private:
    Cont2d* pData;

public:
    class iterator
    {
        friend class Iter2d;

    private:
        const Cont2d* pData;
        IJ p;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = IJ;
        using pointer = IJ;
        using reference = IJ;
        using difference_type = std::ptrdiff_t;

    public:
        /* needs to be default contructible to be compatible with std::ranges */
        iterator() : pData(nullptr), p({0, 0}) {}

        reference operator*() const
        {
            return p;
        }
        iterator& operator++()
        {
            if (++p.second < (*pData)[p.first].size())
                return *this;
            while (++p.first < pData->size() && (*pData)[p.first].empty())
                ;
            p.second = 0;
            if (p.first >= pData->size())
                p.first = size_t(-1);
            return *this;
        }
        iterator operator++(int)
        {
            iterator copy(*this);
            ++(*this);
            return copy;
        }
        auto operator==(const iterator& rhs) const
        {
            return p == rhs.p;
        }
        auto operator!=(const iterator& rhs) const
        {
            return p != rhs.p;
        }

    protected:
        iterator(const Cont2d* pData) : pData(pData), p({0, 0})
        {
            while (p.first < pData->size() && (*pData)[p.first].empty())
                ++p.first;
            if (p.first >= pData->size())
                p.first = size_t(-1);
        }

        iterator(size_t pData_size) : pData(nullptr), p({pData_size, 0}) {}
    };

    Iter2d(Cont2d* pData) : pData(pData) {}

    iterator begin() const
    {
        return iterator(pData);
    }

    /* As we access the elements by indices,
     * we allow the container to be extended during iteration.
     */
    iterator end() const
    {
        return iterator(size_t(-1));
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
    if (tmp.capacity() * 2 <= tmp.size() + tmp.size())
        std::swap(dest, tmp);
    else {
        dest.clear();
        dest.insert(dest.end(), tmp.cbegin(), tmp.cend());
    }
}

/* The container `vec` maps a key to a collection of type T */
template <typename T>
void extend(std::vector<T>& vec, size_t sz)
{
    if (vec.size() < sz)
        vec.resize(sz);
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

template <typename Maps, typename K> //// TODO: Remove
const auto& GetRecentValue(Maps& maps, const K& key)
{
    for (auto p = maps.rbegin(); p != maps.rend(); ++p)
        if (auto result = p->find(key); result != p->end())
            return result->second;
    throw std::out_of_range("Recent Value not found");
}

template <typename T>
int IndexOf(const std::vector<T>& vec, const T& key)
{
    auto it = std::find(vec.begin(), vec.end(), key);
    return it != vec.end() ? int(it - vec.begin()) : -1;
}

template <typename T, typename FnEq>
int IndexOf(const std::vector<T>& vec, FnEq fnEq)
{
    auto it = std::find_if(vec.begin(), vec.end(), fnEq);
    return it != vec.end() ? int(it - vec.begin()) : -1;
}

template <typename T>
int IndexOfInSorted(const std::vector<T>& sorted, const T& key)
{
    auto it = std::lower_bound(sorted.begin(), sorted.end(), key);
    return it != sorted.end() && *it == key ? int(it - sorted.begin()) : -1;
}

/* The container `map` maps a key to a collection of type T */
template <typename T>
T& get(std::vector<T>& map, size_t k)
{
    if (map.size() <= k)
        map.resize(k + 1);
    return map[k];
}

/* T should be a small type */
template <typename T>
const T get(const std::vector<T>& map, size_t k, T default_)
{
    if (map.size() <= k)
        return default_;
    return map[k];
}

/* Obtain a vector of keys */
template <typename Map>
auto get_keys(const Map& map)
{
    std::vector<std::remove_cv_t<decltype(map.begin()->first)>> result;
    for (auto p = map.begin(); p != map.end(); ++p)
        result.push_back(p->first);
    return result;
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

/* If n = 2^k1 + ... + 2^kn,
 * return the array k1, ..., kn. */
inline ut::vector<int, 32> two_exp(unsigned n)
{
    ut::vector<int, 32> result;
    int k = 0;
    while (n > 0) {
        if (n & 1)
            result.push_back(k);
        n >>= 1;
        ++k;
    }
    return result;
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

inline int popcount(uint8_t i)
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
template <size_t NumThreads, typename Fn>
void for_each_par(size_t n, Fn f)
{
    if (n == 0)
        return;

    std::vector<std::future<void>> futures;
    size_t nThreads = std::min(NumThreads, n);
    for (size_t t = 0; t < nThreads; ++t)
        futures.push_back(std::async(std::launch::async, [&f, t, n]() {
            /* If the threads throw exceptions the program exits with return value 1111 */
            try {
                for (size_t i = t; i < n; i += NumThreads)
                    f(i);
            }
            catch (const std::exception& e) {
                printf("Exception in for_each_par: %s\n", e.what());
                std::exit(1111);
            }
        }));
    for (auto& fu : futures)
        fu.wait();
}

/**
 * For i=0,...,n-1, execute f(i) in parallel.
 */
template <typename Fn>
void for_each_par32(size_t n, Fn f)
{
    for_each_par<32>(n, f);
}

/**
 * For i=0,...,n-1, execute f(i) in parallel.
 */
template <typename Fn>
void for_each_par128(size_t n, Fn f)
{
    for_each_par<128>(n, f);
}

}  // namespace ut

#endif
