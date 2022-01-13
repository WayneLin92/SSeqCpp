#ifndef UTILITY_INCLUDED
#define UTILITY_INCLUDED

#include <chrono>
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
 * An iterator of [begin, end) generated on the fly.
 */
class Range
{
public:
    class iterator
    {
        friend class Range;

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
        iterator operator+(int other) const
        {
            return iterator(i_ + other);
        }
        int operator-(const iterator& other) const
        {
            return i_ - other.i_;
        }

    protected:
        iterator(int start) : i_(start) {}

    private:
        int i_;
    };

    Range(int begin, int end) : begin_(begin), end_(end) {}

    iterator begin() const
    {
        return begin_;
    }
    iterator end() const
    {
        return end_;
    }

private:
    iterator begin_;
    iterator end_;
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
 * Remove elements of `cont` which are empty containers
 */
template <typename Container1d>
inline void RemoveEmptyElements(Container1d& cont)
{
    cont.erase(std::remove_if(cont.begin(), cont.end(), [](const typename Container1d::value_type& g) { return g.empty(); }), cont.end());
}

/**
 * Remove elements of `cont` which are zero
 */
template <typename Container>
inline void RemoveZeroElements(Container& cont)
{
    cont.erase(std::remove_if(cont.begin(), cont.end(), [](const typename Container::value_type& g) { return !g; }), cont.end());
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

// default = "YYYY-MM-DD HH:MM:SS"
inline std::string get_time(const std::string& fmt = "%F %T")
{
    auto bt = detail::localtime_xp(std::time(0));
    char buf[64];
    return std::string{buf, std::strftime(buf, sizeof(buf), fmt.c_str(), &bt)};
}

}  // namespace ut

#endif