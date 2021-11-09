#ifndef UTILITY_INCLUDED
#define UTILITY_INCLUDED

#include <iterator>
#include <vector>

/**
 * The namespace `ut` provides basic utility functions
 */
namespace ut {

/*
 * An iterator of [begin, end) computed on the fly.
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

}  // namespace ut

#endif