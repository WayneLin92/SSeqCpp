/** \file steenrod.h
 * A Component for monomials and polynomials over $F_2$.
 * All are incapsulated in the `namespace alg`.
 */

#ifndef STEENROD_H
#define STEENROD_H

#include "myexception.h"
#include "utility.h"
#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <vector>

namespace steenrod {
using array = std::vector<int>;
using array2d = std::vector<array>;
using array3d = std::vector<array2d>;

struct Milnor;
class MMay;
struct May;

inline constexpr size_t MILNOR_BUFFER_SIZE = 9;                     /* Support up to \xi_8 */
inline constexpr int DEG_MAX = (1 << (MILNOR_BUFFER_SIZE + 1)) - 1; /* Maximum degree supported in A */

struct MMilnor
{
    using TypeData = std::array<int, MILNOR_BUFFER_SIZE>;
    TypeData data = {};

    MMilnor() = default;
    MMay ToMMay() const;

    static MMilnor P(int i, int j) /* P_{ij} */
    {
        MMilnor result;
        result.data[size_t(j - i - 1)] = (1 << i);
        return result;
    }

    bool operator<(const MMilnor& rhs) const
    {
        int w1 = weight(), w2 = rhs.weight();
        if (w1 < w2)
            return true;
        if (w2 < w1)
            return false;
        if (data < rhs.data)
            return true;
        return false;
    };

    int weight() const
    {
        int result = 0;
        for (size_t i = 0; i < MILNOR_BUFFER_SIZE; ++i)
            result += (2 * (int)i + 1) * std::_Popcount(unsigned(data[i]));
        return result;
    }

    int deg() const
    {
        int result = 0;
        for (size_t i = 0; i < MILNOR_BUFFER_SIZE; ++i)
            result += ((1 << (i + 1)) - 1) * data[i];
        return result;
    }

    Milnor operator*(const MMilnor& rhs) const;
};
using MonMilnor1d = std::vector<MMilnor>;

struct Milnor
{
    MonMilnor1d data;

    Milnor() {}
    Milnor(const MMilnor& r) : data({r}) {}
    static Milnor P(int i, int j) /* P_{ij} */
    {
        return Milnor(MMilnor::P(i, j));
    }
    May ToMay() const;

    static Milnor Unit()
    {
        return Milnor(MMilnor{{0}});
    }

    Milnor operator+(const Milnor& rhs) const
    {
        Milnor result;
        std::set_symmetric_difference(data.begin(), data.end(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data));
        return result;
    }
    Milnor& operator+=(const Milnor& rhs)
    {
        Milnor tmp;
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data));
        return *this;
    }
    Milnor operator*(const Milnor& rhs) const
    {
        Milnor result;
        for (auto& r : data)
            for (auto& r1 : rhs.data)
                result += r * r1;
        return result;
    }
};
std::ostream& operator<<(std::ostream& sout, const Milnor& x);

inline constexpr size_t MMAY_INDEX_NUM = (MILNOR_BUFFER_SIZE + 1) * (MILNOR_BUFFER_SIZE + 2) / 2 - 1;
static_assert(MMAY_INDEX_NUM < sizeof(uint64_t) * 8);

namespace detail {
    inline constexpr std::array<int, MMAY_INDEX_NUM> MonSteenrodMayI()
    {
        std::array<int, MMAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)MILNOR_BUFFER_SIZE + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMAY_INDEX_NUM)
                    result[n++] = i;
        return result;
    }
    inline constexpr std::array<int, MMAY_INDEX_NUM> MonSteenrodMayJ()
    {
        std::array<int, MMAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)MILNOR_BUFFER_SIZE + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMAY_INDEX_NUM)
                    result[n++] = j;
        return result;
    }
}  // namespace detail

/* May basis for A
** 
** Each element is represented by a 64-bit unsigned integer
*/
inline constexpr std::array<int, MMAY_INDEX_NUM> MMAY_GEN_I = detail::MonSteenrodMayI();
inline constexpr std::array<int, MMAY_INDEX_NUM> MMAY_GEN_J = detail::MonSteenrodMayJ();
inline constexpr uint64_t MMAY_ONE = uint64_t(1) << (MMAY_INDEX_NUM - 1);
inline constexpr uint64_t MMAY_LEFT_BIT = uint64_t(1) << 63;
inline constexpr uint64_t MMAY_MASK_W = 0xffc0000000000000;
inline constexpr uint64_t MMAY_MASK_M = 0x003fffffffffffff;
inline constexpr uint64_t MMAY_NULL = 0xffffffffffffffff;

class MMay
{
private:
    uint64_t data_;

public:
    /**
    * The first 10 bits are used to store the weight.
    * The rest are used to store the exterior monomial.
    */
    MMay() : data_(0) {}
    static uint64_t rawP(int i, int j)
    {
        return MMAY_ONE >> (j * (j + 1) / 2 - i - 1);
    }
    explicit MMay(uint64_t data) : data_(data) {}
    static MMay FromIndex(int i)
    {
        return MMay((MMAY_ONE >> i) + (uint64_t(2 * (MMAY_GEN_J[i] - MMAY_GEN_I[i]) - 1) << MMAY_INDEX_NUM));
    }
    static MMay P(int i, int j)
    {
        return MMay((MMAY_ONE >> (j * (j + 1) / 2 - i - 1)) + (uint64_t(2 * (j - i) - 1) << MMAY_INDEX_NUM));
    }

    MMilnor ToMMilnor() const
    {
        MMilnor result;
        for (int i : *this)
            result.data[size_t(MMAY_GEN_J[i] - MMAY_GEN_I[i] - 1)] += 1 << MMAY_GEN_I[i];
        return result;
    }

    bool operator<(MMay rhs) const
    {
        return data_ < rhs.data_;
    };

    bool operator==(MMay rhs) const
    {
        return data_ == rhs.data_;
    };

    explicit operator bool() const
    {
        return data_;
    }

    MMay divLF(MMay rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        uint64_t m1 = data_ & MMAY_MASK_M;
        uint64_t m2 = rhs.data_ & MMAY_MASK_M;
        if (m1 < m2 || m2 & (m1 - m2))
            throw MyException(0x7ed0a8U, "Not divisible: m1 / m2");
#endif
        return MMay(((data_ ^ rhs.data_) & MMAY_MASK_M) + ((data_ & MMAY_MASK_W) - (rhs.data_ & MMAY_MASK_W)));
    }

    bool divisibleLF(MMay rhs)
    {
        uint64_t m1 = data_ & MMAY_MASK_M;
        uint64_t m2 = rhs.data_ & MMAY_MASK_M;
        return m2 > m1 && !(m1 & (m2 - m1));
    }

    MMay gcdLF(MMay rhs)
    {
        uint64_t m = data_ & rhs.data_ & MMAY_MASK_M;
        return MMay(m + WeightBits(m));
    }

    MMay lcmLF(MMay rhs)
    {
        uint64_t m = (data_ & MMAY_MASK_M) | (rhs.data_ & MMAY_MASK_M);
        return MMay(m + WeightBits(m));
    }

    int weight() const
    {
        return (int)(data_ >> MMAY_INDEX_NUM);
    }

    static uint64_t WeightBits(uint64_t data)
    {
        uint64_t result = 0;
        int i = 0;
        for (uint64_t m = data << (64 - MMAY_INDEX_NUM); m; m <<= 1, ++i)
            if (m & MMAY_LEFT_BIT)
                result += uint64_t(2 * (MMAY_GEN_J[i] - MMAY_GEN_I[i]) - 1);
        return result << MMAY_INDEX_NUM;
    }

    int deg() const
    {
        int result = 0;
        for (int i : *this)
            result += (1 << MMAY_GEN_J[i]) - (1 << MMAY_GEN_I[i]);
        return result;
    }

public:
    class iterator
    {
        friend class MMay;

    private:
        uint64_t m_;
        int i_;

    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = int;
        using pointer = int*;
        using reference = int;

    public:
        iterator() : m_(0), i_(0) {}

        int operator*() const
        {
            return i_;
        }
        const iterator& operator++()
        {
            for (m_ <<= 1, ++i_; m_; m_ <<= 1, ++i_)
                if (m_ & MMAY_LEFT_BIT)
                    break;
            return *this;
        }
        iterator operator++(int)
        {
            iterator copy(*this);
            for (; m_; m_ <<= 1, ++i_)
                if (m_ & MMAY_LEFT_BIT)
                    break;
            return copy;
        }
        bool operator==(const iterator& rhs) const
        {
            return m_ == rhs.m_;
        }
        bool operator!=(const iterator& rhs) const
        {
            return m_ != rhs.m_;
        }

    protected:
        iterator(uint64_t m) : m_(m), i_(0)
        {
            for (; m_; m_ <<= 1, ++i_)
                if (m_ & MMAY_LEFT_BIT)
                    break;
        }
    };

    iterator begin() const
    {
        return iterator(data_ << (64 - MMAY_INDEX_NUM));
    }
    iterator end() const
    {
        return iterator(0);
    }
};
using MMay1d = std::vector<MMay>;
using MMay2d = std::vector<MMay1d>;


/* Elements of A as linear combinations of May basis 
*/
struct May
{
    MMay1d data;
    May() {}
    explicit May(MMay m) : data({m}) {}
    Milnor ToMilnor() const;

    /* Convert from Milnor basis to May basis */
    explicit May(const Milnor& x)
    {
        for (auto& m : x.data)
            data.push_back(m.ToMMay());
        std::sort(data.begin(), data.end());
    }

    static May P(int i, int j)
    {
        return May(MMay::P(i, j));
    }

    May operator+(const May& rhs) const
    {
        May result;
        std::set_symmetric_difference(data.begin(), data.end(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data));
        return result;
    }
    May& operator+=(const May& rhs)
    {
        May tmp;
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data));
        return *this;
    }
    May operator*(const May& rhs) const;
    May mul(MMay rhs) const;
};

std::ostream& operator<<(std::ostream& sout, const May& x);

inline MMay divLF(MMay m1, MMay m2)
{
    return m1.divLF(m2);
}

inline bool divisibleLF(MMay m1, MMay m2)
{
    return m1.divisibleLF(m2);
}

inline MMay gcdLF(MMay m1, MMay m2)
{
    return m1.gcdLF(m2);
}

inline MMay lcmLF(MMay m1, MMay m2)
{
    return m1.lcmLF(m2);
}

/********************************************************
 *                      Modules
 ********************************************************/
/* Modules over A
*/
struct MMod
{
    MMay m;
    int v;
    bool operator<(MMod rhs) const
    {
        if (v > rhs.v)
            return true;
        if (v < rhs.v)
            return false;
        if (m < rhs.m)
            return true;
        return false;
    };
    bool operator==(MMod rhs) const
    {
        return m == rhs.m && v == rhs.v;
    };
    int deg(const array& basis_degs)
    {
        return m.deg() + basis_degs[v];
    }
};
using MMod1d = std::vector<MMod>;

struct Mod
{
    MMod1d data;
    Mod() {}
    Mod(MMod mv) : data({mv}) {}
    Mod(const May& a, int v) {
        for (MMay m : a.data)
            data.push_back(MMod{m, v});
    }

    MMod GetLead() const
    {
#ifndef NDEBUG
        if (data.empty())
            throw MyException(0x900cee0fU, "Mod empty");
#endif
        return data[0];
    }
    
    explicit operator bool() const
    {
        return !data.empty();
    }
    Mod operator+(const Mod& rhs) const
    {
        Mod result;
        std::set_symmetric_difference(data.begin(), data.end(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data));
        return result;
    }
    Mod& operator+=(const Mod& rhs)
    {
        Mod tmp;
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data));
        return *this;
    }
};
using Mod1d = std::vector<Mod>;

Mod operator*(const May& a, const Mod& x);
std::ostream& operator<<(std::ostream& sout, const Mod& x);

/**
 * Replace v_i with `map[i]`.
 */
inline Mod subs(const Mod& x, const Mod1d map)
{
    Mod result;
    for (const MMod& mv : x.data)
        result += May(mv.m) * map[mv.v];
    return result;
}

}  // namespace steenrod



#endif