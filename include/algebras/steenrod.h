/** \file steenrod.h
 * A Component for monomials and polynomials over $F_2$.
 * All are incapsulated in the `namespace alg`.
 */

#ifndef STEENROD_H
#define STEENROD_H

#include "myexception.h"
#include "utility.h"
#include "benchmark.h" ////
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

inline constexpr size_t XI_MAX = 8;                     /* Support up to \xi_8 */
inline constexpr int DEG_MAX = (1 << (XI_MAX + 1)) - 1; /* Maximum degree supported in A */

struct MMilnor
{
    using TypeData = std::array<int, XI_MAX>;
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
        for (size_t i = 0; i < XI_MAX; ++i)
            result += (2 * (int)i + 1) * ut::popcount(unsigned(data[i]));
        return result;
    }

    int deg() const
    {
        int result = 0;
        for (size_t i = 0; i < XI_MAX; ++i)
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


/** May basis for A
 *
 * Each element is represented by a 64-bit unsigned integer
*/

inline constexpr size_t MMAY_INDEX_NUM = (XI_MAX + 1) * (XI_MAX + 2) / 2 - 1;
inline constexpr uint64_t MMAY_ONE = uint64_t(1) << (MMAY_INDEX_NUM - 1);
namespace detail {
    inline constexpr std::array<int, MMAY_INDEX_NUM> MMayGenI()
    {
        std::array<int, MMAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMAY_INDEX_NUM)
                    result[n++] = i;
        return result;
    }
    inline constexpr std::array<int, MMAY_INDEX_NUM> MMayGenJ()
    {
        std::array<int, MMAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMAY_INDEX_NUM)
                    result[n++] = j;
        return result;
    }
    inline constexpr std::array<int, MMAY_INDEX_NUM> MMayGenDeg()
    {
        std::array<int, MMAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMAY_INDEX_NUM)
                    result[n++] = (1 << j) - (1 << i);
        return result;
    }
    inline constexpr std::array<int, MMAY_INDEX_NUM> MMayGenWeight()
    {
        std::array<int, MMAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMAY_INDEX_NUM)
                    result[n++] = 2 * (j - i) - 1;
        return result;
    }
}  // namespace detail
inline constexpr std::array<int, MMAY_INDEX_NUM> MMAY_GEN_I = detail::MMayGenI();
inline constexpr std::array<int, MMAY_INDEX_NUM> MMAY_GEN_J = detail::MMayGenJ();
inline constexpr std::array<int, MMAY_INDEX_NUM> MMAY_GEN_DEG = detail::MMayGenDeg();
inline constexpr std::array<int, MMAY_INDEX_NUM> MMAY_GEN_WEIGHT = detail::MMayGenWeight();
inline constexpr uint64_t MMAY_LEFT_BIT = uint64_t(1) << 63;
inline constexpr uint64_t MMAY_MASK_M = (uint64_t(1) << MMAY_INDEX_NUM) - 1;
inline constexpr uint64_t MMAY_MASK_W = ~MMAY_MASK_M;
inline constexpr uint64_t MMAY_NULL = 0xffffffffffffffff;

class MMay
{
private:
    uint64_t data_;

private:
    static MMay AddWeight(uint64_t data)
    {
        uint64_t weight = 0;
        int i = 0;
        for (uint64_t m = data << (64 - MMAY_INDEX_NUM); m; m <<= 1, ++i)
            if (m & MMAY_LEFT_BIT)
                weight += uint64_t(MMAY_GEN_WEIGHT[i]);
        return MMay(data + (weight << MMAY_INDEX_NUM));
    }

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
    static MMay FromIndex(size_t i)
    {
        return MMay((MMAY_ONE >> i) + (uint64_t(MMAY_GEN_WEIGHT[i]) << MMAY_INDEX_NUM));
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
    uint64_t data() const
    {
        return data_;
    }

    MMay mulLF(MMay rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        if (gcdLF(rhs))
            throw MyException(0x170c454aU, "gcd(m1,m2)!=1");
#endif
        return MMay(((data_ | rhs.data_) & MMAY_MASK_M) + ((data_ & MMAY_MASK_W) + (rhs.data_ & MMAY_MASK_W)));
    }

    MMay divLF(MMay rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        if (!rhs.divisibleLF(*this))
            throw MyException(0x7ed0a8U, "Not divisible: m1 / m2");
#endif
        return MMay(((data_ ^ rhs.data_) & MMAY_MASK_M) + ((data_ & MMAY_MASK_W) - (rhs.data_ & MMAY_MASK_W)));
    }

    bool divisibleLF(MMay rhs)
    {
        uint64_t m1 = data_ & MMAY_MASK_M;
        uint64_t m2 = rhs.data_ & MMAY_MASK_M;
        return m2 >= m1 && !(m1 & (m2 - m1));
    }

    MMay gcdLF(MMay rhs) const
    {
        return AddWeight(data_ & rhs.data_ & MMAY_MASK_M);
    }

    MMay lcmLF(MMay rhs) const
    {
        return AddWeight((data_ | rhs.data_) & MMAY_MASK_M);
    }

    int weight() const
    {
        return (int)(data_ >> MMAY_INDEX_NUM);
    }

    int deg() const
    {
        int result = 0;
        for (int i : *this)
            result += MMAY_GEN_DEG[i];
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

inline MMay mulLF(MMay m1, MMay m2)
{
    return m1.mulLF(m2);
}

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
    Mod(const May& a, int v)
    {
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
    bool operator==(const Mod& rhs) const
    {
        return data == rhs.data;
    };
};
using Mod1d = std::vector<Mod>;
using Mod2d = std::vector<Mod1d>;

std::ostream& operator<<(std::ostream& sout, const Mod& x);
Mod operator*(const May& a, const Mod& x);

template <typename FnMap>
Mod TplSubs(const Mod& x, FnMap map)
{
    Mod result;
    for (const MMod& mv : x.data)
        result += May(mv.m) * map(mv.v);
    return result;
}

/**
 * Replace v_i with `map[i]`.
 */
inline Mod subs(const Mod& x, const Mod1d& map)
{
    Mod result;
    for (const MMod& mv : x.data)
        result += May(mv.m) * map[mv.v];
    return result;
}

/********************************************************
 *          Modules with Compact data structure
 ********************************************************/

/* The left 12 bits will be used to store the basis */
inline constexpr unsigned MMOD_BASIS_BITS = 12;
inline constexpr uint64_t MMOD_MASK_M = (uint64_t(1) << (64 - MMOD_BASIS_BITS)) - 1;
inline constexpr uint64_t MMOD_MASK_V = ~MMOD_MASK_M;
/* Modules over A */
class MModCpt
{
private:
    uint64_t data_;

public:
    MModCpt() : data_(0) {}
    MModCpt(uint64_t data) : data_(data) {}
    MModCpt(MMay m, int v) : data_(m.data() + ((~uint64_t(v)) << (64 - MMOD_BASIS_BITS))) {}

    bool operator<(MModCpt rhs) const
    {
        return data_ < rhs.data_;
    };
    bool operator==(MModCpt rhs) const
    {
        return data_ == rhs.data_;
    };
    explicit operator bool() const
    {
        return data_;
    }
    /*int deg(const array& basis_degs)
    {
        return MMay(data_).deg() + basis_degs[(~data_) >> (64 - MMOD_BASIS_BITS)];
    }*/
    MMay m() const
    {
        return MMay(data_ & MMOD_MASK_M);
    }
    int v() const
    {
        return int((~data_) >> (64 - MMOD_BASIS_BITS));
    }
};
using MModCpt1d = std::vector<MModCpt>;
using MModCpt2d = std::vector<MModCpt1d>;

struct ModCpt
{
    MModCpt1d data;
    ModCpt() {}
    ModCpt(MModCpt mv) : data({mv}) {}
    ModCpt(const May& a, int v)
    {
        for (MMay m : a.data)
            data.push_back(MModCpt(m, v));
    }

    MModCpt GetLead() const
    {
#ifndef NDEBUG
        if (data.empty())
            throw MyException(0x900cee0fU, "ModCpt empty");
#endif
        return data[0];
    }

    explicit operator bool() const
    {
        return !data.empty();
    }
    ModCpt operator+(const ModCpt& rhs) const
    {
        ModCpt result;
        std::set_symmetric_difference(data.begin(), data.end(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data));
        return result;
    }
    ModCpt& operator+=(const ModCpt& rhs)
    {
        ModCpt tmp;
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data));
        return *this;
    }
    bool operator==(const ModCpt& rhs) const
    {
        return data == rhs.data;
    };
};
using ModCpt1d = std::vector<ModCpt>;
using ModCpt2d = std::vector<ModCpt1d>;
using ModCpt3d = std::vector<ModCpt2d>;

std::ostream& operator<<(std::ostream& sout, const ModCpt& x);
ModCpt operator*(const May& a, const ModCpt& x);
/* Compute the product in the associated graded algebra */
ModCpt mulLF(MMay m, const ModCpt& x);

void MulMilnorV3(MMay lhs, MMay rhs, May& result);////

}  // namespace steenrod

#endif