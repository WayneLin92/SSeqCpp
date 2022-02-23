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
using MMay = uint64_t;
struct May;

inline constexpr size_t MILNOR_BUFFER_SIZE = 9;                     /* Support up to \xi_8 */
inline constexpr int DEG_MAX = (1 << (MILNOR_BUFFER_SIZE + 1)) - 1; /* Maximum degree supported in A */

struct MMilnor
{
    using TypeData = std::array<int, MILNOR_BUFFER_SIZE>;
    TypeData data;

    MMilnor() : data({}) {}
    explicit MMilnor(MMay m);

    static MMilnor P(int i, int j) /* P_{ij} */
    {
        MMilnor result{{}};
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

    /* May leading form according to weight */
    MMay ToMMay() const;
};
using MonMilnor1d = std::vector<MMilnor>;

struct Milnor
{
    MonMilnor1d data;

    Milnor() {}
    Milnor(const MMilnor& r) : data({r}) {}
    explicit Milnor(MMay m) : data({MMilnor(m)}) {}
    explicit Milnor(const May& x);

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

inline constexpr size_t MAY_INDEX_NUM = (MILNOR_BUFFER_SIZE + 1) * (MILNOR_BUFFER_SIZE + 2) / 2 - 1;
static_assert(MAY_INDEX_NUM < sizeof(uint64_t) * 8);

namespace detail {
    inline constexpr std::array<int, MAY_INDEX_NUM> MonSteenrodMayI()
    {
        std::array<int, MAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)MILNOR_BUFFER_SIZE + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MAY_INDEX_NUM)
                    result[n++] = i;
        return result;
    }
    inline constexpr std::array<int, MAY_INDEX_NUM> MonSteenrodMayJ()
    {
        std::array<int, MAY_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)MILNOR_BUFFER_SIZE + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MAY_INDEX_NUM)
                    result[n++] = j;
        return result;
    }
}  // namespace detail

/* May basis for A
** 
** Each element is represented by a 64-bit unsigned integer
*/
using MMay1d = std::vector<MMay>;
using MMay2d = std::vector<MMay1d>;
inline constexpr std::array<int, MAY_INDEX_NUM> MAY_GEN_I = detail::MonSteenrodMayI();
inline constexpr std::array<int, MAY_INDEX_NUM> MAY_GEN_J = detail::MonSteenrodMayJ();

inline constexpr MMay MMayFromIndex(int index)
{
    return MMay(1) << index;
}

inline constexpr int MayIndexP(int i, int j)
{
#ifndef NDEBUG
    if (!(0 <= i && i < j))
        throw MyException(0x4bfc0d6fU, "Invalid input");
#endif  // !NDEBUG
    return j * (j + 1) / 2 - i - 1;
}

inline constexpr MMay MMayP(int i, int j) /* P_{ij} */
{
    return MMayFromIndex(MayIndexP(i, j));
}

inline constexpr int Weight(MMay m)
{
    int result = 0;
    int i = 0;
    while (m) {
        if (m & 1)
            result += 2 * (MAY_GEN_J[i] - MAY_GEN_I[i]) - 1;
        m >>= 1;
        ++i;
    }
    return result;
}

inline constexpr int GetDeg(MMay m)
{
    int result = 0;
    int i = 0;
    while (m) {
        if (m & 1)
            result += (1 << MAY_GEN_J[i]) - (1 << MAY_GEN_I[i]);
        m >>= 1;
        ++i;
    }
    return result;
}

inline bool CmpMMay(MMay lhs, MMay rhs)
{
    int w1 = Weight(lhs), w2 = Weight(rhs);
    if (w1 < w2)
        return true;
    if (w2 < w1)
        return false;
    if (ut::Reverse(lhs) < ut::Reverse(rhs)) /* We want Revlex here */
        return true;
    return false;
};


/* Elements of A as linear combinations of May basis 
*/
struct May
{
    MMay1d data;
    May() {}
    explicit May(MMay m) : data({m}) {}

    /* Convert from Milnor basis to May basis */
    explicit May(const Milnor& x)
    {
        for (auto& m : x.data)
            data.push_back(m.ToMMay());
        std::sort(data.begin(), data.end(), CmpMMay);
    }

    static May P(int i, int j)
    {
        return May(MMayP(i, j));
    }

    May operator+(const May& rhs) const
    {
        May result;
        std::set_symmetric_difference(data.begin(), data.end(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data), CmpMMay);
        return result;
    }
    May& operator+=(const May& rhs)
    {
        May tmp;
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data), CmpMMay);
        return *this;
    }
    May mul(MMay rhs) const
    {
        return May(Milnor(*this) * Milnor(rhs));
    }
    May operator*(const May& rhs) const
    {
        return May(Milnor(*this) * Milnor(rhs));
    }
};

std::ostream& operator<<(std::ostream& sout, const May& x);

inline MMay divLF(MMay m1, MMay m2)
{
#ifndef NDEBUG /* DEBUG */
    if (m1 < m2 || m2 & (m1 - m2))
        throw MyException(0x7ed0a8U, "Not divisible: m1 / m2");
#endif
    return m1 ^ m2;
}

inline bool divisibleLF(MMay m1, MMay m2)
{
    return m2 > m1 && !(m1 & (m2 - m1));
}

inline MMay gcdLF(MMay m1, MMay m2)
{
    return m1 & m2;
}

inline MMay lcmLF(MMay m1, MMay m2)
{
    return m1 | m2;
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
        if (CmpMMay(m, rhs.m))
            return true;
        return false;
    };
    bool operator==(MMod rhs) const
    {
        return m == rhs.m && v == rhs.v;
    };
    int deg(const array& basis_degs)
    {
        return GetDeg(m) + basis_degs[v];
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