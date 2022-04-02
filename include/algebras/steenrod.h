/** \file steenrod.h
 * A Component for monomials and polynomials over $F_2$.
 * All are encapsulated in the `namespace alg`.
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

inline constexpr size_t XI_MAX = 8;                     /* Support up to \xi_8 */
inline constexpr int DEG_MAX = (1 << (XI_MAX + 1)) - 1; /* Maximum degree supported in A */

inline constexpr size_t MMILNOR_INDEX_NUM = (XI_MAX + 1) * (XI_MAX + 2) / 2 - 1;
inline constexpr uint64_t MMILNOR_ONE = uint64_t(1) << (MMILNOR_INDEX_NUM - 1);
namespace detail {
    inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMilnorGenI()
    {
        std::array<int, MMILNOR_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMILNOR_INDEX_NUM)
                    result[n++] = i;
        return result;
    }
    inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMilnorGenJ()
    {
        std::array<int, MMILNOR_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMILNOR_INDEX_NUM)
                    result[n++] = j;
        return result;
    }
    inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMilnorGenDeg()
    {
        std::array<int, MMILNOR_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMILNOR_INDEX_NUM)
                    result[n++] = (1 << j) - (1 << i);
        return result;
    }
    inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMilnorGenWeight()
    {
        std::array<int, MMILNOR_INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)XI_MAX + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < MMILNOR_INDEX_NUM)
                    result[n++] = 2 * (j - i) - 1;
        return result;
    }
}  // namespace detail
inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMILNOR_GEN_I = detail::MMilnorGenI();
inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMILNOR_GEN_J = detail::MMilnorGenJ();
inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMILNOR_GEN_DEG = detail::MMilnorGenDeg();
inline constexpr std::array<int, MMILNOR_INDEX_NUM> MMILNOR_GEN_WEIGHT = detail::MMilnorGenWeight();
inline constexpr uint64_t MMILNOR_LEFT_BIT = uint64_t(1) << 63;
inline constexpr uint64_t MMILNOR_MASK_M = (uint64_t(1) << MMILNOR_INDEX_NUM) - 1;
inline constexpr uint64_t MMILNOR_MASK_W = ~MMILNOR_MASK_M;
inline constexpr uint64_t MMILNOR_NULL = 0xffffffffffffffff;

/** Milnor basis for A ordered by May filtration $w(xi_j^{2^i})=2j-1$
 *
 * Each element is represented by a 64-bit unsigned integer
 */

class MMilnor
{
private:
    uint64_t data_;

private:
    static MMilnor AddWeight(uint64_t data)
    {
        uint64_t weight = 0;
        int i = 0;
        for (uint64_t m = data << (64 - MMILNOR_INDEX_NUM); m; m <<= 1, ++i)
            if (m & MMILNOR_LEFT_BIT)
                weight += uint64_t(MMILNOR_GEN_WEIGHT[i]);
        return MMilnor(data + (weight << MMILNOR_INDEX_NUM));
    }

public:
    /**
     * The first 10 bits are used to store the weight.
     * The rest are used to store the exterior monomial.
     */
    MMilnor() : data_(0) {}
    explicit MMilnor(uint64_t data) : data_(data) {}
    static MMilnor FromIndex(size_t i)
    {
        return MMilnor((MMILNOR_ONE >> i) + (uint64_t(MMILNOR_GEN_WEIGHT[i]) << MMILNOR_INDEX_NUM));
    }
    static MMilnor P(int i, int j)
    {
        return MMilnor((MMILNOR_ONE >> (j * (j + 1) / 2 - i - 1)) + (uint64_t(2 * (j - i) - 1) << MMILNOR_INDEX_NUM));
    }
    static uint64_t rawP(int i, int j)
    {
        return MMILNOR_ONE >> (j * (j + 1) / 2 - i - 1);
    }
    static MMilnor Xi(const std::array<int, XI_MAX>& xi)
    {
        uint64_t result = 0;
        uint64_t weight = 0;
        for (int d = 1; d <= xi.size(); ++d) {
            for (int n = xi[size_t(d - 1)], i = 0; n; n >>= 1, ++i) {
                if (n & 1) {
                    int j = i + d;
                    result |= MMilnor::rawP(i, j);
                    weight += 2 * uint64_t(d) - 1;
                }
            }
        }
        return MMilnor(result + (weight << MMILNOR_INDEX_NUM));
    }

public:
    std::array<int, XI_MAX> ToXi() const
    {
        std::array<int, XI_MAX> result = {};
        for (int i : *this)
            result[size_t(MMILNOR_GEN_J[i] - MMILNOR_GEN_I[i] - 1)] += 1 << MMILNOR_GEN_I[i];
        return result;
    }

    bool operator<(MMilnor rhs) const
    {
        return data_ < rhs.data_;
    };

    bool operator==(MMilnor rhs) const
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

    MMilnor mulLF(MMilnor rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        if (gcdLF(rhs))
            throw MyException(0x170c454aU, "gcd(m1,m2)!=1");
#endif
        return MMilnor(((data_ | rhs.data_) & MMILNOR_MASK_M) + ((data_ & MMILNOR_MASK_W) + (rhs.data_ & MMILNOR_MASK_W)));
    }

    MMilnor divLF(MMilnor rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        if (!rhs.divisibleLF(*this))
            throw MyException(0x7ed0a8U, "Not divisible: m1 / m2");
#endif
        return MMilnor(((data_ ^ rhs.data_) & MMILNOR_MASK_M) + ((data_ & MMILNOR_MASK_W) - (rhs.data_ & MMILNOR_MASK_W)));
    }

    bool divisibleLF(MMilnor rhs)
    {
        uint64_t m1 = data_ & MMILNOR_MASK_M;
        uint64_t m2 = rhs.data_ & MMILNOR_MASK_M;
        return m2 >= m1 && !(m1 & (m2 - m1));
    }

    MMilnor gcdLF(MMilnor rhs) const
    {
        return AddWeight(data_ & rhs.data_ & MMILNOR_MASK_M);
    }

    MMilnor lcmLF(MMilnor rhs) const
    {
        return AddWeight((data_ | rhs.data_) & MMILNOR_MASK_M);
    }

    int weight() const
    {
        return (int)(data_ >> MMILNOR_INDEX_NUM);
    }

    int deg() const
    {
        int result = 0;
        for (int i : *this)
            result += MMILNOR_GEN_DEG[i];
        return result;
    }

    std::string StrXi();

public:
    class iterator
    {
        friend class MMilnor;

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
                if (m_ & MMILNOR_LEFT_BIT)
                    break;
            return *this;
        }
        iterator operator++(int)
        {
            iterator copy(*this);
            for (; m_; m_ <<= 1, ++i_)
                if (m_ & MMILNOR_LEFT_BIT)
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
                if (m_ & MMILNOR_LEFT_BIT)
                    break;
        }
    };

    iterator begin() const
    {
        return iterator(data_ << (64 - MMILNOR_INDEX_NUM));
    }
    iterator end() const
    {
        return iterator(0);
    }
};
using MMilnor1d = std::vector<MMilnor>;
using MMilnor2d = std::vector<MMilnor1d>;

/* Elements of A as linear combinations of Milnor basis
 */
struct Milnor
{
    MMilnor1d data;
    Milnor() {}
    explicit Milnor(MMilnor m) : data({m}) {}

    static Milnor P(int i, int j)
    {
        return Milnor(MMilnor::P(i, j));
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
    Milnor operator*(const Milnor& rhs) const;
    Milnor mul(MMilnor rhs) const;
};

std::ostream& operator<<(std::ostream& sout, const Milnor& x);

inline MMilnor mulLF(MMilnor m1, MMilnor m2)
{
    return m1.mulLF(m2);
}

inline MMilnor divLF(MMilnor m1, MMilnor m2)
{
    return m1.divLF(m2);
}

inline bool divisibleLF(MMilnor m1, MMilnor m2)
{
    return m1.divisibleLF(m2);
}

inline MMilnor gcdLF(MMilnor m1, MMilnor m2)
{
    return m1.gcdLF(m2);
}

inline MMilnor lcmLF(MMilnor m1, MMilnor m2)
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
    MMilnor m;
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
    Mod(const Milnor& a, int v)
    {
        for (MMilnor m : a.data)
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
Mod operator*(const Milnor& a, const Mod& x);

template <typename FnMap>
Mod TplSubs(const Mod& x, FnMap map)
{
    Mod result;
    for (const MMod& mv : x.data)
        result += Milnor(mv.m) * map(mv.v);
    return result;
}

/**
 * Replace v_i with `map[i]`.
 */
inline Mod subs(const Mod& x, const Mod1d& map)
{
    Mod result;
    for (const MMod& mv : x.data)
        result += Milnor(mv.m) * map[mv.v];
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
    MModCpt(MMilnor m, int v) : data_(m.data() + ((~uint64_t(v)) << (64 - MMOD_BASIS_BITS))) {}

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
    MMilnor m() const
    {
        return MMilnor(data_ & MMOD_MASK_M);
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
    ModCpt(const Milnor& a, int v)
    {
        for (MMilnor m : a.data)
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
ModCpt operator*(const Milnor& a, const ModCpt& x);
/* Compute the product in the associated graded algebra */
ModCpt mulLF(MMilnor m, const ModCpt& x);

void MulMilnorV3(MMilnor lhs, MMilnor rhs, Milnor& result);  ////

}  // namespace steenrod

#endif