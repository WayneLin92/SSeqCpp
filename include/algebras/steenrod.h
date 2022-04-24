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

/********************************************************
 *                  class Milnor
 ********************************************************/

inline constexpr size_t XI_MAX = 7;                     /* Support up to \xi_8 */
inline constexpr int DEG_MAX = (1 << (XI_MAX + 1)) - 2; /* Maximum degree supported in A */

inline constexpr size_t MMILNOR_M_BITS = (XI_MAX + 1) * (XI_MAX + 2) / 2 - 1;
inline constexpr size_t MMILNOR_W_BITS = 8;
inline constexpr uint64_t MMILNOR_ONE = uint64_t(1) << (MMILNOR_M_BITS - 1);
namespace detail {
    inline constexpr std::array<size_t, MMILNOR_M_BITS> MMilnorGenI()
    {
        std::array<size_t, MMILNOR_M_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_M_BITS)
                    result[n++] = i;
        return result;
    }
    inline constexpr std::array<size_t, MMILNOR_M_BITS> MMilnorGenJ()
    {
        std::array<size_t, MMILNOR_M_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_M_BITS)
                    result[n++] = j;
        return result;
    }
    inline constexpr std::array<int, MMILNOR_M_BITS> MMilnorGenDeg()
    {
        std::array<int, MMILNOR_M_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_M_BITS)
                    result[n++] = (1 << j) - (1 << i);
        return result;
    }
    inline constexpr std::array<uint64_t, MMILNOR_M_BITS> MMilnorGenWeight()
    {
        std::array<uint64_t, MMILNOR_M_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_M_BITS)
                    result[n++] = (2 * (j - i) - 1) << (64 - MMILNOR_W_BITS);
        return result;
    }
}  // namespace detail
inline constexpr std::array<size_t, MMILNOR_M_BITS> MMILNOR_GEN_I = detail::MMilnorGenI();
inline constexpr std::array<size_t, MMILNOR_M_BITS> MMILNOR_GEN_J = detail::MMilnorGenJ();
inline constexpr std::array<uint64_t, MMILNOR_M_BITS> MMILNOR_GEN_WEIGHT = detail::MMilnorGenWeight();
inline constexpr std::array<int, MMILNOR_M_BITS> MMILNOR_GEN_DEG = detail::MMilnorGenDeg();
inline constexpr uint64_t UINT64_LEFT_BIT = uint64_t(1) << 63;
inline constexpr uint64_t MMILNOR_MASK_M = (uint64_t(1) << MMILNOR_M_BITS) - 1;
inline constexpr uint64_t MMILNOR_MASK_W = ((uint64_t(1) << MMILNOR_W_BITS) - 1) << (64 - MMILNOR_W_BITS);
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
        for (uint64_t m = data << (64 - MMILNOR_M_BITS); m; m <<= 1, ++i)
            if (m & UINT64_LEFT_BIT)
                weight += MMILNOR_GEN_WEIGHT[i];
        return MMilnor(data + (weight << MMILNOR_M_BITS));
    }

public:
    /**
     * The first 10 bits are used to store the weight.
     * The rest are used to store the exterior monomial.
     */
    MMilnor() : data_(0) {}
    constexpr explicit MMilnor(uint64_t data) : data_(data) {}
    static MMilnor FromIndex(size_t i)
    {
        return MMilnor((MMILNOR_ONE >> i) | MMILNOR_GEN_WEIGHT[i]);
    }
    static uint64_t RawP(int i, int j)
    {
        size_t index = size_t(j * (j + 1) / 2 - i - 1);
        return (MMILNOR_ONE >> index) + MMILNOR_GEN_WEIGHT[index];
    }
    static MMilnor P(int i, int j)
    {
        return MMilnor(MMilnor::RawP(i, j));
    }
    static MMilnor Xi(const int* xi)
    {
        uint64_t result = 0;
        for (int d = 1; d <= XI_MAX; ++d)
            for (int n = xi[size_t(d - 1)], i = 0; n; n >>= 1, ++i)
                if (n & 1)
                    result += MMilnor::RawP(i, i + d);
        return MMilnor(result);
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
        return MMilnor(((data_ | rhs.data_) & MMILNOR_MASK_M) | ((data_ & MMILNOR_MASK_W) + (rhs.data_ & MMILNOR_MASK_W)));
    }

    MMilnor divLF(MMilnor rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        if (!rhs.divisibleLF(*this))
            throw MyException(0x7ed0a8U, "Not divisible: m1 / m2");
#endif
        return MMilnor(((data_ ^ rhs.data_) & MMILNOR_MASK_M) | ((data_ & MMILNOR_MASK_W) - (rhs.data_ & MMILNOR_MASK_W)));
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

    uint64_t weight() const
    {
        return data_ >> (64 - MMILNOR_W_BITS);
    }

    int deg() const
    {
        int result = 0;
        for (int i : *this)
            result += MMILNOR_GEN_DEG[i];
        return result;
    }

    std::string StrXi() const;

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
                if (m_ & UINT64_LEFT_BIT)
                    break;
            return *this;
        }
        iterator operator++(int)
        {
            iterator copy(*this);
            for (; m_; m_ <<= 1, ++i_)
                if (m_ & UINT64_LEFT_BIT)
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
                if (m_ & UINT64_LEFT_BIT)
                    break;
        }
    };

    iterator begin() const
    {
        return iterator(data_ << (64 - MMILNOR_M_BITS));
    }
    iterator end() const
    {
        return iterator(0);
    }
};
using MMilnor1d = std::vector<MMilnor>;
using MMilnor2d = std::vector<MMilnor1d>;

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

    std::string StrXi() const;
};

std::ostream& operator<<(std::ostream& sout, const Milnor& x);
inline Milnor operator*(MMilnor m1, MMilnor m2) ////
{
    return Milnor(m1) * Milnor(m2);
}

/********************************************************
 *                    class Mod
 ********************************************************/

inline constexpr uint64_t MMOD_MASK_M = MMILNOR_MASK_W + MMILNOR_MASK_M;
inline constexpr uint64_t MMOD_MASK_V = ~MMOD_MASK_M;
/* Modules over A */
class MMod
{
private:
    uint64_t data_;

public:
    MMod() : data_(0) {}
    constexpr explicit MMod(uint64_t data) : data_(data) {}
    MMod(MMilnor m, uint64_t v) : data_(m.data() + ((~v) << MMILNOR_M_BITS)) {}

    static MMod FromRaw(MMilnor m, uint64_t v)
    {
        return MMod(m.data() + v);
    }

    MMilnor m() const
    {
        return MMilnor(data_ & MMOD_MASK_M);
    }
    uint64_t v() const
    {
        return (~data_ & MMOD_MASK_V) >> MMILNOR_M_BITS;
    }
    uint64_t v_raw() const
    {
        return data_ & MMOD_MASK_V;
    }

    bool operator==(MMod rhs) const
    {
        return data_ == rhs.data_;
    };
    bool operator<(MMod rhs) const
    {
        //return data_ < rhs.data_;
        auto vr1 = v_raw(), vr2 = rhs.v_raw();
        if (vr1 < vr2)
            return true;
        if (vr1 > vr2)
            return false;
        if (data_ < rhs.data_)
            return true;
        return false;
    };
    explicit operator bool() const
    {
        return data_;
    }
    std::string Str() const;
    std::string StrXi() const;
};
using MMod1d = std::vector<MMod>;
using MMod2d = std::vector<MMod1d>;

inline MMod mulLF(MMilnor m, MMod x)
{
    return MMod::FromRaw(mulLF(m, x.m()), x.v_raw());
}

struct Mod
{
    MMod1d data;
    Mod() {}
    Mod(MMod mv) : data({mv}) {}

    MMod GetLead() const
    {
#ifndef NDEBUG
        if (data.empty())
            throw MyException(0x900cee0fU, "Trying to GetLead() for empty Mod.");
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
    /* `*this += m * x` */
    Mod& iaddmul(MMilnor m, const Mod& x, Milnor& tmp_a, Mod& tmp_m1, Mod& tmp_m2);
    bool operator==(const Mod& rhs) const
    {
        return data == rhs.data;
    };

    std::string StrXi() const;
    std::string Str() const;
};
using Mod1d = std::vector<Mod>;
using Mod2d = std::vector<Mod1d>;
using Mod3d = std::vector<Mod2d>;

Mod operator*(MMilnor m, const Mod& x);

inline std::ostream& operator<<(std::ostream& sout, const Mod& x)
{
    return sout << x.Str();
}

/* Compute the product in the associated graded algebra */
Mod mulLF(MMilnor m, const Mod& x);

template <typename FnMap>
Mod TplSubs(const Mod& x, FnMap map)
{
    Mod result;
    for (const MMod& mv : x.data)
        result += Milnor(mv.m()) * map(mv.v());
    return result;
}

/**
 * Replace v_i with `map[i]`.
 */
//inline Mod subs(const Mod& x, const Mod1d& map)
//{
//    Mod result;
//    for (const MMod& mv : x.data)
//        result += Milnor(mv.m()) * map[mv.v()];
//    return result;
//}

}  // namespace steenrod

#endif