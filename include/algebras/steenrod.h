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

inline constexpr size_t XI_MAX = 8;      /* Support up to \xi_8 */
inline constexpr size_t XI_MAX_MULT = 7; /* Multiplication support up to \xi_7 */

inline constexpr size_t MMILNOR_E_BITS = 37;
inline constexpr size_t MMILNOR_W_BITS = 9;
inline constexpr uint64_t MMILNOR_ONE = uint64_t(1) << (MMILNOR_E_BITS - 1);
namespace detail {
    inline constexpr std::array<size_t, MMILNOR_E_BITS> MMilnorGenI()
    {
        std::array<size_t, MMILNOR_E_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_E_BITS)
                    result[n++] = i;
        return result;
    }
    inline constexpr std::array<size_t, MMILNOR_E_BITS> MMilnorGenJ()
    {
        std::array<size_t, MMILNOR_E_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_E_BITS)
                    result[n++] = j;
        return result;
    }
    inline constexpr std::array<int, MMILNOR_E_BITS> MMilnorGenDeg()
    {
        std::array<int, MMILNOR_E_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_E_BITS)
                    result[n++] = (1 << j) - (1 << i);
        return result;
    }
    inline constexpr std::array<uint64_t, MMILNOR_E_BITS> MMilnorGenWeight()
    {
        std::array<uint64_t, MMILNOR_E_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_E_BITS)
                    result[n++] = (2 * (j - i) - 1) << MMILNOR_E_BITS;
        return result;
    }
}  // namespace detail
inline constexpr std::array<size_t, MMILNOR_E_BITS> MMILNOR_GEN_I = detail::MMilnorGenI();
inline constexpr std::array<size_t, MMILNOR_E_BITS> MMILNOR_GEN_J = detail::MMilnorGenJ();
inline constexpr std::array<uint64_t, MMILNOR_E_BITS> MMILNOR_GEN_WEIGHT = detail::MMilnorGenWeight();
inline constexpr std::array<int, MMILNOR_E_BITS> MMILNOR_GEN_DEG = detail::MMilnorGenDeg();
inline constexpr uint64_t UINT64_LEFT_BIT = uint64_t(1) << 63;
inline constexpr uint64_t MMILNOR_MASK_E = (uint64_t(1) << MMILNOR_E_BITS) - 1;
inline constexpr uint64_t MMILNOR_MASK_W = ((uint64_t(1) << MMILNOR_W_BITS) - 1) << MMILNOR_E_BITS;
inline constexpr uint64_t MMILNOR_NULL = ~uint64_t(0);

/* Maximum degree supported in A if `MMILNOR_E_BITS == 37` */
inline constexpr int DEG_MAX = 383;
inline constexpr int DEG_MAX_MULT = (1 << (XI_MAX_MULT + 1)) - 2;

/** Milnor basis for A ordered by May filtration $w(xi_j^{2^i})=2j-1$
 *
 * Each element is represented by a 64-bit unsigned integer
 * (0, w, e): (19 bits, 9 bits, 37 bits)
 */

class MMilnor
{
private:
    uint64_t data_;

public:
    MMilnor() : data_(0) {}
    constexpr explicit MMilnor(uint64_t data) : data_(data) {}

    static MMilnor FromIndex(size_t i)
    {
        return MMilnor((MMILNOR_ONE >> i) | MMILNOR_GEN_WEIGHT[i]);
    }

    static uint64_t dataP(int i, int j)
    {
        size_t index = size_t(j * (j + 1) / 2 - i - 1);
        return (MMILNOR_ONE >> index) | MMILNOR_GEN_WEIGHT[index];
    }

    static MMilnor P(int i, int j)
    {
        return MMilnor(MMilnor::dataP(i, j));
    }

    static MMilnor Xi(const uint32_t* xi)
    {
        uint64_t result = 0;
        for (int d = 1; d <= XI_MAX_MULT; ++d)  ////
            for (uint32_t n = xi[size_t(d - 1)], i = 0; n; n >>= 1, ++i)
                if (n & 1)
                    result += MMilnor::dataP(i, i + d);  ////
        return MMilnor(result);
    }

    static MMilnor FromE(uint64_t e)
    {
        return MMilnor(e + WRaw(e));
    }

    static uint64_t WRaw(uint64_t e)
    {
        uint64_t w_raw = 0;
        int i = 0;
        for (uint64_t b = e << (64 - MMILNOR_E_BITS); b; b <<= 1, ++i)
            if (b & UINT64_LEFT_BIT)
                w_raw += MMILNOR_GEN_WEIGHT[i];
        return w_raw;
    }

public:
    std::array<uint32_t, XI_MAX> ToXi() const
    {
        std::array<uint32_t, XI_MAX> result = {};
        for (int i : *this)
            result[size_t(MMILNOR_GEN_J[i] - MMILNOR_GEN_I[i] - 1)] += uint32_t(1) << MMILNOR_GEN_I[i];
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
    uint64_t w_raw() const
    {
        return data_ & MMILNOR_MASK_W;
    }
    uint64_t w() const
    {
        return data_ >> MMILNOR_E_BITS;
    }
    uint64_t e() const
    {
        return data_ & MMILNOR_MASK_E;
    }
    uint64_t w_may() const
    {
        return (w() + ut::popcount(e())) / 2;
    }

    MMilnor mulLF(MMilnor rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        if (gcdLF(rhs))
            throw MyException(0x170c454aU, "gcd(m1,m2)!=1");
#endif
        return MMilnor(((data_ | rhs.data_) & MMILNOR_MASK_E) | (w_raw() + rhs.w_raw()));
    }

    MMilnor divLF(MMilnor rhs) const
    {
#ifndef NDEBUG /* DEBUG */
        if (!rhs.divisibleLF(*this))
            throw MyException(0x7ed0a8U, "Not divisible: m1 / m2");
#endif
        return MMilnor(((data_ ^ rhs.data_) & MMILNOR_MASK_E) | (w_raw() - rhs.w_raw()));
    }

    bool divisibleLF(MMilnor rhs)
    {
        uint64_t m1 = e();
        uint64_t m2 = rhs.e();
        return m2 >= m1 && !(m1 & (m2 - m1));
    }

    MMilnor gcdLF(MMilnor rhs) const
    {
        return FromE(data_ & rhs.data_ & MMILNOR_MASK_E);
    }

    MMilnor lcmLF(MMilnor rhs) const
    {
        return FromE((data_ | rhs.data_) & MMILNOR_MASK_E);
    }

    int deg() const
    {
        int result = 0;
        for (int i : *this)
            result += MMILNOR_GEN_DEG[i];
        return result;
    }

    std::string StrXi() const;
    std::string Str() const;

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
        return iterator(data_ << (64 - MMILNOR_E_BITS));
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

inline std::ostream& operator<<(std::ostream& sout, MMilnor m)
{
    return std::cout << m.Str();
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
    std::string Str() const;
};

inline std::ostream& operator<<(std::ostream& sout, const Milnor& x)
{
    return std::cout << x.Str();
}
inline Milnor operator*(MMilnor m1, MMilnor m2)  ////
{
    return Milnor(m1) * Milnor(m2);
}

/********************************************************
 *                    class Mod
 ********************************************************/

inline constexpr uint64_t MMOD_M_BITS = MMILNOR_W_BITS + MMILNOR_E_BITS;
inline constexpr uint64_t MMOD_MASK_V = ~(MMILNOR_MASK_W | MMILNOR_MASK_E);
/* Modules over A */
class MMod
{
private:
    uint64_t data_;

public:
    constexpr MMod() : data_(MMOD_MASK_V) {}
    constexpr explicit MMod(uint64_t data) : data_(data) {}
    explicit MMod(MMilnor m, uint64_t v) : data_(m.data() | (~v << MMOD_M_BITS)) {}

    uint64_t data() const
    {
        return data_;
    }
    uint64_t w_raw() const
    {
        return data_ & MMILNOR_MASK_W;
    }
    uint64_t v_raw() const
    {
        return data_ & MMOD_MASK_V;
    }
    uint64_t e() const
    {
        return data_ & MMILNOR_MASK_E;
    }
    MMilnor m_no_weight() const
    {
        return MMilnor(e());
    }
    uint64_t w() const
    {
        return w_raw() >> MMILNOR_E_BITS;
    }
    uint64_t v() const
    {
        return ~data_ >> MMOD_M_BITS;
    }
    MMilnor m() const
    {
        return MMilnor::FromE(e());
    }
    uint64_t w_may() const
    {
        return (w() + ut::popcount(e())) / 2;
    }

    bool operator==(MMod rhs) const
    {
        return data_ == rhs.data_;
    };
    bool operator<(MMod rhs) const
    {
        return data_ < rhs.data_;
    };
    explicit operator bool() const
    {
        return data_;
    }
    /* If lhs.v == rhs.v, use this function to test divisibility. */
    bool divisibleLF(MMod rhs)
    {
#ifndef NDEBUG /* DEBUG */
        if (v_raw() != rhs.v_raw())
            throw MyException(0x3bc29cceU, "use MMod::divisible only when the two v's agree.");
#endif
        uint64_t e1 = e();
        uint64_t e2 = rhs.e();
        return e2 >= e1 && !(e1 & (e2 - e1));
    }
    MMilnor divLF(MMod rhs)
    {
#ifndef NDEBUG /* DEBUG */
        if (!rhs.divisibleLF(*this))
            throw MyException(0x4357a92cU, "Not divisible: lhs / rhs");
#endif
        return MMilnor(((data_ ^ rhs.data_) & MMILNOR_MASK_E) | (w_raw() - rhs.w_raw()));
    }

    int deg_m() const
    {
        return m_no_weight().deg();
    }

    std::string Str() const;
    std::string StrXi() const;
};
using MMod1d = std::vector<MMod>;
using MMod2d = std::vector<MMod1d>;

inline MMod mulLF(MMilnor m, MMod x)
{
    return MMod(((m.data() | x.data()) & ~MMILNOR_MASK_W) | (m.w_raw() + x.w_raw()));
}

inline bool divisibleLF(MMod lhs, MMod rhs)
{
    return lhs.divisibleLF(rhs);
}

inline MMilnor divLF(MMod lhs, MMod rhs)
{
    return lhs.divLF(rhs);
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

    Mod LF() const
    {
        Mod result;
        if (!*this)
            return result;
        result = GetLead();
        for (size_t i = 1; i < data.size(); ++i)
            if (data[0].v_raw() == data[i].v_raw() && data[0].w_raw() == data[i].w_raw())
                result.data.push_back(data[i]);
        return result;
    }

    Mod LFMay() const
    {
        Mod result;
        if (!*this)
            return result;
        result = GetLead();
        for (size_t i = 1; i < data.size(); ++i)
            if (data[0].v_raw() == data[i].v_raw() && data[0].w_may() == data[i].w_may())
                result.data.push_back(data[i]);
        return result;
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
    Mod& iadd(const Mod& rhs, Mod& tmp)
    {
        tmp.data.clear();
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data));
        return *this;
    }
    Mod& operator+=(const Mod& rhs)
    {
        Mod tmp;
        return iadd(rhs, tmp);
    }
    /* `*this += m * x` */
    Mod& iaddmul(MMilnor m, const Mod& x, Milnor& tmp_a, Mod& tmp_m1, Mod& tmp_m2);
    Mod& iaddmulMay(MMilnor m, const Mod& x, Mod& tmp);
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
Mod MulMay(MMilnor m, const Mod& x);

inline std::ostream& operator<<(std::ostream& sout, const Mod& x)
{
    return sout << x.Str();
}

inline void Reduce(Mod& x, const Mod1d& y, Mod& tmp)
{
    for (size_t i = 0; i < y.size(); ++i)
        if (std::binary_search(x.data.begin(), x.data.end(), y[i].GetLead()))
            x.iadd(y[i], tmp);
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
// inline Mod subs(const Mod& x, const Mod1d& map)
//{
//     Mod result;
//     for (const MMod& mv : x.data)
//         result += Milnor(mv.m()) * map[mv.v()];
//     return result;
// }

}  // namespace steenrod

#endif