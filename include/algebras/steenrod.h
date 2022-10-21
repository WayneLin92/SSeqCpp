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
using int1d = std::vector<int>;
using int2d = std::vector<int1d>;
using int3d = std::vector<int2d>;
using int4d = std::vector<int3d>;

/********************************************************
 *                  class Milnor
 ********************************************************/

inline constexpr int DEG_MAX = 383;      /* Maximum degree supported in A if `MMILNOR_E_BITS == 37` */
inline constexpr size_t XI_MAX = 8;      /* Support up to \xi_8 */
inline constexpr size_t XI_MAX_MULT = 8; /* Multiplication support up to \xi_8 */

inline constexpr size_t MMILNOR_E_BITS = 37;
inline constexpr size_t MMILNOR_W_BITS = 9;
inline constexpr uint64_t MMILNOR_ONE = uint64_t(1) << (MMILNOR_E_BITS - 1);
namespace detail {
    inline constexpr std::array<uint8_t, MMILNOR_E_BITS> MMilnorGenI()
    {
        std::array<uint8_t, MMILNOR_E_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_E_BITS)
                    result[n++] = (uint8_t)i;
        return result;
    }
    inline constexpr std::array<uint8_t, MMILNOR_E_BITS> MMilnorGenJ()
    {
        std::array<uint8_t, MMILNOR_E_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_E_BITS)
                    result[n++] = (uint8_t)j;
        return result;
    }
    inline constexpr std::array<uint8_t, MMILNOR_E_BITS> MMilnorGenJmIm1()
    {
        std::array<uint8_t, MMILNOR_E_BITS> result = {};
        size_t n = 0;
        for (size_t j = 1; j <= XI_MAX + 1; ++j)
            for (size_t i = j; i-- > 0;)
                if (n < MMILNOR_E_BITS)
                    result[n++] = (uint8_t)(j - i - 1);
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
    template <size_t d>
    constexpr std::array<uint64_t, DEG_MAX / ((1 << d) - 1)> MMilnorXiD()
    {
        constexpr uint32_t N = DEG_MAX / ((1 << d) - 1);
        std::array<uint64_t, N> result = {};
        for (uint32_t m = 0; m < N; ++m) {
            for (uint32_t n = m, i = 0; n; n >>= 1, ++i) {
                if (n & 1) {
                    uint32_t j = i + d;
                    uint32_t index = j * (j + 1) / 2 - i - 1;
                    result[m] |= MMILNOR_ONE >> index;
                }
            }
        }
        return result;
    }
}  // namespace detail
inline constexpr std::array<uint8_t, MMILNOR_E_BITS> MMILNOR_GEN_I = detail::MMilnorGenI();
inline constexpr std::array<uint8_t, MMILNOR_E_BITS> MMILNOR_GEN_J = detail::MMilnorGenJ();
inline constexpr std::array<uint8_t, MMILNOR_E_BITS> MMILNOR_GEN_JMIM1 = detail::MMilnorGenJmIm1();
inline constexpr std::array<uint64_t, MMILNOR_E_BITS> MMILNOR_GEN_WEIGHT = detail::MMilnorGenWeight();
inline constexpr std::array<int, MMILNOR_E_BITS> MMILNOR_GEN_DEG = detail::MMilnorGenDeg();
inline constexpr uint64_t UINT64_LEFT_BIT = uint64_t(1) << 63;
inline constexpr uint64_t MMILNOR_MASK_E = (uint64_t(1) << MMILNOR_E_BITS) - 1;
inline constexpr uint64_t MMILNOR_MASK_W = ((uint64_t(1) << MMILNOR_W_BITS) - 1) << MMILNOR_E_BITS;
inline constexpr uint64_t MMILNOR_NULL = ~uint64_t(0);

inline constexpr auto MMILNOR_Xi0 = detail::MMilnorXiD<1>();
inline constexpr auto MMILNOR_Xi1 = detail::MMilnorXiD<2>();
inline constexpr auto MMILNOR_Xi2 = detail::MMilnorXiD<3>();
inline constexpr auto MMILNOR_Xi3 = detail::MMilnorXiD<4>();
inline constexpr auto MMILNOR_Xi4 = detail::MMilnorXiD<5>();
inline constexpr auto MMILNOR_Xi5 = detail::MMilnorXiD<6>();
inline constexpr auto MMILNOR_Xi6 = detail::MMilnorXiD<7>();
inline constexpr auto MMILNOR_Xi7 = detail::MMilnorXiD<8>();

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

    /* This is a critial function */
    static MMilnor Xi(const uint32_t* xi)
    {
        const uint64_t w_raw = uint64_t(1 * ut::popcount(xi[0]) + 3 * ut::popcount(xi[1]) + 5 * ut::popcount(xi[2]) + 7 * ut::popcount(xi[3]) + 9 * ut::popcount(xi[4]) + 11 * ut::popcount(xi[5]) + 13 * ut::popcount(xi[6]) + 15 * ut::popcount(xi[7]))
                               << MMILNOR_E_BITS;
        const uint64_t e = MMILNOR_Xi0[xi[0]] | MMILNOR_Xi1[xi[1]] | MMILNOR_Xi2[xi[2]] | MMILNOR_Xi3[xi[3]] | MMILNOR_Xi4[xi[4]] | MMILNOR_Xi5[xi[5]] | MMILNOR_Xi6[xi[6]] | MMILNOR_Xi7[xi[7]];
        return MMilnor(w_raw | e);
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
            result[MMILNOR_GEN_JMIM1[i]] += uint32_t(1) << MMILNOR_GEN_I[i];
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
    Milnor(MMilnor m) : data({m}) {}

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
using Milnor1d = std::vector<Milnor>;

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
    Mod& iaddP(const Mod& rhs, Mod& tmp)
    {
        tmp.data.clear();
        std::set_symmetric_difference(data.cbegin(), data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(tmp.data));
        ut::copy(tmp.data, data);
        return *this;
    }
    Mod& operator+=(const Mod& rhs)
    {
        Mod tmp;
        return iaddP(rhs, tmp);
    }
    /* `*this += m * x` */
    Mod& iaddmulP(MMilnor m, const Mod& x, Milnor& tmp_a, Mod& tmp_x1, Mod& tmp_x2);
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

void mulP(MMilnor m, const Mod& x, Mod& result, Milnor& tmp);
inline Mod operator*(MMilnor m, const Mod& x)
{
    Mod result;
    Milnor tmp;
    mulP(m, x, result, tmp);
    return result;
}
void MulMayP(MMilnor m, const Mod& x, Mod& result, Milnor& tmp);
Mod MulMay(MMilnor m, const Mod& x);

inline std::ostream& operator<<(std::ostream& sout, const Mod& x)
{
    return sout << x.Str();
}

inline void Reduce(Mod& x, const Mod1d& y, Mod& tmp)
{
    for (size_t i = 0; i < y.size(); ++i)
        if (std::binary_search(x.data.begin(), x.data.end(), y[i].GetLead()))
            x.iaddP(y[i], tmp);
}

/* Compute the product in the associated graded algebra */
Mod mulLF(MMilnor m, const Mod& x);

/**
 * Replace v_i with `map[i]`.
 */
inline void subsP(const Mod& x, const Mod1d& map, Mod& result, Milnor& tmp_a, Mod& tmp_x1, Mod& tmp_x2)
{
    for (const MMod& mv : x.data)
        result.iaddmulP(mv.m(), map[mv.v()], tmp_a, tmp_x1, tmp_x2);
}

/**
 * Replace v_i with `map[i]`.
 */
inline Mod subs(const Mod& x, const Mod1d& map)
{
    Mod result, tmp_x1, tmp_x2;
    Milnor tmp_a;
    subsP(x, map, result, tmp_a, tmp_x1, tmp_x2);
    return result;
}

}  // namespace steenrod

#endif