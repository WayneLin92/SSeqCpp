/** \file algebras.h
 * A Component for monomials and polynomials over $F_2$.
 * All are incapsulated in the `namespace alg`.
 */

#ifndef ALGEBRAS_H
#define ALGEBRAS_H

#include "myexception.h"
#include "utility.h"
#include <array>
#include <queue>

/**
 * The namespace `alg` provides the basis types for monomials and polynomials.
 * and vectors of them
 */
namespace alg {

inline constexpr int DEG_MAX = INT_MAX;

using array = std::vector<int>;
using array2d = std::vector<array>;
using array3d = std::vector<array2d>;
using array4d = std::vector<array3d>;
using pairint = std::pair<int, int>;
using pairint1d = std::vector<pairint>;

/** The 3d grading for the May spectral sequence.
 *
 * v is the complement degree of the May degree.
 * $v(R_{ij})=j-i-1$.
 * Degrees are ordered by t, s, v.
 */
struct MayDeg
{
    int s, t, v;
    bool operator<(const MayDeg& rhs) const
    {
        if (t < rhs.t)
            return true;
        else if (t == rhs.t) {
            if (s < rhs.s)
                return true;
            else if (s == rhs.s)
                if (v < rhs.v)
                    return true;
        }
        return false;
    };
    MayDeg operator+(const MayDeg& rhs) const
    {
        return MayDeg{s + rhs.s, t + rhs.t, v + rhs.v};
    };
    MayDeg operator-(const MayDeg& rhs) const
    {
        return MayDeg{s - rhs.s, t - rhs.t, v - rhs.v};
    };
    MayDeg& operator+=(const MayDeg& rhs)
    {
        s += rhs.s;
        t += rhs.t;
        v += rhs.v;
        return *this;
    };
    MayDeg& operator-=(const MayDeg& rhs)
    {
        s -= rhs.s;
        t -= rhs.t;
        v -= rhs.v;
        return *this;
    };
    bool operator==(const MayDeg& rhs) const
    {
        return v == rhs.v && s == rhs.s && t == rhs.t;
    };
    bool operator!=(const MayDeg& rhs) const
    {
        return t != rhs.t || s != rhs.s || v != rhs.v;
    };
    MayDeg operator*(int rhs) const
    {
        return MayDeg{s * rhs, t * rhs, v * rhs};
    };
};
using MayDeg1d = std::vector<MayDeg>;

/********************************************************
 *                      Monomials
 ********************************************************/
/** \defgroup Monomials Monomials
 * ------------------------------------------- @{
 */

/**
 * A factor of a monomial which represents a power of a generator.
 */
struct GenPow
{
    int gen; /**< ID for the generator. */
    int exp; /**< Exponent. */
    GenPow() = default;
    GenPow(int g, int e) : gen(g), exp(e) {} /**< The constructor. */
    /**
     * The operator< helps defining a monomial ordering.
     */
    bool operator<(const GenPow& rhs) const
    {
        return gen > rhs.gen || (gen == rhs.gen && exp < rhs.exp);
    }
    /** The operator == */
    bool operator==(const GenPow& rhs) const
    {
        return gen == rhs.gen && exp == rhs.exp;
    }
    static std::vector<GenPow> Mon(int g, int e = 1)
    {
        return std::vector<GenPow>{{g, e}};
    }
};
using Mon = std::vector<GenPow>;
using Mon1d = std::vector<Mon>;
using Mon2d = std::vector<Mon1d>;
using MonInd = Mon::const_iterator;
using MonTrace = uint64_t;
using MonTrace1d = std::vector<MonTrace>;

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template<typename FnGenDeg>
inline int TplGetDeg(const Mon& mon, FnGenDeg _gen_deg)
{
    int result = 0;
    for (MonInd p = mon.begin(); p != mon.end(); ++p)
        result += _gen_deg(p->gen) * p->exp;
    return result;
}
/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
inline int GetDeg(const Mon& mon, const array& gen_degs)
{
    return TplGetDeg(mon, [&gen_degs](int i) { return gen_degs[i]; });
}
/**
 * Obtain the t degree of a monomial given the degrees of generators.
 */
inline int GetDegT(const Mon& mon, const MayDeg1d& gen_degs)
{  
    return TplGetDeg(mon, [&gen_degs](int i) { return gen_degs[i].t; });
}
/**
 * Obtain the May degree of a monomial given the May degrees of generators.
 */
inline MayDeg GetMayDeg(const Mon& mon, const MayDeg1d& gen_maydegs)
{
    MayDeg result{0, 0, 0};
    for (MonInd p = mon.begin(); p != mon.end(); ++p)
        result += gen_maydegs[p->gen] * p->exp;
    return result;
}

Mon mul(const Mon& mon1, const Mon& mon2);

/**
 * The funtions returns the quotient of the monomials.
 * It requires that mon2 divides mon1.
 */
Mon div(const Mon& mon1, const Mon& mon2);

Mon pow(const Mon& m, int e);

/**
 * Return if m1 divides m2.
 */
bool divisible(const Mon& m1, const Mon& m2);

/* A function to accelerate the computation of `divisible(m1, m2)` */
inline MonTrace Trace(const Mon& lead)
{
    MonTrace result = 0;
    for (size_t i = 0; i < lead.size(); ++i) {
        const int bits_exp1 = 56;
        const int bits_exp2 = 64 - bits_exp1;
        result |= (MonTrace(1) << (lead[i].gen % bits_exp1));
        if (lead[i].exp >= 2)
            result |= (MonTrace(1) << ((lead[i].gen % bits_exp2) + bits_exp1));
    }
    return result;
}

inline bool divisible(const Mon& m1, const Mon& m2, MonTrace t1, MonTrace t2)
{
    return t2 >= t1 && !(t1 & (t2 - t1)) && divisible(m1, m2);
}

/**
 *  Return the largest integer e where m1 = m2^e * r.
 */
int log(const Mon& m1, const Mon& m2);

/**
 * Greatest common divisor.
 */
Mon GCD(const Mon& m1, const Mon& m2);
/**
 * Least common multiple.
 */
Mon LCM(const Mon& m1, const Mon& m2);

namespace detail {
    /*
     * This version of `Mon` is designed to avoid the allocations of memory.
     */
    class MonOnStack
    {
        static constexpr size_t BUFFER_SIZE = 30;
    public:
        using const_iterator = std::array<GenPow, BUFFER_SIZE>::const_iterator;

    private:
        std::array<GenPow, BUFFER_SIZE> data_;
        size_t size_;

    public:
        MonOnStack() : size_(0){};
        explicit operator Mon() const
        {
            Mon result;
            result.reserve(size_);
            for (auto p = begin(); p < end(); ++p)
                result.push_back(*p);
            return result;
        }
        size_t size() const
        {
            return size_;
        }
        const_iterator begin() const
        {
            return data_.begin();
        }
        const_iterator end() const
        {
            return data_.begin() + size_;
        }
        GenPow back() const
        {
            return data_[size_ - 1];
        }
        void clear()
        {
            size_ = 0;
        }
        void push_back(GenPow gp)
        {
            data_[size_++] = gp;
        }
        void push_back(MonInd begin, MonInd end)
        {
            for (auto p = begin; p < end; ++p)
                data_[size_++] = *p;
        }
        void push_back(const_iterator begin, const_iterator end)
        {
            for (auto p = begin; p < end; ++p)
                data_[size_++] = *p;
        }
        void emplace_back(int gen, int exp)
        {
            data_[size_++] = {gen, exp};
        }
    };

    void mul(const Mon& mon1, const MonOnStack& mon2, Mon& result);
}  // namespace detail

/** @} ---------------------------------------- */

/********************************************************
 *                      Polynomials
 ********************************************************/
/** \defgroup Polynomials Polynomials
 * ------------------------------------------- @{
 */

/**
 * Lexicographical monomial ordering
 */
struct CmpLex
{
    static constexpr std::string_view name = "Lex";
    template <typename TypeIter1, typename TypeIter2>
    static bool cmp_ranges(TypeIter1 m1begin, TypeIter1 m1end, TypeIter2 m2begin, TypeIter2 m2end)
    {
        return std::lexicographical_compare(m2begin, m2end, m1begin, m1end); /* m1 > m2 */
    }
    template <typename Type1, typename Type2>
    static bool cmp(const Type1& m1, const Type2& m2)
    {
        return cmp_ranges(m1.begin(), m1.end(), m2.begin(), m2.end());
    }
};

/**
 * Reversed lexicographical monomial ordering
 */
struct CmpRevlex
{
    static constexpr std::string_view name = "Revlex";
    template <typename TypeIter1, typename TypeIter2>
    static bool cmp_ranges(TypeIter1 m1begin, TypeIter1 m1end, TypeIter2 m2begin, TypeIter2 m2end)
    {
        return std::lexicographical_compare(m1begin, m1end, m2begin, m2end); /* m1 < m2 */
    }
    template <typename Type1, typename Type2>
    static bool cmp(const Type1& m1, const Type2& m2)
    {
        return cmp_ranges(m1.begin(), m1.end(), m2.begin(), m2.end());
    }
};

/**
 * Polynomial with a monomial ordering as the template argument
 */
template <typename FnCmp>
struct Polynomial
{
    Mon1d data;

    static Polynomial<FnCmp> Unit()
    {
        return Polynomial<FnCmp>{{{}}};
    }

    static Polynomial<FnCmp> Gen(int index)
    {
        return Polynomial<FnCmp>{{{{index, 1}}}};
    }

    static Polynomial<FnCmp> GenExp(int index, int exp)
    {
        if (exp > 0)
            return Polynomial<FnCmp>{{{{index, exp}}}};
        else if (exp == 0)
            return Polynomial<FnCmp>{{
                {},
            }};
        else
            throw MyException(0x20eb6831U, "Negative exponent.");
    }

    static Polynomial<FnCmp> Mon_(Mon mon)
    {
        return Polynomial<FnCmp>{{std::move(mon)}};
    }

    static Polynomial<FnCmp> Sort(Mon1d data)
    {
        Polynomial<FnCmp> result = {std::move(data)};
        std::sort(result.data.begin(), result.data.end(), FnCmp::template cmp<Mon, Mon>);
        return result;
    }

    const Mon& GetLead() const
    {
#ifndef NDEBUG /* DEBUG */
        if (data.empty())
            throw MyException(0x21eab89e, "Trying to get leading monomial from zero poly.");
#endif
        return data.front();
    }

    int GetDeg(const array& gen_degs) const
    {
        if (data.empty())
            return -10000;
        else
            return alg::GetDeg(data.front(), gen_degs);
    }

    MayDeg GetMayDeg(const MayDeg1d& gen_degs) const
    {
        if (data.empty())
            return MayDeg{-10000, -10000, -10000};
        else
            return alg::GetMayDeg(data.front(), gen_degs);
    }

    int GetMayDegT(const MayDeg1d& gen_degs) const
    {
        if (data.empty())
            return -10000;
        else
            return alg::GetDegT(data.front(), gen_degs);
    }

    Polynomial<FnCmp> Square() const
    {
        Polynomial<FnCmp> result;
        for (const Mon& m : data)
            result.data.push_back(pow(m, 2));
        return result;
    }

    bool operator==(const Polynomial<FnCmp>& rhs) const
    {
        return data == rhs.data;
    }
    bool operator!=(const Polynomial<FnCmp>& rhs) const
    {
        return !(data == rhs.data);
    }
    explicit operator bool() const
    {
        return !data.empty();
    }
    Polynomial<FnCmp> operator+(const Polynomial<FnCmp>& rhs) const
    {
        Polynomial<FnCmp> result;
        std::set_symmetric_difference(data.cbegin(), data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data), FnCmp::template cmp<Mon, Mon>);
        return result;
    }
    Polynomial<FnCmp>& operator+=(const Polynomial<FnCmp>& rhs)
    {
        Polynomial<FnCmp> tmp;
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data), FnCmp::template cmp<Mon, Mon>);
        return *this;
    }
    Polynomial<FnCmp> operator*(const Mon& rhs) const
    {
        Polynomial<FnCmp> result;
        result.data.reserve(data.size());
        for (const Mon& m : data)
            result.data.push_back(mul(m, rhs));
        return result;
    }
    Polynomial<FnCmp> operator*(const detail::MonOnStack& rhs) const
    {
        Polynomial<FnCmp> result;
        result.data.resize(data.size());
        for (size_t i = 0; i < data.size(); ++i)
            mul(data[i], rhs, result.data[i]);
        return result;
    }
    Polynomial<FnCmp> operator*(const Polynomial<FnCmp>& rhs) const
    {
        Polynomial<FnCmp> result;
        for (int k = 0; k <= (int)data.size() + (int)rhs.data.size() - 2; ++k) {
            int i_min = std::max(k - (int)rhs.data.size() + 1, 0);
            int i_max = std::min((int)data.size() - 1, k);
            for (int i = i_min; i <= i_max; ++i)
                result.data.push_back(mul(data[i], rhs.data[k - i]));
        }
        std::sort(result.data.begin(), result.data.end(), FnCmp::template cmp<Mon, Mon>);
        for (size_t i = 1; i < result.data.size(); ++i) {
            if (result.data[i] == result.data[i - 1]) {
                result.data[i].clear();
                result.data[i - 1].clear();
                ++i;
            }
        }
        ut::RemoveEmptyElements(result.data);
        return result;
    }
};
using PolyLex = Polynomial<CmpLex>;
using PolyLex1d = std::vector<PolyLex>;
using PolyRevlex = Polynomial<CmpRevlex>;
using PolyRevlex1d = std::vector<PolyRevlex>;
template<typename FnCmp>
using Polynomial1d = std::vector<Polynomial<FnCmp>>;

/**
 * A fast algorithm that computes
 * `poly ** n`
 */
template <typename FnCmp>
Polynomial<FnCmp> pow(const Polynomial<FnCmp>& poly, int n)
{
    using Poly = Polynomial<FnCmp>;
    Poly result = Poly::Unit();
    if (n == 0)
        return result;
    Poly power = poly;
    while (n) {
        if (n % 2 != 0)
            result = result * power;
        n >>= 1;
        if (n)
            power = power.Square();
    }
    return result;
}

/** @} ---------------------------------------- */

/**
 * Hash a monomial.
 */
inline size_t hash(const Mon& mon)
{
    std::size_t seed = 0;
    for (auto& ge : mon) {
        seed ^= (size_t)ge.gen + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= (size_t)ge.exp + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

/**
 * Convert a polynomial to an array of hashes of the monomials.
 */
inline array hash1d(const Mon1d& poly)
{
    array result;
    for (auto& mon : poly)
        result.push_back((int)hash(mon));
    std::sort(result.begin(), result.end());
    return result;
}

/**
 * Compute the differential of a monomial.
 * `diffs` is the array $(dg_i)$.
 */
template <typename FnCmp>
Polynomial<FnCmp> GetDiff(const Mon& mon, const std::vector<Polynomial<FnCmp>>& diffs)
{
    using Poly = Polynomial<FnCmp>;
    Poly result;
    for (MonInd k = mon.begin(); k != mon.end(); ++k) {
        if (k->exp % 2)
            result += diffs[k->gen] * div(mon, GenPow::Mon(k->gen));
    }
    return result;
}
/**
 * Compute the differential of a polynomial.
 * `diffs` is the array $(dg_i)$.
 */
template <typename FnCmp>
Polynomial<FnCmp> GetDiff(const Polynomial<FnCmp>& poly, const std::vector<Polynomial<FnCmp>>& diffs)
{
    using Poly = Polynomial<FnCmp>;
    Poly result;
    for (const Mon& mon : poly.data) {
        for (MonInd k = mon.begin(); k != mon.end(); ++k) {
            if (k->exp % 2)
                result += diffs[k->gen] * div(mon, GenPow::Mon(k->gen));
        }
    }
    return result;
}

/**
 * Replace the generators in `poly` with elements given in `map`.
 * @param poly The polynomial to be substituted.
 * @param map `map[i]` is the polynomial that substitutes the generator of id `i`.
 */
template <typename FnCmp, typename FnMap>
Polynomial<FnCmp> subs(const Mon1d& data, FnMap map)
{
    using Poly = Polynomial<FnCmp>;
    Poly result;
    for (const Mon& m : data) {
        Poly fm = Poly::Unit();
        for (MonInd p = m.begin(); p != m.end(); ++p)
            fm = fm * pow(map(p->gen), p->exp);
        result += fm;
    }
    return result;
}
template <typename FnCmp>
Polynomial<FnCmp> subs(const Mon1d& data, const std::vector<Polynomial<FnCmp>>& map)
{
    return subs<FnCmp>(data, [&map](int i) { return map[i]; });
}
/**
 * Replace the generators in `poly` with new id given in `map_gen_id`.
 * @param poly The polynomial to be substituted.
 * @param map_gen_id `map_gen_id[i]` is the new id that substitutes the old id `i`.
 */
template <typename FnCmp>
Polynomial<FnCmp> subs(const Mon1d& data, const array& map_gen_id)
{
    using Poly = Polynomial<FnCmp>;
    Poly result;
    for (const Mon& m : data) {
        Mon m1;
        for (GenPow ge : m)
            m1.push_back({map_gen_id[ge.gen], ge.exp});
        std::sort(m1.begin(), m1.end(), [](const GenPow& lhs, const GenPow& rhs) { return lhs.gen < rhs.gen; });
        result.data.push_back(m1);
    }
    std::sort(result.data.begin(), result.data.end(), FnCmp::template cmp<Mon, Mon>);
    return result;
}

array Poly2Indices(const Mon1d& poly, const Mon1d& basis);
Mon1d Indices2Poly(const array& indices, const Mon1d& basis);

}  // namespace alg

#endif /* ALGEBRAS_H */