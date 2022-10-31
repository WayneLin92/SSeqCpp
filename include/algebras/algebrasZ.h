/** \file algebras.h
 * A Component for monomials and polynomials over $F_2$.
 * All are incapsulated in the `namespace alg`.
 */

#ifndef ALGEBRASZ_H
#define ALGEBRASZ_H

#include "algebras.h"
#include <numeric>
#include <unordered_set>

/**
 * The namespace `alg` provides the basis types for monomials and polynomials.
 * and vectors of them
 */
namespace algZ {
using namespace alg;

/********************************************************
 *                      Monomials
 ********************************************************/

inline constexpr int FIL_MAX = 200;

/* Monomial for the algebra over Z(2)
 * c_ = -1 means an unknown term in filtration `fil_`
 */
class Mon
{
private:
    int fil_, c_;
    alg2::Mon m0_, m1_;

public:
    constexpr Mon() : c_(0), fil_(0) {}
    constexpr Mon(int c, const alg2::Mon& m0, const alg2::Mon& m1, int fil) : c_(c), m0_(m0), m1_(m1), fil_(fil) {}

    /* index=0 is reserved for the coefficient 2 */
    static Mon Gen(uint32_t index, uint32_t exp, int fil, bool bEven)
    {
        if (index == 0)
            return twoTo(exp);
        if (bEven)
            return Mon(0, GE(index, exp), {}, fil);
        else
            return Mon(0, {}, GE(index, exp), fil);
    }

    /* Dummy 2x^2 for groebner basis */
    static Mon Trivial(uint32_t index, int fil)
    {
#ifndef NDEBUG
        if (index == 0)
            throw MyException(0xa5b6f1a6U, "index=0 should be reserved for the coefficient 2");
#endif
        return Mon(1, {}, GE(index, 2), fil * 2 + 1);
    }

    /* An unknown term in filtration fil */
    static Mon O(int fil)
    {
        return Mon(-1, {}, {}, fil);
    }

    static Mon twoTo(int exp)
    {
        return Mon(exp, {}, {}, exp);
    }

    void SetFil(const int1d& gen_fils);

    bool operator<(const Mon& rhs) const
    {
        if (fil_ < rhs.fil_)
            return true;
        if (fil_ > rhs.fil_)
            return false;
        if (c_ < rhs.c_)
            return true;
        if (c_ > rhs.c_)
            return false;
        if (m1_ < rhs.m1_)
            return true;
        if (rhs.m1_ < m1_)
            return false;
        if (m0_ < rhs.m0_)
            return true;
        return false;
    };

    bool operator==(const Mon& rhs) const
    {
        return fil_ == rhs.fil_ && c_ == rhs.c_ && m0_ == rhs.m0_ && m1_ == rhs.m1_;
    };

    /* Return true if it is not 1 */
    explicit operator bool() const
    {
        return bool(m0_) || bool(m1_) || c_ != 0;
    }

    int c() const
    {
        return c_;
    }

    int fil() const
    {
        return fil_;
    }

    auto& m0() const
    {
        return m0_;
    }

    auto& m1() const
    {
        return m1_;
    }

    void imul2()
    {
        ++c_;
        ++fil_;
    }

    bool Is2Torsion() const
    {
        return std::any_of(m1().begin(), m1().end(), [](GE ge) { return ge.e() > 1; });
    }

    bool IsUnKnown() const
    {
        return c_ < 0;
    }

    /* return the first generator > 0 */
    uint32_t frontg() const
    {
        uint32_t result = 0;
        if (m0_)
            result = m0_.begin()->g();
        if (m1_ && result > m1_.begin()->g())
            result = m1_.begin()->g();
        return result;
    }

    /* return the largest generator */
    uint32_t backg() const
    {
        uint32_t result = 0;
        if (m0_)
            result = m0_.backg();
        if (m1_ && result < m1_.backg())
            result = m1_.backg();
        return result;
    }

    MonTrace Trace() const
    {
        return m0_.Trace() | m1_.Trace() | (c_ > 0 ? (MonTrace(1) << 63) : 0);
    }

    std::string Str() const;
};

using Mon1d = std::vector<Mon>;
using Mon2d = std::vector<Mon1d>;

inline std::ostream& operator<<(std::ostream& sout, const Mon& x)
{
    return std::cout << x.Str();
}

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template <typename FnGenDeg>
auto GetDegTpl(const Mon& mon, const FnGenDeg& _gen_deg)
{
    return _gen_deg(0) * mon.c() + alg2::GetDegTpl(mon.m0(), _gen_deg) + alg2::GetDegTpl(mon.m1(), _gen_deg);
}

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template <typename T>
auto GetDeg(const Mon& mon, const std::vector<T>& gen_degs)
{
    return GetDegTpl(mon, [&gen_degs](int i) { return gen_degs[i]; });
}

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template <typename T>
auto GetDegT(const Mon& mon, const std::vector<T>& gen_degs)
{
    return GetDegTpl(mon, [&gen_degs](int i) { return gen_degs[i].t; });
}

inline Mon mul_unsigned(const Mon& mon1, const Mon& mon2)
{
    return Mon(mon1.c() + mon2.c(), mon1.m0() * mon2.m0(), mon1.m1() * mon2.m1(), mon1.fil() + mon2.fil());
}

/* get the sign of the product */
int get_sign(const Mon& mon1, const Mon& mon2, const Mon& prod);

inline Mon div_unsigned(const Mon& mon1, const Mon& mon2)
{
#ifndef NDEBUG /* DEBUG */
    if (mon1.c() < mon2.c() || mon1.fil() < mon2.fil())
        throw MyException(0x32112d9eU, "mon1/mon2 not divisible!\n");
#endif
    return Mon(mon1.c() - mon2.c(), mon1.m0() / mon2.m0(), mon1.m1() / mon2.m1(), mon1.fil() - mon2.fil());
}
inline Mon pow(const Mon& mon, int e)
{
    return Mon(mon.c() * e, alg2::pow(mon.m0(), e), alg2::pow(mon.m1(), e), mon.fil() * e);
}

/**
 * Return if m1 divides m2.
 */
inline bool divisible(const Mon& mon1, const Mon& mon2)
{
    return mon2.c() >= mon1.c() && alg2::divisible(mon1.m0(), mon2.m0()) && alg2::divisible(mon1.m1(), mon2.m1());
}

inline bool divisible(const Mon& mon1, const Mon& mon2, MonTrace t1, MonTrace t2)
{
    return t2 >= t1 && !(t1 & (t2 - t1)) && divisible(mon1, mon2);
}

/* Warning: fil might need to be modified later */
inline Mon GCD(const Mon& mon1, const Mon& mon2)
{
    return Mon(std::min(mon1.c(), mon2.c()), alg2::GCD(mon1.m0(), mon2.m0()), alg2::GCD(mon1.m1(), mon2.m1()), -1024);
}

/* Warning: fil might need to be modified later */
inline Mon LCM(const Mon& mon1, const Mon& mon2)
{
    return Mon(std::max(mon1.c(), mon2.c()), alg2::LCM(mon1.m0(), mon2.m0()), alg2::LCM(mon1.m1(), mon2.m1()), -1024);
}

/********************************************************
 *                      Polynomials
 ********************************************************/

struct Poly;
Poly operator-(const Mon& mon);
Poly operator*(const Mon& mon1, const Mon& mon2);
void mulP(const Mon& mon, const Poly& poly, Poly& result);  ////
void mulP(const Poly& p1, const Poly& p2, Poly& result);    ////

struct Poly
{
    Mon1d data;

    Poly() {}
    Poly(const Mon& m) : data({m}) {}

    static Poly Unit()
    {
        return Mon();
    }

    static Poly Gen(uint32_t index, uint32_t exp, int fil, bool bEven)
    {
        return Mon::Gen(index, exp, fil, bEven);
    }

    static Poly twoTo(int exp)
    {
        return Mon::twoTo(exp);
    }

    const Mon& GetLead() const
    {
#ifndef NDEBUG /* DEBUG */
        if (data.empty())
            throw MyException(0x612c2027U, "Trying to get leading monomial from zero poly.");
#endif
        return data.front();
    }

    int GetDeg(const int1d& gen_degs) const
    {
        if (data.empty())
            return -1024;
        else
            return algZ::GetDeg(data.front(), gen_degs);
    }

    /* Return the filtration differences between the leading term and the big O term */
    int EffNum() const
    {
        if (*this && data.back().IsUnKnown())
            return data.back().fil() - data.front().fil();
        else
            return FIL_MAX;
    }

    /* Return the filtration of the big O term.
     * Return FIL_MAX if no such term.
     */
    int UnknownFil() const
    {
        if (*this && data.back().IsUnKnown())
            return data.back().fil();
        else
            return FIL_MAX;
    }

    bool operator==(const Poly& rhs) const
    {
        return data == rhs.data;
    }
    bool operator!=(const Poly& rhs) const
    {
        return !(data == rhs.data);
    }
    explicit operator bool() const
    {
        return !data.empty();
    }
    Poly operator-() const
    {
        Poly result, tmp;
        for (auto& m : data)
            result.iaddP(-m, tmp);
        return result;
    }
    Poly negP(Poly& tmp, const int1d& gen_2tor_degs) const;
    Poly& iaddP(const Poly& rhs, Poly& tmp);
    Poly& operator+=(const Poly& rhs)
    {
        Poly tmp;
        return iaddP(rhs, tmp);
    }
    Poly& isubP(const Poly& rhs, Poly& tmp, const int1d& gen_2tor_degs)
    {
        return iaddP(rhs.negP(tmp, gen_2tor_degs), tmp);
    }
    /* this += mon * poly */
    void iaddmulP(const Mon& mon, const Poly& poly, Poly& tmp, const int1d& gen_2tor_degs);
    /* this -= mon * poly */
    void isubmulP(const Mon& mon, const Poly& poly, Poly& tmp, const int1d& gen_2tor_degs);

    Poly operator*(const Mon& rhs) const
    {
        Poly result;
        mulP(*this, rhs, result);
        return result;
    }
    Poly operator*(const Poly& rhs) const
    {
        Poly result;
        mulP(*this, rhs, result);
        return result;
    }
    Poly& imulP(const Poly& rhs, Poly& tmp)
    {
        mulP(*this, rhs, tmp);
        std::swap(*this, tmp);
        return *this;
    }
    std::string Str() const;
};

using Poly1d = std::vector<Poly>;

inline std::ostream& operator<<(std::ostream& sout, const Poly& x)
{
    return std::cout << x.Str();
}

inline Poly operator+(const Poly& p1, const Poly& p2)
{
    Poly result(p1), tmp;
    result.iaddP(p2, tmp);
    return result;
}

void powP(const Poly& poly, uint32_t n, Poly& result, Poly& tmp);
inline Poly pow(const Poly& poly, uint32_t n)
{
    Poly result, tmp;
    powP(poly, n, result, tmp);
    return result;
}

/**
 * Replace the generators in `poly` with elements given in `map`.
 * @param poly The polynomial to be substituted.
 * @param map `map[i]` is the polynomial that substitutes the generator of id `i`.
 */
template <typename FnMap>
Poly subsTpl(const Poly& poly, const FnMap& map)
{
    Poly result, tmp_prod, tmp;
    for (const Mon& m : poly.data) {
        Poly fm = Poly::twoTo(m.c());
        for (auto p = m.m0().begin(); p != m.m0().end(); ++p) {
            powP(map(p->g()), p->e(), tmp_prod, tmp);
            fm.imulP(tmp_prod, tmp);
        }
        for (auto p = m.m1().begin(); p != m.m1().end(); ++p) {
            powP(map(p->g()), p->e(), tmp_prod, tmp);
            fm.imulP(tmp_prod, tmp);
        }
        result.iaddP(fm, tmp);
    }
    return result;
}

inline Poly subs(const Poly& poly, const std::vector<Poly>& map)
{
    return subsTpl(poly, [&map](size_t i) { return map[i]; });
}

/**
 * Replace the generators in `poly` with new id given in `map_gen_id`.
 * @param poly The polynomial to be substituted.
 * @param map_gen_id `map_gen_id[i]` is the new id that substitutes the old id `i`.
 */
inline Poly subs(const Poly& poly, const uint1d& map_gen_id, const uint1d& gen_fils, const std::unordered_set<int> gen_ids_even)
{
    return subsTpl(poly, [&map_gen_id, &gen_fils, &gen_ids_even](size_t i) { return Poly::Gen(map_gen_id[i], 1, gen_fils[i], gen_ids_even.find((int)i) != gen_ids_even.end()); });
}

int1d Poly2Indices(const Poly& poly, const Mon1d& basis);
Poly Indices2Poly(const int1d& indices, const Mon1d& basis);

/********************************************************
 *                      Modules
 ********************************************************/

struct MMod
{
    uint32_t v;
    Mon m; /* m.fil_ should include the filtration of v */

    MMod(const Mon& m_, uint32_t v_, int fil_v) : m(m_.c(), m_.m0(), m_.m1(), m_.fil() + fil_v), v(v_) {}
    MMod(int c, const alg2::Mon& m0, const alg2::Mon& m1, int fil, int v_) : m(c, m0, m1, fil), v(v_) {}

    /* An unknown term in filtration fil */
    static MMod O(int fil)
    {
        return MMod(Mon::O(fil), UINT32_MAX, 0);
    }

    bool operator<(const MMod& rhs) const
    {
        if (fil() < rhs.fil())
            return true;
        if (fil() > rhs.fil())
            return false;
        if (v > rhs.v)
            return true;
        if (v < rhs.v)
            return false;
        if (m < rhs.m)
            return true;
        return false;
    }

    bool operator==(const MMod& rhs) const
    {
        return m == rhs.m && v == rhs.v;
    }

    void imul2()
    {
        m.imul2();
    }

    int fil() const
    {
        return m.fil();
    }

    bool Is2Torsion() const
    {
        return m.Is2Torsion();
    }

    bool IsUnKnown() const
    {
        return m.IsUnKnown();
    }

    std::string Str() const;
};

using MMod1d = std::vector<MMod>;
using MMod2d = std::vector<MMod1d>;

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template <typename T>
auto GetDeg(const MMod& mon, const std::vector<T>& gen_degs, const std::vector<T>& v_degs)
{
    return GetDeg(mon.m, gen_degs) + v_degs[mon.v];
}

struct Mod;
Mod operator-(const MMod& mon);
Mod operator*(const Mon& m, const MMod& x);

struct Mod
{
    MMod1d data;
    Mod() {}
    Mod(MMod mv) : data({mv}) {}
    Mod(const Poly& poly, uint32_t v, int fil_v)
    {
        data.reserve(poly.data.size());
        for (const Mon& m : poly.data)
            data.emplace_back(m, v, fil_v);
        ut::RemoveIf(data, [](const MMod& mv) { return mv.m.fil() > FIL_MAX; });
    }

    const MMod& GetLead() const
    {
#ifndef NDEBUG
        if (data.empty())
            throw MyException(0x24b0f95cU, "Trying to GetLead() for empty Mod.");
#endif
        return data[0];
    }

    /* Return the filtration differences between the leading term and the big O term */
    int EffNum() const
    {
        if (*this && data.back().IsUnKnown())
            return data.back().fil() - data.front().fil();
        else
            return FIL_MAX;
    }

    int UnknownFil() const
    {
        if (*this && data.back().IsUnKnown())
            return data.back().fil();
        else
            return FIL_MAX;
    }

    Mod& iaddP(const Mod& rhs, Mod& tmp);

    Mod operator-() const;
    Mod& isubP(const Mod& rhs, Mod& tmp)
    {
        return iaddP(-rhs, tmp);
    }

    void iaddmulP(const Mon& m, const Mod& x, Mod& tmp, const int1d& gen_2tor_degs);
    void isubmulP(const Mon& m, const Mod& x, Mod& tmp, const int1d& gen_2tor_degs);

    void ReduceSizeByChangingSignP(Mod& tmp1, Mod& tmp2);
    void ReduceSizeByChangingSign()
    {
        Mod tmp1, tmp2;
        ReduceSizeByChangingSignP(tmp1, tmp2);
    }

    explicit operator bool() const
    {
        return !data.empty();
    }

    bool operator==(const Mod& rhs) const
    {
        return data == rhs.data;
    }

    std::string Str() const;
};

using Mod1d = std::vector<Mod>;

inline std::ostream& operator<<(std::ostream& sout, const Mod& x)
{
    return std::cout << x.Str();
}

inline Mod operator+(const Mod& x, const Mod& y)
{
    Mod result(x), tmp;
    result.iaddP(y, tmp);
    return result;
}

inline Mod operator*(const Mon& m, const Mod& x)
{
    Mod result, tmp;
    for (auto& xm : x.data)
        result.iaddP(m * xm, tmp);
    return result;
}

inline Mod operator*(const Poly& p, const Mod& x)
{
    Mod result, tmp;
    for (auto& m : p.data)
        result.iaddP(m * x, tmp);
    return result;
}

/**
 * Replace the generators in `poly` with elements given in `map`.
 * @param poly The polynomial to be substituted.
 * @param map `map[i]` is the polynomial that substitutes the generator of id `i`.
 */
template <typename FnMap>
Poly subsModTpl(const Mod& x, const FnMap& map)
{
    Poly result, tmp_prod, tmp;
    for (const MMod& m : x.data)
        result.iaddP(map(m.v) * m.m, tmp);
    return result;
}

inline Poly subsMod(const Mod& x, const std::vector<Poly>& map)
{
    return subsModTpl(x, [&map](size_t i) { return map[i]; });
}

int1d Mod2Indices(const Mod& x, const MMod1d& basis);
Mod Indices2Mod(const int1d& indices, const MMod1d& basis);

}  // namespace algZ

#endif /* ALGEBRAS_H */