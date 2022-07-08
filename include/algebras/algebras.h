/** \file algebras.h
 * A Component for monomials and polynomials over $F_2$.
 * All are incapsulated in the `namespace alg`.
 */

#ifndef ALGEBRAS_H
#define ALGEBRAS_H

#include "myexception.h"
#include "utility.h"
#include <array>
#include <climits>
#include <queue>

/**
 * The namespace `alg` provides the basis types for monomials and polynomials.
 * and vectors of them
 */
namespace alg {

inline constexpr int DEG_MAX = INT_MAX;

using int1d = std::vector<int>;
using int2d = std::vector<int1d>;
using int3d = std::vector<int2d>;
using int4d = std::vector<int3d>;
using pairii = std::pair<int, int>;
using pairii1d = std::vector<pairii>;

using uint1d = std::vector<uint32_t>;

/** The 3d grading for the May spectral sequence.
 *
 * v is the complement degree of the May degree.
 * $v(R_{ij})=j-i-1$.
 * Degrees are ordered by t, s, v.
 */
struct MayDeg
{
    int s, t, v;
    MayDeg() : s(0), t(0), v(0) {}
    MayDeg(int s_, int t_, int v_) : s(s_), t(t_), v(v_) {}

    bool operator<(const MayDeg& rhs) const
    {
        if (t < rhs.t)
            return true;
        if (t > rhs.t)
            return false;
        if (s < rhs.s)
            return true;
        if (s > rhs.s)
            return false;
        if (v < rhs.v)
            return true;
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

/** The 2d grading for the Adams spectral sequence.
 *
 * v is the complement degree of the May degree.
 * $v(R_{ij})=j-i-1$.
 * Degrees are ordered by t, s, v.
 */
struct AdamsDeg
{
    int s, t;
    AdamsDeg() : s(0), t(0) {}
    AdamsDeg(int s_, int t_) : s(s_), t(t_) {}
    bool operator<(const AdamsDeg& rhs) const
    {
        if (t < rhs.t)
            return true;
        if (t > rhs.t)
            return false;
        if (s < rhs.s)
            return true;
        return false;
    };
    AdamsDeg operator+(const AdamsDeg& rhs) const
    {
        return AdamsDeg{s + rhs.s, t + rhs.t};
    };
    AdamsDeg operator-(const AdamsDeg& rhs) const
    {
        return AdamsDeg{s - rhs.s, t - rhs.t};
    };
    AdamsDeg& operator+=(const AdamsDeg& rhs)
    {
        s += rhs.s;
        t += rhs.t;
        return *this;
    };
    AdamsDeg& operator-=(const AdamsDeg& rhs)
    {
        s -= rhs.s;
        t -= rhs.t;
        return *this;
    };
    bool operator==(const AdamsDeg& rhs) const
    {
        return s == rhs.s && t == rhs.t;
    };
    bool operator!=(const AdamsDeg& rhs) const
    {
        return !operator==(rhs);
    };
    AdamsDeg operator*(int rhs) const
    {
        return AdamsDeg{s * rhs, t * rhs};
    };
};
using AdamsDeg1d = std::vector<AdamsDeg>;

/** \defgroup Monomials Monomials
 * ------------------------------------------- @{
 */
/********************************************************
 *                      Monomials
 ********************************************************/

struct GE
{
    uint32_t data;
    constexpr GE() : data(0xffff << 16) {}
    constexpr explicit GE(uint32_t data_) : data(data_) {}
    constexpr GE(uint32_t g, uint32_t e) : data((~g << 16) | e) {}

    uint32_t g() const
    {
        return ~data >> 16;
    }
    uint32_t g_raw() const
    {
        return data & 0xffff0000;
    }
    uint32_t e() const
    {
        return data & 0xffff;
    }

    bool operator==(GE rhs) const
    {
        return data == rhs.data;
    };
    bool operator<(GE rhs) const
    {
        return data < rhs.data;
    };

    std::string Str() const;
};

using GE1d = std::vector<GE>;

using MonTrace = uint64_t;
using MonTrace1d = std::vector<MonTrace>;

struct Mon
{
    GE1d data;

    Mon() {}
    Mon(GE p) : data({p}) {}

    static Mon Gen(uint32_t index, uint32_t exp = 1)
    {
        return GE(index, exp);
    }

    bool operator==(const Mon& rhs) const
    {
        return data == rhs.data;
    };
    bool operator<(const Mon& rhs) const
    {
        return data < rhs.data;
    };
    explicit operator bool() const
    {
        return !data.empty();
    }
    auto& operator[](size_t i) const
    {
        return data[i];
    }
    auto begin() const
    {
        return data.begin();
    }
    auto end() const
    {
        return data.end();
    }
    auto size() const
    {
        return data.size();
    }
    void push_back(GE p)
    {
        data.push_back(p);
    }

    MonTrace Trace() const;
    std::string Str() const;
};

using Mon1d = std::vector<Mon>;
using Mon2d = std::vector<Mon1d>;

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template <typename FnGenDeg>
inline auto GetDegTpl(const Mon& mon, const FnGenDeg& _gen_deg)
{
    using TypeReturn = decltype(_gen_deg(0));
    auto result = TypeReturn();
    for (auto p = mon.begin(); p != mon.end(); ++p)
        result += _gen_deg(p->g()) * p->e();
    return result;
}

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template <typename T>
inline auto GetDeg(const Mon& mon, const std::vector<T>& gen_degs)
{
    return GetDegTpl(mon, [&gen_degs](int i) { return gen_degs[i]; });
}

/**
 * Obtain the t degree of a monomial given the degrees of generators.
 */
template <typename T>
inline int GetDegT(const Mon& mon, const std::vector<T>& gen_degs)
{
    return GetDegTpl(mon, [&gen_degs](int i) { return gen_degs[i].t; });
}

void mulP(const Mon& mon1, const Mon& mon2, Mon& result);
inline Mon mul(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    mulP(mon1, mon2, result);
    return result;
}

void divP(const Mon& mon1, const Mon& mon2, Mon& result);
inline Mon div(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    divP(mon1, mon2, result);
    return result;
}

void powP(const Mon& mon, int e, Mon& result);
inline Mon pow(const Mon& mon, int e)
{
    Mon result;
    powP(mon, e, result);
    return result;
}

/**
 * Return if m1 divides m2.
 */
bool divisible(const Mon& mon1, const Mon& mon2);

inline bool divisible(const Mon& mon1, const Mon& mon2, MonTrace t1, MonTrace t2)
{
    return t2 >= t1 && !(t1 & (t2 - t1)) && divisible(mon1, mon2);
}

/**
 *  Return the largest integer e where m1 = m2^e * r.
 */
int log(const Mon& mon1, const Mon& mon2);


void GcdP(const Mon& mon1, const Mon& mon2, Mon& result);
inline Mon GCD(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    GcdP(mon1, mon2, result);
    return result;
}

void LcmP(const Mon& mon1, const Mon& mon2, Mon& result);
inline Mon LCM(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    LcmP(mon1, mon2, result);
    return result;
}

/** @} ---------------------------------------- */
/** \defgroup Polynomials Polynomials
 * ------------------------------------------- @{
 */

/********************************************************
 *                      Polynomials
 ********************************************************/

struct Poly;
void mulP(const Poly& p1, const Poly& p2, Poly& result);
void mulP(const Poly& poly, const Mon& mon, Poly& result);

struct Poly  // TODO: change to Poly
{
    Mon1d data;

    Poly() {}
    Poly(Mon m) : data({std::move(m)}) {}

    static Poly Unit()
    {
        return Mon();
    }

    static Poly Gen(uint32_t index, uint32_t exp = 1)
    {
        return Mon::Gen(index, exp);
    }

    const Mon& GetLead() const
    {
#ifndef NDEBUG /* DEBUG */
        if (data.empty())
            throw MyException(0x21eab89e, "Trying to get leading monomial from zero poly.");
#endif
        return data.front();
    }

    int GetDeg(const int1d& gen_degs) const
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
            return alg::GetDeg(data.front(), gen_degs);
    }

    AdamsDeg GetMayDeg(const AdamsDeg1d& gen_degs) const
    {
        if (data.empty())
            return AdamsDeg{-10000, -10000};
        else
            return alg::GetDeg(data.front(), gen_degs);
    }

    int GetMayDegT(const MayDeg1d& gen_degs) const
    {
        if (data.empty())
            return -10000;
        else
            return alg::GetDegT(data.front(), gen_degs);
    }

    int GetAdamsDegT(const AdamsDeg1d& gen_degs) const
    {
        if (data.empty())
            return -10000;
        else
            return alg::GetDegT(data.front(), gen_degs);
    }

    void frobP(Poly& result) const
    {
        result.data.clear();
        for (const Mon& m : data)
            result.data.push_back(pow(m, 2));
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
    Poly operator+(const Poly& rhs) const
    {
        Poly result;
        std::set_symmetric_difference(data.cbegin(), data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data));
        return result;
    }
    Poly& iaddP(const Poly& rhs, Poly& tmp)
    {
        tmp.data.clear();
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data));
        return *this;
    }
    Poly& operator+=(const Poly& rhs)
    {
        Poly tmp;
        return iaddP(rhs, tmp);
    }
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

void powP(const Poly& poly, uint32_t n, Poly& result, Poly& tmp);
inline Poly pow(const Poly& poly, uint32_t n) {
    Poly result, tmp;
    powP(poly, n, result, tmp);
    return result;
}

/** @} ---------------------------------------- */

/**
 * Compute the differential of a monomial.
 * `diffs` is the array $(dg_i)$.
 */
Poly GetDiff(const Mon& mon, const Poly1d& diffs);

/**
 * Compute the differential of a polynomial.
 * `diffs` is the array $(dg_i)$.
 */
Poly GetDiff(const Poly& poly, const Poly1d& diffs);

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
        Poly fm = Poly::Unit();
        for (auto p = m.begin(); p != m.end(); ++p) {
            powP(map(p->g()), p->e(), tmp_prod, tmp);
            fm.imulP(tmp_prod, tmp);
        }
        result += fm;
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
Poly subs(const Poly& poly, const uint1d& map_gen_id)
{
    return subsTpl(poly, [&map_gen_id](size_t i) { return Poly::Gen(map_gen_id[i]); });
}

uint1d Poly2Indices(const Mon1d& poly, const Mon1d& basis);
Mon1d Indices2Poly(const uint1d& indices, const Mon1d& basis);

}  // namespace alg

#endif /* ALGEBRAS_H */