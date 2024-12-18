/** \file algebras.h
 * A Component for monomials and polynomials over $F_2$.
 * All are incapsulated in the `namespace alg`.
 */

#ifndef ALGEBRAS_H
#define ALGEBRAS_H

#include "myexception.h"
#include "utility.h"
//#include "benchmark.h"////
#include <algorithm>
#include <array>
#include <climits>
#include <iostream>

/**
 * The namespace `alg` provides some basis types for differential alg namespaces.
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
    bool operator==(const MayDeg& rhs) const = default;
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
 * Degrees are ordered by t, s.
 */
struct AdamsDeg
{
    int s, t;
    constexpr AdamsDeg() : s(0), t(0) {}
    constexpr AdamsDeg(int s_, int t_) : s(s_), t(t_) {}
    int stem() const
    {
        return t - s;
    }

    constexpr static AdamsDeg Null()
    {
        return AdamsDeg(-1024, -1024);
    }
    bool IsNull() const
    {
        return s == -1024 && t == -1024;
    }
    constexpr bool operator<(AdamsDeg rhs) const
    {
        if (t < rhs.t)
            return true;
        if (t > rhs.t)
            return false;
        if (s < rhs.s)
            return true;
        return false;
    };
    constexpr AdamsDeg operator+(const AdamsDeg& rhs) const
    {
        return AdamsDeg{s + rhs.s, t + rhs.t};
    };
    constexpr AdamsDeg operator-(const AdamsDeg& rhs) const
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
    bool operator==(const AdamsDeg& rhs) const = default;
    constexpr AdamsDeg operator*(int rhs) const
    {
        return AdamsDeg{s * rhs, t * rhs};
    };
    constexpr AdamsDeg operator%(const AdamsDeg& rhs) const
    {
        int q = rhs.s == 0 ? t / rhs.t : (rhs.t == 0 ? s / rhs.s : std::min(s / rhs.s, t / rhs.t));
        return *this - rhs * q;
    };
    std::string Str() const
    {
        return '(' + std::to_string(stem()) + ',' + std::to_string(s) + ')';
    }
};
using AdamsDeg1d = std::vector<AdamsDeg>;
using AdamsDeg2d = std::vector<AdamsDeg1d>;

inline std::ostream& operator<<(std::ostream& sout, const AdamsDeg& deg)
{
    return sout << deg.Str();
}

struct GE
{
    uint32_t data;
    constexpr GE() : data(0xffffu << 16) {} /* default value x_0^0 */
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
    uint32_t e_masked() const
    {
        return data & 0x7fff;
    }

    bool operator==(GE rhs) const
    {
        return data == rhs.data;
    };
    bool operator!=(GE rhs) const
    {
        return data != rhs.data;
    };
    bool operator<(GE rhs) const
    {
        return data < rhs.data;
    };

    std::string Str() const;
};

using GE1d = std::vector<GE>;

/* For fast divisibility check */
using MonTrace = uint64_t;
using MonTrace1d = std::vector<MonTrace>;

}  // namespace alg

/**
 * The namespace `alg2` provides the basis types for monomials and polynomials.
 * and vectors of them
 */
namespace alg2 {
using namespace alg;

/** \defgroup Monomials Monomials
 * ------------------------------------------- @{
 */
/********************************************************
 *                      Monomials
 ********************************************************/
constexpr size_t MONSIZE = 12;

class Mon
{
public:
    using It = std::array<GE, MONSIZE>::const_iterator;

private:
    uint32_t size_;
    std::array<GE, MONSIZE> data_;

public:
    constexpr Mon() : size_(0) {}
    constexpr Mon(GE p) : size_(1)
    {
        data_[0] = p;
    }

    /* This is made for `SortMod2()` */
    constexpr static Mon Null()
    {
        Mon result;
        result.size_ = 0xffffffff;
        return result;
    }

    static Mon Gen(uint32_t index, uint32_t exp = 1)
    {
        return GE(index, exp);
    }

    bool operator==(const Mon& rhs) const
    {
        if (size_ != rhs.size_)
            return false;
        for (size_t i = 0; i < size_; ++i)
            if (data_[i] != rhs.data_[i])
                return false;
        return true;
    };

    explicit operator bool() const
    {
        return size_;
    }

    bool IsGen() const
    {
        return size_ == 1 && data_[0].e() == 1;
    }

    /* This is made for `SortMod2()` */
    bool IsNull() const
    {
        return size_ == 0xffffffff;
    }

    auto& front() const
    {
#ifndef NDEBUG
        if (size_ == 0)
            throw ErrorIdMsg(0xfebc8802U, "Calling front() on empty");
#endif
        return data_[0];
    }

    auto& back() const
    {
#ifndef NDEBUG
        if (size_ == 0)
            throw ErrorIdMsg(0xfebc8802U, "Calling back() on empty");
#endif
        return data_[size_t(size_ - 1)];
    }
    uint32_t backg() const
    {
        return back().g();
    }

    auto begin() const
    {
        return data_.begin();
    }

    auto end() const
    {
        return data_.begin() + size_;
    }

    constexpr auto size() const
    {
        return size_;
    }
    void push_back(GE p)
    {
        if (size_ >= MONSIZE) {
            throw RunTimeError("alg2::Mon overflow");
            /*printf("alg2::Mon overflow\n");
            std::exit(79018145);*/
        }
        data_[size_++] = p;
    }
    auto& operator[](size_t i) const
    {
        return data_[i];
    }
    void insert(It it_begin, It it_end)
    {
        for (It it = it_begin; it != it_end; ++it)
            push_back(*it);
    }
    bool operator<(const Mon& rhs) const
    {
        return std::lexicographical_compare(begin(), end(), rhs.begin(), rhs.end());
    };
    void clear()
    {
        size_ = 0;
    }

    MonTrace Trace() const;
    std::string Str() const;
};

using Mon1d = std::vector<Mon>;
using Mon2d = std::vector<Mon1d>;

inline std::ostream& operator<<(std::ostream& sout, const Mon& x)
{
    return sout << x.Str();
}

/**
 * Obtain the degree of a monomial given the degrees of generators.
 */
template <typename FnGenDeg>
inline auto GetDegTpl(const Mon& mon, FnGenDeg _gen_deg)
{
    using TypeReturn = decltype(_gen_deg(0));
    auto result = TypeReturn();
    for (auto p = mon.begin(); p != mon.end(); ++p)
        result += _gen_deg(p->g()) * p->e_masked();
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

Mon operator*(const Mon& mon1, const Mon& mon2);
Mon operator/(const Mon& mon1, const Mon& mon2);
Mon pow(const Mon& mon, int e);

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

Mon GCD(const Mon& mon1, const Mon& mon2);
Mon LCM(const Mon& mon1, const Mon& mon2);

/********************************************************
 *                      Polynomials
 ********************************************************/

struct Poly;
void mulP(const Poly& poly, const Mon& mon, Poly& result);
void mulP(const Poly& p1, const Poly& p2, Poly& result);

struct Poly
{
    Mon1d data;

    Poly() {}
    Poly(const Mon& m) : data({m}) {}

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
            throw ErrorIdMsg(0x21eab89e, "Trying to get leading monomial from zero poly.");
#endif
        return data.front();
    }

    int GetDeg(const int1d& gen_degs) const
    {
        if (data.empty())
            return -10000;
        else
            return alg2::GetDeg(data.front(), gen_degs);
    }

    MayDeg GetMayDeg(const MayDeg1d& gen_degs) const
    {
        if (data.empty())
            return MayDeg{-10000, -10000, -10000};
        else
            return alg2::GetDeg(data.front(), gen_degs);
    }

    AdamsDeg GetAdamsDeg(const AdamsDeg1d& gen_degs) const
    {
        if (data.empty())
            return AdamsDeg::Null();
        else
            return alg2::GetDeg(data.front(), gen_degs);
    }

    int GetMayDegT(const MayDeg1d& gen_degs) const
    {
        if (data.empty())
            return -10000;
        else
            return alg2::GetDegT(data.front(), gen_degs);
    }

    int GetAdamsDegT(const AdamsDeg1d& gen_degs) const
    {
        if (data.empty())
            return -10000;
        else
            return alg2::GetDegT(data.front(), gen_degs);
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
        std::set_symmetric_difference(data.cbegin(), data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(tmp.data));
        ut::copy(tmp.data, data);
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
using Poly2d = std::vector<Poly1d>;

inline Poly operator*(const Mon& m, const Poly& x)
{
    return x * m;
}

inline std::ostream& operator<<(std::ostream& sout, const Poly& x)
{
    return sout << x.Str();
}

void powP(const Poly& poly, uint32_t n, Poly& result, Poly& tmp);
inline Poly pow(const Poly& poly, uint32_t n)
{
    Poly result, tmp;
    powP(poly, n, result, tmp);
    return result;
}

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
            powP(map(p->g()), p->e_masked(), tmp_prod, tmp);
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
inline Poly subs(const Poly& poly, const uint1d& map_gen_id)
{
    return subsTpl(poly, [&map_gen_id](size_t i) { return Poly::Gen(map_gen_id[i]); });
}

/********************************************************
 *                      Modules
 ********************************************************/

struct MMod
{
    uint32_t v;
    Mon m;

    constexpr MMod(const Mon& m_, uint32_t v_) : m(m_), v(v_) {}

    /* This is made for `SortMod2()` */
    constexpr static MMod Null()
    {
        MMod result(Mon::Null(), 0);
        return result;
    }

    /* This is made for `SortMod2()` */
    bool IsNull() const
    {
        return m.IsNull();
    }

    bool operator<(const MMod& rhs) const
    {
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

    bool IsGen() const
    {
        return !m;
    }

    std::string Str() const;
};

using MMod1d = std::vector<MMod>;
using MMod2d = std::vector<MMod1d>;

inline MMod operator*(const Mon& m, const MMod& x)
{
    return MMod(m * x.m, x.v);
}

struct Mod
{
    MMod1d data;
    Mod() {}
    Mod(MMod mv) : data({mv}) {}
    Mod(const Poly& poly, uint32_t v)
    {
        data.reserve(poly.data.size());
        for (const Mon& m : poly.data)
            data.emplace_back(m, v);
    }

    MMod GetLead() const
    {
#ifndef NDEBUG
        if (data.empty())
            throw ErrorIdMsg(0x488793f0U, "Trying to GetLead() for empty Mod.");
#endif
        return data[0];
    }

    Mod& iaddP(const Mod& rhs, Mod& tmp)
    {
        tmp.data.clear();
        std::set_symmetric_difference(data.cbegin(), data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(tmp.data));
        ut::copy(tmp.data, data);
        return *this;
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
    return sout << x.Str();
}

inline Mod operator+(const Mod& x, const Mod& y)
{
    Mod result;
    std::set_symmetric_difference(x.data.cbegin(), x.data.cend(), y.data.cbegin(), y.data.cend(), std::back_inserter(result.data));
    return result;
}

inline Mod operator*(const Mon& m, const Mod& x)
{
    Mod result;
    for (auto& xm : x.data)
        result.data.emplace_back(m * xm.m, xm.v);
    return result;
}

inline Mod operator*(const Poly& p, const Mod& x)
{
    Mod result, tmp;
    for (auto& m : p.data)
        result.iaddP(m * x, tmp);
    return result;
}

void mulP(const Mon& mon, const Mod& mod, Mod& result);
void mulP(const Poly& poly, const Mod& mod, Mod& result);

/**
 * Replace the generators in `poly` with elements given in `map`.
 * @param poly The polynomial to be substituted.
 * @param map `map[i]` is the polynomial that substitutes the generator of id `i`.
 */
template <typename FnMap>
auto subsModTpl(const Mod& x, const FnMap& map)
{
    using T = decltype(map(0));
    T result{}, tmp{};
    for (const MMod& m : x.data)
        result.iaddP(m.m * map(m.v), tmp);
    return result;
}

inline Poly subs(const Mod& x, const std::vector<Poly>& map)
{
    return subsModTpl(x, [&map](size_t i) { return map[i]; });
}

inline Mod subs(const Mod& x, const std::vector<Mod>& map)
{
    return subsModTpl(x, [&map](size_t i) { return map[i]; });
}

int1d Poly2Indices(const Poly& poly, const Mon1d& basis);
Poly Indices2Poly(const int1d& indices, const Mon1d& basis);
int1d Mod2Indices(const Mod& x, const MMod1d& basis);
Mod Indices2Mod(const int1d& indices, const MMod1d& basis);

template<typename Alg, typename TypeMon>
void Alg2IndicesP(const Alg& alg, const std::vector<TypeMon>& basis, int1d& result)
{
    result.clear();
    for (const auto& m : alg.data) {
        auto p = std::lower_bound(basis.begin(), basis.end(), m);
#ifdef MYDEBUG
        if (p == basis.end() || m < (*p)) {
            throw ErrorIdMsg(0x57f14e21U, "Index not found");
        }
#endif
        result.push_back(int(p - basis.begin()));
    }
}

template <typename Alg, typename Basis>
void Alg2IndicesP(const Alg& alg, const Basis& basis, AdamsDeg deg, int1d& result)
{
    if (alg)
        Alg2IndicesP(alg, basis.at(deg), result);
	else
		result.clear();
}

template <typename Alg, typename TypeMon>
void Indices2AlgP(const int1d& indices, const std::vector<TypeMon>& basis, Alg& result)
{
    result.data.clear();
    for (int i : indices)
        result.data.push_back(basis[i]);
}

template <typename Alg, typename Basis>
void Indices2AlgP(const int1d& indices, const Basis& basis, AdamsDeg deg, Alg& result)
{
    if (indices.empty())
		result.data.clear();
	else
		Indices2AlgP(indices, basis.at(deg), result);
}

template <typename T>
class Basis
{
private:
    std::vector<std::vector<T>> data_;

public:
    bool has(AdamsDeg deg) const {
        int n = deg.stem();
        return n < (int)data_.size() && deg.s < (int)data_[n].size();
    }
    auto& at(AdamsDeg deg) const
    {
        return data_[deg.stem()][deg.s];
    }
    auto& at(AdamsDeg deg, const T& default_) const
    {
        auto n = deg.stem();
        if (n >= (int)data_.size() || deg.s >= (int)data_[n].size())
            return default_;
        return data_[deg.stem()][deg.s];
    }
    auto& operator[](AdamsDeg deg)
    {
        return ut::get(ut::get(data_, deg.stem()), deg.s);
    }

    AdamsDeg1d degs() const
    {
        AdamsDeg1d result;
        for (int n = 0; n < (int)data_.size(); ++n)
            for (int s = 0; s < (int)data_[n].size(); ++s)
                result.push_back(AdamsDeg(s, n + s));
        return result;
    }
};

using BasisMon = Basis<Mon1d>;
using BasisMMod = Basis<MMod1d>;

}  // namespace alg2

#endif /* ALGEBRAS_H */