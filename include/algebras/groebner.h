/** \file groebner.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
#include "linalg.h"
#include <execution>
#include <map>
#include <unordered_map>
#include <memory>

/* extension of namespace alg in algebras.h */
namespace alg {

namespace detail {
    void mul(const Mon& mon1, const MonOnStack& mon2, MonOnStack& result);
    void mul(const Mon& mon1, const Mon& mon2, MonOnStack& result);
    void div(const Mon& mon1, const Mon& mon2, MonOnStack& result);

    /*
     * Assuming that f[i] is divisable by g[0],
     * this function replaces f with f + g * (f[i]/g[0]).
     *
     * We tried our best to reduce the memory allocations.
     * TODO: use log
     */
    template <typename FnCmp>
    void Reduce(Polynomial<FnCmp>& f, const Polynomial<FnCmp>& g, const size_t index)
    {
        Mon1d h;
        h.reserve(f.data.size() + g.data.size() - (size_t)index - 2);
        MonOnStack q;
        div(f.data[index], g.data[0], q);
        size_t i = index + 1, j = 1;
        MonOnStack prod;
        bool updateProd = true;
        while (j < g.data.size()) {
            if (i == f.data.size()) {
                f.data.resize(index + h.size() + g.data.size() - j);
                for (size_t k = 0; k < h.size(); ++k)
                    f.data[index + k] = std::move(h[k]);
                size_t index1 = index + h.size();
                for (size_t k = j; k < g.data.size(); ++k)
                    mul(g.data[k], q, f.data[index1 + k - j]);
                return;
            }
            if (updateProd)
                mul(g.data[j], q, prod);
            if (FnCmp::template cmp(f.data[i], prod)) {
                updateProd = false;
                h.push_back(std::move(f.data[i++]));
            }
            else {
                updateProd = true;
                ++j;
                if (FnCmp::template cmp(prod, f.data[i]))
                    h.push_back(Mon(prod));
                else
                    ++i;
            }
        }
        size_t index1 = index + h.size();
        size_t new_size = index1 + f.data.size() - i;
        if (new_size > f.data.size()) {
            size_t old_size = f.data.size();
            f.data.resize(new_size);
            for (size_t k = old_size; k-- > i;)
                f.data[index1 + k - i] = std::move(f.data[k]);
        }
        else if (new_size < f.data.size()) {
            for (size_t k = i; k < f.data.size(); ++k)
                f.data[index1 + k - i] = std::move(f.data[k]);
            f.data.erase(f.data.begin() + new_size, f.data.end());
        }
        for (size_t k = 0; k < h.size(); ++k)
            f.data[index + k] = std::move(h[k]);
    }

    /*
     * Return if `mon1` and `mon2` have a nontrivial common factor.
     */
    bool HasGCD(const Mon& mon1, const Mon& mon2);

    /*
     * If GCD(mon1, mon2) has d_min <= deg <= d_max, then mon_out = GCD() and return true.
     * Otherwise return false and mon_out is undefined.
     */
    bool GCD(const Mon& mon1, const Mon& mon2, MonOnStack& result, const array& gen_degs, int& deg, int d_min, int d_max);

}  // namespace detail

template <typename FnCmp>
class Groebner
{
private:
    using TypeIndexKey = int;
    using TypeIndex = std::unordered_map<TypeIndexKey, std::vector<int>>;
    using TypeTrace = uint64_t;
public:
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;

public:
    Poly1d data; // TODO: make it private

private:
    std::vector<TypeTrace> traces_;
    TypeIndex index_;

public:
    Groebner() = default;
    Groebner(Poly1d polys) : data(std::move(polys))
    {
        for (int i = 0; i < (int)data.size(); ++i) {
            index_[Key(data[i].GetLead())].push_back(i);
            traces_.push_back(Trace(data[i].GetLead()));
        }
    }

    template <typename FnCmp1>
    Groebner<FnCmp1> ExtendMO() const /* Extend the monomial ordering */
    {
        Groebner<FnCmp1> gb1;
        for (const auto& p : data)
            gb1.push_back(typename Groebner<FnCmp1>::Poly{p.data});
        return gb1;
    }

private:
    static TypeIndexKey Key(const Mon& lead)
    {
       return TypeIndexKey{lead.size() == 1 ? lead.back().gen : lead.back().gen + 0x10000 * lead[lead.size() - 2].gen };
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(const Mon& mon) const
    {
        auto t = Trace(mon);
        for (int i = 0; i < (int)mon.size(); ++i) {
            for (int j = -1; j < i; ++j) {
                auto key = TypeIndexKey{j == -1 ? mon[i].gen : mon[i].gen + 0x10000 * mon[j].gen};
                auto p = index_.find(key);
                if (p != index_.end()) {
                    for (int k : p->second) {
                        if (t >= traces_[k] && !(traces_[k] & (t - traces_[k])))
                            if (divisible(data[k].GetLead(), mon))
                                return (int)k;
                    }
                }
            }
        }
        return -1;
    }

public:
    Poly Reduce(Poly poly) const
    {
        size_t index = 0;
        while (index < poly.data.size()) {
            int gb_index = IndexOfDivisibleLeading(poly.data[index]);
            if (gb_index != -1) {
                alg::detail::Reduce(poly, data[gb_index], index);
            }
            else
                ++index;
        }
        return poly;
    }

    /* Leadings[i] is the monomials that end with generator i. */
    Mon2d GetLeadings(size_t gens_size) const
    {
        alg::Mon2d result;
        result.resize(gens_size);
        for (size_t i = 0; i < data.size(); ++i)
            result[data[i].GetLead().back().gen].push_back(data[i].GetLead());
        return result;
    }

public: /* Convenient interfaces for the member `gb` */
    auto size() const
    {
        return data.size();
    }
    auto& operator[](size_t index) const
    {
        return data[index];
    }
    void push_back(Poly g)
    {
        index_[Key(g.GetLead())].push_back((int)data.size());
        traces_.push_back(Trace(g.GetLead()));
        data.push_back(std::move(g));
    }
    bool operator==(const Groebner<FnCmp>& rhs) const
    {
        return data == rhs.data;
    }
    static TypeTrace Trace(const Mon& lead)
    {
        TypeTrace result = 0;
        for (size_t i = 0; i < lead.size(); ++i) {
            const int bits_exp1 = 56;
            const int bits_exp2 = 64 - bits_exp1;
            result |= (TypeTrace(1) << (lead[i].gen % bits_exp1));
            if (lead[i].exp >= 2)
                result |= (TypeTrace(1) << ((lead[i].gen % bits_exp2) + bits_exp1));
        }
        return result;
    }
    TypeTrace traces(size_t i) const
    {
        return traces_[i];
    }
};

namespace detail {
    struct GbBufferEle
    {
        Mon t1, t2;
        int i1, i2;

        /* Return t1 * gb[i1] + t2 * gb[i2]  */
        template <typename FnCmp>
        Polynomial<FnCmp> GetPoly(const Groebner<FnCmp>& gb) const
        {
            Polynomial<FnCmp> result;
            MonOnStack prod1, prod2;
            result.data.reserve(gb[i1].data.size() + gb[i2].data.size() - 2);
            auto p1 = gb[i1].data.begin() + 1, p2 = gb[i2].data.begin() + 1;
            while (p1 != gb[i1].data.end() && p2 != gb[i2].data.end()) {
                mul(*p1, t1, prod1);
                mul(*p2, t2, prod2);
                if (FnCmp::template cmp(prod1, prod2)) {
                    result.data.push_back(Mon(prod1));
                    ++p1;
                }
                else if (FnCmp::template cmp(prod2, prod1)) {
                    result.data.push_back(Mon(prod2));
                    ++p2;
                }
                else {
                    ++p1;
                    ++p2;
                }
            }
            while (p1 != gb[i1].data.end()) {
                mul(*p1++, t1, prod1);
                result.data.push_back(Mon(prod1));
            }
            while (p2 != gb[i2].data.end()) {
                mul(*p2++, t2, prod2);
                result.data.push_back(Mon(prod2));
            }
            return result;
        }
    };
}
using GroebnerLex = Groebner<CmpLex>;
using GroebnerRevlex = Groebner<CmpRevlex>;

using GbBuffer = std::map<int, std::vector<detail::GbBufferEle>>;

/**
 * Comsume relations from `buffer` in degree `<= deg`
 * while adding new relations back to `buffer` in degree `<= deg_max`.
 * `deg=-1` or `deg_max=-1` means infinity.
 */
template <typename FnCmp, typename FnPred, typename FnDeg>
void AddRels(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, FnPred pred, FnDeg _gen_deg, int deg, int deg_max)
{
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using PPoly1d = std::vector<const Poly*>;

    std::map<int, PPoly1d> rels_graded;
    for (const auto& rel : rels) {
        if (rel) {
            int d = TplGetDeg(rel.GetLead(), _gen_deg);
            if (d <= deg || deg == -1) {
                rels_graded[d].push_back(&rel);
                if (buffer.find(d) == buffer.end())
                    buffer[d] = {};
            }
        }
    }
    std::vector<std::pair<int, detail::GbBufferEle>> buffer_new_rels;
    auto p_buffer = buffer.begin();
    for (; p_buffer != buffer.end() && (p_buffer->first <= deg || deg == -1); ++p_buffer) {
        /* Reduce relations from buffer in degree `p_buffer->first` */
        const int d = p_buffer->first;
        // std::cout << "d=" << d << '\n';  ////
        Poly1d rels_d;
        bool bRels = rels_graded.find(d) != rels_graded.end();
        Poly1d rels_tmp(p_buffer->second.size() + (bRels ? rels_graded.at(d).size() : 0));
        int buffer_size = (int)p_buffer->second.size();
        ut::Range range(0, buffer_size);
        std::for_each(std::execution::par_unseq, range.begin(), range.end(), [&gb, p_buffer, &rels_tmp](int i) { rels_tmp[i] = gb.Reduce(p_buffer->second[i].GetPoly(gb)); });
        if (bRels) {
            auto& p_rels = rels_graded.at(d);
            ut::Range range1(0, (int)p_rels.size());
            std::for_each(std::execution::seq, range1.begin(), range1.end(), [&gb, &p_rels, buffer_size, &rels_tmp](int i) { rels_tmp[buffer_size + (size_t)i] = gb.Reduce(*p_rels[i]); });
        }

        /* Triangulate these relations */
        for (auto& rel : rels_tmp) {
            for (const Poly& rel1 : rels_d)
                if (std::binary_search(rel.data.begin(), rel.data.end(), rel1.GetLead(), FnCmp::template cmp<Mon, Mon>))
                    rel += rel1;
            if (rel)
                rels_d.push_back(std::move(rel));
        }

        /* Add these relations */
        buffer_new_rels.resize(gb.size() + rels_d.size());
        for (auto& rel : rels_d) {
            ut::Range range2(0, (int)gb.size());
            auto trace_rel = gb.Trace(rel.GetLead());
            std::for_each(std::execution::par_unseq, range2.begin(), range2.end(), [&gb, &buffer_new_rels, &rel, trace_rel, pred, _gen_deg, d, deg_max](int i) {
                if ((trace_rel & gb.traces(i)) && detail::HasGCD(rel.GetLead(), gb[i].GetLead()) && pred(rel.GetLead(), gb[i].GetLead())) {
                    /* detail::MonOnStack gcd;
                    detail::GCD(rel.GetLead(), gb[i].GetLead(), gcd, ) */
                    Mon gcd = GCD(rel.GetLead(), gb[i].GetLead());
                    Mon t1 = div(rel.GetLead(), gcd);
                    Mon t2 = div(gb[i].GetLead(), gcd);
                    int deg_new_rel = d + TplGetDeg(gb[i].GetLead(), _gen_deg) - TplGetDeg(gcd, _gen_deg);
                    if (deg_new_rel <= deg_max || deg_max == -1)
                        buffer_new_rels[i] = std::make_pair(deg_new_rel, detail::GbBufferEle{std::move(t1), std::move(t2), i, (int)gb.size()});
                    else
                        buffer_new_rels[i].first = -1;
                }
                else
                    buffer_new_rels[i].first = -1;
            });
            for (size_t i = 0; i < gb.size(); ++i) {
                if (buffer_new_rels[i].first != -1)
                    buffer[buffer_new_rels[i].first].push_back(std::move(buffer_new_rels[i].second));
            }
            gb.push_back(std::move(rel));
        }
    }
    buffer.erase(buffer.begin(), p_buffer);
}

template <typename FnCmp>
void AddRels(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, const array& gen_degs, int deg, int deg_max)
{
    AddRels(
        gb, buffer, rels, [](const Mon&, const Mon&) { return true; }, [&gen_degs](int i) { return gen_degs[i]; }, deg, deg_max);
}

template <typename FnCmp>
void AddRels(Groebner<FnCmp>& gb, const std::vector<Polynomial<FnCmp>>& rels, const array& gen_degs, int deg)
{
    GbBuffer buffer;
    AddRels(gb, buffer, rels, gen_degs, deg, deg);
}

/**
 * A fast algorithm that computes
 * `poly ** n` modulo `gb`
 */
template <typename FnCmp>
Polynomial<FnCmp> pow(const Polynomial<FnCmp>& poly, int n, const Groebner<FnCmp>& gb)
{
    using Poly = Polynomial<FnCmp>;
    Poly result = Poly::Unit();
    if (n == 0)
        return result;
    Poly power = poly;
    while (n) {
        if (n % 2 != 0)
            result = gb.Reduce(result * power);
        n >>= 1;
        if (n)
            power = gb.Reduce(power.Square());
    }
    return result;
}

/**
 * Replace the generators in `poly` with elements given in `map`
 * and evaluate modulo `gb`.
 * @param poly The polynomial to be substituted.
 * @param map `map(i)` is the polynomial that substitutes the generator of id `i`.
 * @param gb A Groebner basis.
 */
template <typename FnCmp, typename FnMap>
Polynomial<FnCmp> SubsMGb(const Mon1d& data, FnMap map, const Groebner<FnCmp>& gb)
{
    using Poly = Polynomial<FnCmp>;
    Poly result;
    for (const Mon& m : data) {
        Poly fm = Poly::Unit();
        for (MonInd p = m.begin(); p != m.end(); ++p)
            fm = gb.Reduce(fm * pow(map(p->gen), p->exp, gb));
        result += fm;
    }
    return result;
}

/**
 * Specialization.
 * @see SubsMGb(const Poly&, FnType, const GbType&)
 */
template <typename FnCmp>
Polynomial<FnCmp> SubsMGb(const Mon1d& data, const std::vector<Polynomial<FnCmp>>& map, const Groebner<FnCmp>& gb)
{
    return SubsMGb(
        data, [&map](int i) { return map[i]; }, gb);
}

/**
 * Generate buffer in degree `t_min <= t <= t_max`
 */
template <typename FnCmp>
GbBuffer GenerateBuffer(const alg::Groebner<FnCmp>& gb, const array& gen_degs, int d_min, int d_max)
{
    GbBuffer buffer;
    Mon gcd;
    int deg_gcd;
    for (size_t i = 0; i < gb.data.size(); ++i) {
        int deg1 = gb[i].GetDeg(gen_degs);
        for (size_t j = i + 1; j < gb.data.size(); ++j) {
            int deg2 = gb[j].GetDeg(gen_degs);
            if (d_min < deg1 + deg2 && std::max(deg1, deg2) <= d_max) {
                if (detail::GCD(gb[i].GetLead(), gb[j].GetLead(), gcd, gen_degs, deg_gcd, std::max(deg1 + deg2 - d_max, 1), deg1 + deg2 - d_min)) {
                    buffer[deg_gcd].push_back(detail::GbBufferEle{gcd, (int)i, (int)j});
                }
            }
        }
    }
    return buffer;
}

/**********************************************************
 * Algorithms that use Groebner basis
 **********************************************************/

inline constexpr int GEN_IDEAL = 0x4000000;
inline constexpr int GEN_SUB = 0x2000000;

namespace detail {
    /* For monomials in the homological groebner basis get the degree of the first part consisting of orginal generators and set `p` to be
     * the end location of the first part.
     */
    template <typename TypeMonIter>
    int GetPart1Deg(const TypeMonIter mon_begin, const TypeMonIter mon_end, const array& gen_degs, TypeMonIter& p)
    {
        int result = 0;
        for (p = mon_begin; p != mon_end && !(p->gen & GEN_SUB); ++p)
            result += gen_degs[p->gen] * p->exp;
        return result;
    };

}
/**
 * Revlex on extra generators and FnCmp on the rest
 */
template <typename FnCmp>
struct CmpIdeal
{
    using submo = FnCmp;
    static constexpr std::string_view name = "Ideal";
    template <typename Type1, typename Type2>
    static bool cmp(const Type1& m1, const Type2& m2)
    {
        auto mid1 = std::lower_bound(m1.begin(), m1.end(), GEN_IDEAL, [](const GenPow& p, int g) { return p.gen < g; });
        auto mid2 = std::lower_bound(m2.begin(), m2.end(), GEN_IDEAL, [](const GenPow& p, int g) { return p.gen < g; });
        if (CmpRevlex::cmp_ranges(mid1, m1.end(), mid2, m2.end()))
            return true;
        else if (std::equal(mid1, m1.end(), mid2, m2.end())) {
            if (FnCmp::template cmp_ranges(m1.begin(), mid1, m2.begin(), mid2))
                return true;
        }
        return false;
    }
};

/*
 * The monomial ordering that computes the homology.
 */
template <typename FnCmp>
struct CmpHomology
{
    using submo = FnCmp;
    static constexpr std::string_view name = "Homology";
    inline static array gen_degs;
    template <typename Type1, typename Type2>
    static bool cmp(const Type1& m1, const Type2& m2)
    {
        typename Type1::const_iterator m1_mid1;
        typename Type2::const_iterator m2_mid1;
        int d1 = detail::template GetPart1Deg(m1.begin(), m1.end(), gen_degs, m1_mid1);
        int d2 = detail::template GetPart1Deg(m2.begin(), m2.end(), gen_degs, m2_mid1);
        if (d1 > d2)
            return true;
        if (d1 < d2)
            return false;
        
        if (FnCmp::template cmp_ranges(m1.begin(), m1_mid1, m2.begin(), m2_mid1))
            return true;
        if (FnCmp::template cmp_ranges(m2.begin(), m2_mid1, m1.begin(), m1_mid1))
            return false;

        auto m1_mid2 = std::lower_bound(m1_mid1, m1.end(), GEN_IDEAL + GEN_SUB, [](const GenPow& p, int g) { return p.gen < g; });
        auto m2_mid2 = std::lower_bound(m2_mid1, m2.end(), GEN_IDEAL + GEN_SUB, [](const GenPow& p, int g) { return p.gen < g; });
        if (CmpRevlex::cmp_ranges(m1_mid2, m1.end(), m2_mid2, m2.end()))
            return true;
        if (CmpRevlex::cmp_ranges(m2_mid2, m2.end(), m1_mid2, m1.end()))
            return false;

        if (CmpRevlex::cmp_ranges(m1_mid1, m1_mid2, m2_mid1, m2_mid2))
            return true;
        return false;
    }
};

namespace detail {
    template <typename FnCmp>
    void AddRelsIdeal(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, const array& gen_degs, const array& gen_degs_y, int deg, int deg_max)
    {
        AddRels(
            gb, buffer, rels, [](const Mon& m1, const Mon& m2) { return !(m1.back().gen & GEN_IDEAL) && !(m2.back().gen & GEN_IDEAL); },
            [&gen_degs, &gen_degs_y](int i) {
                return ((i & GEN_IDEAL) ? gen_degs_y[i - GEN_IDEAL] : gen_degs[i]);
            },
            deg, deg_max);
    }
    template <typename FnCmp>
    void AddRelsModule(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, const array& gen_degs, const array& gen_degs_y, int deg, int deg_max)
    {
        AddRels(
            gb, buffer, rels, [](const Mon& m1, const Mon& m2) { return !((m1.back().gen & GEN_IDEAL) && (m2.back().gen & GEN_IDEAL) && (m1.back().gen != m2.back().gen)); },
            [&gen_degs, &gen_degs_y](int i) {
                return ((i & GEN_IDEAL) ? gen_degs_y[i - GEN_IDEAL] : gen_degs[i]);
            },
            deg, deg_max);
    }
    template <typename FnCmp>
    void AddRelsHomology(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, const MayDeg1d& gen_degs, const MayDeg1d& gen_degs_x, const MayDeg1d& gen_degs_b, int deg, int deg_max)
    {
        AddRels(
            gb, buffer, rels, [](const Mon& m1, const Mon& m2) { return true; },
            [&gen_degs, &gen_degs_x, &gen_degs_b](int i) {
                if (i & GEN_IDEAL)
                    return gen_degs_b[i - GEN_IDEAL - GEN_SUB].t;
                else if (i & GEN_SUB)
                    return gen_degs_x[i - GEN_SUB].t;
                else
                    return gen_degs[i].t;
            },
            deg, deg_max);
    }
}  // namespace detail

/**
 * Compute the minimal generating set of `vectors` inplace.
 *
 * Each element of vectors is considered as an element of $R^n$ where $R$ is the algebra
 * determined by the Groebner basis `gb`.
 *
 * The vectors should all be nontrivial.
 *
 * The degree of the basis of $R^n$ is determined by `basis_degs`.
 */
template <typename FnCmp>
std::vector<std::vector<Polynomial<FnCmp>>>& Indecomposables(const Groebner<FnCmp>& gb, std::vector<std::vector<Polynomial<FnCmp>>>& vectors, const array& gen_degs, const array& basis_degs)
{
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner<FnCmp>;
    using PolyI = Polynomial<CmpIdeal<FnCmp>>;
    using PolyI1d = std::vector<PolyI>;
    using GbI = alg::Groebner<CmpIdeal<FnCmp>>;

    if (vectors.empty())
        return vectors;

    /* Convert each vector v into a relation \\sum v_iy_i */
    PolyI1d rels;
    array degs;
    for (auto& v : vectors) {
        PolyI rel;
        for (size_t i = 0; i < basis_degs.size(); ++i)
            if (v[i])
                rel += PolyI::Sort((v[i] * GenPow::Mon(GEN_IDEAL + (int)i)).data);
        for (size_t i = 0; i < basis_degs.size(); ++i) {
            if (v[i]) {
                degs.push_back(v[i].GetDeg(gen_degs) + basis_degs[i]);
                break;
            }
        }
        if (rel)
            rels.push_back(std::move(rel));
        else
            v.clear();
    }
    ut::RemoveEmptyElements(vectors);
    if (!vectors.empty()) {
        array indices = ut::range((int)vectors.size());
        std::sort(indices.begin(), indices.end(), [&degs](int i, int j) { return degs[i] < degs[j]; });

        /* Add relations ordered by degree to gb1 */
        GbI gbI = gb.ExtendMO<CmpIdeal<FnCmp>>();
        GbBuffer buffer;
        int deg_max = degs[indices.back()];
        for (int i : indices) {
            detail::AddRelsModule(gbI, buffer, {}, gen_degs, basis_degs, degs[i], deg_max);
            PolyI rel = gbI.Reduce(rels[i]);
            if (rel)
                detail::AddRelsModule(gbI, buffer, {std::move(rel)}, gen_degs, basis_degs, degs[i], deg_max);
            else
                vectors[i].clear();
        }

        /* Keep only the indecomposables in `vectors` */
        ut::RemoveEmptyElements(vectors);
    }
    return vectors;
}

/**
 * Compute the generating set of linear relations among `polys`.
 *
 * The result is truncated by `deg<=deg_max`.
 */
template <typename FnCmp>
std::vector<std::vector<Polynomial<FnCmp>>> AnnSeq(const Groebner<FnCmp>& gb, const std::vector<Polynomial<FnCmp>>& polys, const array& gen_degs, int deg_max)
{
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Poly2d = std::vector<Poly1d>;
    using Gb = alg::Groebner<FnCmp>;
    using PolyI = Polynomial<CmpIdeal<FnCmp>>;
    using PolyI1d = std::vector<PolyI>;
    using GbI = alg::Groebner<CmpIdeal<FnCmp>>;

    Poly2d result;
    if (polys.empty())
        return result;
    PolyI1d rels;
    array gen_degs_y;
    int n = (int)polys.size();

    /* Add relations Xi=polys[i] to gb to obtain gb1 */
    for (int i = 0; i < n; ++i) {
        PolyI p{polys[i].data};
        gen_degs_y.push_back(p.GetDeg(gen_degs));
        p += PolyI::Gen(GEN_IDEAL + i);
        rels.push_back(std::move(p));
    }
    GbI gbI = gb.ExtendMO<CmpIdeal<FnCmp>>();
    GbBuffer buffer;
    detail::AddRelsIdeal(gbI, buffer, rels, gen_degs, gen_degs_y, deg_max, deg_max);

    /* Extract linear relations from gb1 */
    for (const PolyI& g : gbI.data) {
        if (g.GetLead().back().gen & GEN_IDEAL) {
            Poly1d ann;
            ann.resize(n);
            for (const Mon& m : g.data) {
                auto mid = std::lower_bound(m.begin(), m.end(), GEN_IDEAL, [](const GenPow& p, int g) { return p.gen < g; });
                Mon m1(m.begin(), mid), m2(mid, m.end());  // TODO: improve
                ann[m2.front().gen - GEN_IDEAL] += gb.Reduce(SubsMGb(
                                                                 {div(m2, GenPow::Mon(m2.front().gen))}, [&polys](int i) { return polys[i - GEN_IDEAL]; }, gb)
                                                             * m1);
            }
            result.push_back(std::move(ann));
        }
    }

    /* Add commutators to linear relations */
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (deg_max == -1 || gen_degs_y[i] + gen_degs_y[j] <= deg_max) {
                Poly1d result_i;
                result_i.resize(n);
                result_i[i] = polys[j];
                result_i[j] = polys[i];
                result.push_back(std::move(result_i));
            }
        }
    }

    Indecomposables(gb, result, gen_degs, gen_degs_y);
    return result;
}

template <typename FnCmp>
std::map<MayDeg, Mon1d> ExtendBasis(const Mon2d& leadings, const MayDeg1d& gen_degs, std::map<MayDeg, Mon1d>& basis, int t)
{
    std::map<MayDeg, Mon1d> basis_t;
    if (t == 0) {
        basis[MayDeg{0, 0, 0}].push_back({});
        return basis_t;
    }
    for (int gen_id = 0; gen_id < (int)gen_degs.size(); ++gen_id) {
        int t1 = t - gen_degs[gen_id].t;
        if (t1 >= 0) {
            auto basis_t1_begin = basis.lower_bound(MayDeg{0, t1, 0});
            auto basis_t1_end = basis.lower_bound(MayDeg{0, t1 + 1, 0});
            for (auto p_basis = basis_t1_begin; p_basis != basis_t1_end; ++p_basis) {
                for (auto p_m = p_basis->second.begin(); p_m != p_basis->second.end(); ++p_m) {
                    if (p_m->empty() || gen_id >= p_m->back().gen) {
                        Mon mon = mul(*p_m, GenPow::Mon(gen_id));
                        if (std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(), [&mon](const Mon& _m) { return divisible(_m, mon); }))
                            basis_t[p_basis->first + gen_degs[gen_id]].push_back(std::move(mon));
                    }
                }
            }
        }
    }

    for (auto& p : basis_t)
        std::sort(p.second.begin(), p.second.end(), FnCmp::template cmp<Mon, Mon>);
    return basis_t;
}

template <typename FnCmp>
void Homology(const Groebner<FnCmp>& gb, const MayDeg1d& gen_degs, /*std::vector<std::string>& gen_names, */ const std::vector<Polynomial<FnCmp>>& gen_diffs, const MayDeg& deg_diff, GroebnerRevlex& gb_h, MayDeg1d& gen_degs_h,
              std::vector<std::string>& gen_names_h, std::vector<Polynomial<FnCmp>>& gen_repr_h, int t_max)
{
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using PolyH = Polynomial<CmpHomology<FnCmp>>;
    using PolyH1d = std::vector<PolyH>;
    using Gb = Groebner<FnCmp>;
    using GbH = Groebner<CmpHomology<FnCmp>>;

    Gb gb_AoAZ = gb;
    GbBuffer buffer_AoAZ;
    Mon2d leadings_AoAZ = gb_AoAZ.GetLeadings(gen_degs.size());
    size_t gb_AoAZ_old_size = gb_AoAZ.size();
    std::map<MayDeg, Mon1d> basis_AoAZ;
    std::vector<const Mon*> basis_AoAZ_flat;
    std::vector<Poly> diff_basis_AoAZ_flat;
    basis_AoAZ[MayDeg{0, 0, 0}].push_back({});

    CmpHomology<FnCmp>::gen_degs.clear();
    for (const auto& deg : gen_degs)
        CmpHomology<FnCmp>::gen_degs.push_back(deg.t);
    Groebner<CmpHomology<FnCmp>> gb_H = gb.ExtendMO<CmpHomology<FnCmp>>();
    GbBuffer buffer_h;
    MayDeg1d gen_degs_b;
    Poly1d gen_reprs_b; /* dm_i = b_i */

    int deg_t_allX = 0;
    for (auto& deg : gen_degs)
        deg_t_allX += deg.t;
    //std::cout << "deg_t_allX=" << deg_t_allX << '\n';

    for (int t = 1; t_max == -1 || t <= t_max; ++t) {
        //std::cout << "t=" << t << '\n';  //
        PolyH1d rels_h_t;
        Poly1d rels_AoAZ_t;
        if (t <= deg_t_allX) {
            /* Find m_i, b_i and x_i */
            for (size_t i = gb_AoAZ_old_size; i < gb_AoAZ.data.size(); ++i)
                leadings_AoAZ[gb_AoAZ[i].GetLead().back().gen].push_back(gb_AoAZ[i].GetLead());
            gb_AoAZ_old_size = gb_AoAZ.size();
            std::map<MayDeg, Mon1d> basis_AoAZ_t = ExtendBasis<FnCmp>(leadings_AoAZ, gen_degs, basis_AoAZ, t);
            for (auto& [deg, b_deg] : basis_AoAZ_t) {
                array2d map_diff;
                for (const Mon& m : b_deg) {
                    Poly diff_m = gb.Reduce(GetDiff(m, gen_diffs));
                    map_diff.push_back(alg::hash1d(diff_m.data));
                }
                array2d image_diff, kernel_diff, g_diff;
                lina::SetLinearMap(map_diff, image_diff, kernel_diff, g_diff);

                /* Add x_i and the relations for x_i */
                for (const array& k : kernel_diff) {
                    gen_names_h.push_back("x_{" + std::to_string(gen_degs_h.size()) + "}");
                    gen_degs_h.push_back(deg);
                    gen_repr_h.push_back({Indices2Poly(k, b_deg)});
                    rels_AoAZ_t.push_back(gen_repr_h.back());
                    rels_h_t.push_back(PolyH({gen_repr_h.back().data}) + PolyH::Gen(GEN_SUB + (int)gen_degs_h.size() - 1));
                    b_deg[k.front()].clear();
                }
                ut::RemoveEmptyElements(b_deg);

                /* Add b_i, m_i and the relation dm_i = b_i */
                std::vector<Mon1d::iterator> toBeRemoved;
                MayDeg d_dm = deg + deg_diff;
                for (const Mon& m : b_deg) {
                    basis_AoAZ_flat.push_back(&m);
                    diff_basis_AoAZ_flat.push_back(gb.Reduce(GetDiff(m, gen_diffs)));
                    auto& dm = diff_basis_AoAZ_flat.back();
                    if (dm) {
                        rels_AoAZ_t.push_back(dm);
                        auto lead_dm = dm.GetLead();
                        auto ptr = std::lower_bound(basis_AoAZ_t[d_dm].begin(), basis_AoAZ_t[d_dm].end(), lead_dm, FnCmp::template cmp<Mon, Mon>);
                        if (ptr != basis_AoAZ_t[d_dm].end() && (*ptr == lead_dm))
                            toBeRemoved.push_back(ptr);
                    }
                    gen_degs_b.push_back(deg);
                    rels_h_t.push_back(PolyH({dm.data}) + PolyH::Gen(GEN_SUB + GEN_IDEAL + (int)gen_degs_b.size() - 1));
                }
                if (!toBeRemoved.empty()) {
                    for (auto p : toBeRemoved)
                        p->clear();
                    ut::RemoveEmptyElements(basis_AoAZ_t[d_dm]);
                }
            }
            basis_AoAZ.merge(basis_AoAZ_t);
            AddRels(gb_AoAZ, buffer_AoAZ, rels_AoAZ_t, CmpHomology<FnCmp>::gen_degs, t, t_max);
            rels_AoAZ_t.clear();
        }
        else {
            if (buffer_h.empty())
                break;
        }
        /* Find  y_i */
        size_t i_start_t = gb_H.size();
        detail::AddRelsHomology(gb_H, buffer_h, rels_h_t, gen_degs, gen_degs_h, gen_degs_b, t, t_max);
        rels_h_t.clear();

        for (size_t i = i_start_t; i < gb_H.size(); ++i) {
            const auto& lead = gb_H[i].GetLead();
            if ((lead.front().gen & GEN_SUB) && (lead.back().gen & GEN_IDEAL) && (lead.back().exp == 1) && (lead.size() == 1 || !(lead[lead.size() - 2].gen & GEN_IDEAL))) {
                PolyH y;
                for (const Mon& m : gb_H[i].data) {
                    auto mid = std::lower_bound(m.begin(), m.end(), GEN_IDEAL + GEN_SUB, [](const GenPow& p, int g) { return p.gen < g; });
                    Mon m1(m.begin(), mid), m2(mid, m.end());  // TODO: improve if bottlenet
                    y += gb_H.Reduce(PolyH::Mon_(mul(mul(m1, div(m2, GenPow::Mon(m2.front().gen))), *basis_AoAZ_flat[m2.front().gen - GEN_IDEAL - GEN_SUB])));
                }
                if (y && !(y.GetLead().front().gen & GEN_SUB)) {
                    Poly repr_y = SubsMGb(
                        y.data, [&gen_repr_h, &diff_basis_AoAZ_flat](int i) { return ((i & GEN_IDEAL) ? diff_basis_AoAZ_flat[i - GEN_IDEAL - GEN_SUB] : ((i & GEN_SUB) ? gen_repr_h[i - GEN_SUB] : Poly::Gen(i))); }, gb);
                    detail::AddRelsHomology(gb_H, buffer_h, {y + PolyH::Gen(GEN_SUB + (int)gen_degs_h.size())}, gen_degs, gen_degs_h, gen_degs_b, t, t_max);
                    gen_names_h.push_back("y_{" + std::to_string(gen_degs_h.size()) + "}");
                    gen_degs_h.push_back(repr_y.GetMayDeg(gen_degs));
                    gen_repr_h.push_back(std::move(repr_y));
                    rels_AoAZ_t.push_back(gen_repr_h.back());
                }
            }
        }
        AddRels(gb_AoAZ, buffer_AoAZ, rels_AoAZ_t, CmpHomology<FnCmp>::gen_degs, t, t_max);
        rels_AoAZ_t.clear();
    }
    /* Prepare relations for homology */
    for (auto& p : gb_H.data) {
        if (p && (p.GetLead().front().gen & GEN_SUB)) {
            auto p_m = std::lower_bound(p.data.begin(), p.data.end(), 1, [](const Mon& m, int) { return !(m.back().gen & GEN_IDEAL); });
            p.data.erase(p_m, p.data.end());
            if (p)
                gb_h.push_back(alg::subs<alg::CmpRevlex>(p.data, [](int i) { return PolyRevlex::Gen(i - GEN_SUB); }));
        }
    }
}

#define TEMPLATE_EXAMPLES
#ifdef TEMPLATE_EXAMPLES
namespace template_examples {
    using FnCmp = alg::CmpLex;
    using FnPred = bool (*)(alg::Mon, alg::Mon);
    using FnDeg = int (*)(int);

    inline void AddRels_(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, FnPred pred, FnDeg _get_deg, int deg, int deg_max)
    {
        return alg::AddRels(gb, buffer, rels, pred, _get_deg, deg, deg_max);
    }

    inline void Homology_(const Groebner<FnCmp>& gb, const MayDeg1d& gen_degs, /*std::vector<std::string>& gen_names, */ const std::vector<Polynomial<FnCmp>>& gen_diffs, const MayDeg& deg_diff, GroebnerRevlex& gb_h, MayDeg1d& gen_degs_h,
                          std::vector<std::string>& gen_names_h, std::vector<Polynomial<FnCmp>>& gen_repr_h, int t_max)
    {
        Homology(gb, gen_degs, gen_diffs, deg_diff, gb_h, gen_degs_h, gen_names_h, gen_repr_h, t_max);
    }
}  // namespace template_examples
#endif

} /* namespace alg */

#endif /* GROEBNER_H */
