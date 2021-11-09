/** \file groebner.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
#include <execution>
#include <map>
#include <memory>

/* extension of namespace alg in algebras.h */
namespace alg {

namespace detail {
    /*
     * This version of `mul` is designed to reduce the allocations of memory.
     */
    void mul(const Mon& mon1, const Mon& mon2, Mon& mon_out);

    /*
     * Assuming that f[i] is divisable by g[0],
     * this function replaces f with f + g * (f[i]/g[0]).
     *
     * We tried our best to reduce the memory allocations.
     */
    template <typename FnCmp>
    void Reduce(Polynomial<FnCmp>& f, const Polynomial<FnCmp>& g, const size_t index)
    {
        Mon1d h;
        h.reserve(f.data.size() + g.data.size() - (size_t)index - 2);
        const Mon q = div(f.data[index], g.data[0]);
        size_t i = index + 1, j = 1;
        Mon prod;
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
            if (FnCmp::cmp(f.data[i], prod)) {
                updateProd = false;
                h.push_back(std::move(f.data[i++]));
            }
            else {
                updateProd = true;
                ++j;
                if (FnCmp::cmp(prod, f.data[i]))
                    h.push_back(std::move(prod));
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
    bool GCD(const Mon& mon1, const Mon& mon2, Mon& mon_out, const array& gen_degs, int& deg, int d_min, int d_max);

}  // namespace detail

template <typename FnCmp>
class Groebner
{
public:
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;

public:
    Poly1d data;

private:
    std::map<int, array> lc;

public:
    Groebner() = default;
    Groebner(Poly1d polys) : data(std::move(polys))
    {
        for (int i = 0; i < (int)data.size(); ++i)
            lc[data[i].GetLead().back().gen].push_back(i);
    }

    template <typename FnCmp1>
    Groebner<FnCmp1> ExtendMO() const /* Extend the monomial ordering */
    {
        Groebner<FnCmp1> gb1;
        for (const auto& p : data)
            gb1.push_back(Groebner<FnCmp1>::Poly{p.data});
        return gb1;
    }

public:
    Poly Reduce(Poly poly) const
    {
        size_t index = 0;
        while (index != poly.data.size()) {
            int gb_index = [&]() {
                for (const auto& p : poly.data[index]) {
                    if (lc.find(p.gen) != lc.end()) {
                        for (int i : lc.at(p.gen))
                            if (divisible(data[i].GetLead(), poly.data[index]))
                                return i;
                    }
                }
                return -1;
            }();
            if (gb_index != -1)
                alg::detail::Reduce(poly, data[gb_index], index);
            else
                ++index;
        }
        return poly;
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
        lc[g.GetLead().back().gen].push_back((int)data.size());
        data.push_back(std::move(g));
    }
};

struct GbBufferEle
{
    Mon gcd;
    int i1, i2;

    template <typename FnCmp>
    Polynomial<FnCmp> GetPoly(const Groebner<FnCmp>& gb)
    {
        return gb[i1] * div(gb[i2].GetLead(), gcd) + gb[i2] * div(gb[i1].GetLead(), gcd);
    }
};
using GroebnerLex = Groebner<CmpLex>;
using GroebnerRevlex = Groebner<CmpRevlex>;

using GbBuffer = std::map<int, std::vector<GbBufferEle>>;

/**
 * Comsume relations from `buffer` in degree `<= deg`
 * while adding new relations back to `buffer` in degree `<= deg_max`.
 * `deg=-1` or `deg_max=-1` means infinity.
 */
template <typename FnCmp, typename FnPred, typename FnDeg>
void AddRels(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, FnPred pred, FnDeg _get_deg, int deg, int deg_max)
{
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using PPoly1d = std::vector<const Poly*>;

    std::map<int, PPoly1d> rels_graded;
    for (const auto& rel : rels) {
        int d = _get_deg(rel.GetLead());
        if (d <= deg || deg == -1) {
            rels_graded[d].push_back(&rel);
            if (buffer.find(d) == buffer.end())
                buffer[d] = {};
        }
    }
    std::vector<std::pair<int, GbBufferEle>> buffer_new_rels;
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
                if (std::binary_search(rel.data.begin(), rel.data.end(), rel1.GetLead(), FnCmp::cmp))
                    rel += rel1;
            if (rel)
                rels_d.push_back(std::move(rel));
        }

        /* Add these relations */
        buffer_new_rels.resize(gb.size() + rels_d.size());
        for (auto& rel : rels_d) {
            if (rel) {
                ut::Range range2(0, (int)gb.size());
                std::for_each(std::execution::seq, range2.begin(), range2.end(), [&gb, &buffer_new_rels, &rel, pred, _get_deg, d, deg_max](int i) {
                    if (detail::HasGCD(rel.GetLead(), gb[i].GetLead()) && pred(rel.GetLead(), gb[i].GetLead())) {
                        Mon gcd = GCD(rel.GetLead(), gb[i].GetLead());
                        int deg_new_rel = d + _get_deg(gb[i].GetLead()) - _get_deg(gcd);
                        if (deg_new_rel <= deg_max || deg_max == -1)
                            buffer_new_rels[i] = std::make_pair(deg_new_rel, GbBufferEle{std::move(gcd), i, (int)gb.size()});
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
    }
    buffer.erase(buffer.begin(), p_buffer);
}

template <typename FnCmp>
void AddRels(Groebner<FnCmp>& gb, const std::vector<Polynomial<FnCmp>>& rels, const array& gen_degs, int deg)
{
    GbBuffer buffer;
    AddRels(
        gb, buffer, rels, [](const Mon&, const Mon&) { return true; }, [&gen_degs](const Mon& mon) { return GetDeg(mon, gen_degs); }, deg, deg);
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
Polynomial<FnCmp> subs(const Mon1d& data, FnMap map, const Groebner<FnCmp>& gb)
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
 * @see subs(const Poly&, FnType, const GbType&)
 */
template <typename FnCmp>
Polynomial<FnCmp> subs(const Mon1d& data, const std::vector<Polynomial<FnCmp>>& map, const Groebner<FnCmp>& gb)
{
    return subs(
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
                    buffer[deg_gcd].push_back(GbBufferEle{gcd, (int)i, (int)j});
                }
            }
        }
    }
    return buffer;
}

/**********************************************************
 * Algorithms that use Groebner basis
 **********************************************************/

constexpr unsigned int GEN_IDEAL = 0x80000000;
/**
 * Revlex on extra generators and FnCmp on the rest
 */
template <typename FnCmp>
struct CmpIdeal
{
    using submo = FnCmp;
    static constexpr std::string_view name = "CmpIdeal";
    static bool cmp(const Mon& m1, const Mon& m2)
    {
        auto mid1 = std::lower_bound(m1.begin(), m1.end(), GEN_IDEAL, [](const GenPow& p, unsigned int g) { return p.gen < g; });
        auto mid2 = std::lower_bound(m2.begin(), m2.end(), GEN_IDEAL, [](const GenPow& p, unsigned int g) { return p.gen < g; });
        if (CmpRevlex::cmp_ranges(mid1, m1.end(), mid2, m2.end()))
            return true;
        else if (std::equal(mid1, m1.end(), mid2, m2.end())) {
            if (FnCmp::cmp_ranges(m1.begin(), mid1, m2.begin(), mid2))
                return true;
        }
        return false;
    }
};

namespace detail {
    template <typename FnCmp>
    void AddRelsIdeal(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, const array& gen_degs, const array& gen_degs_y, int deg, int deg_max)
    {
        alg::AddRels(
            gb, buffer, rels, [](const Mon& m1, const Mon& m2) { return !((m1.back().gen & GEN_IDEAL) && (m2.back().gen & GEN_IDEAL) && (m1.back().gen != m2.back().gen)); },
            [&gen_degs, &gen_degs_y](const Mon& mon) {
                int result = 0;
                for (MonInd p = mon.begin(); p != mon.end(); ++p)
                    result += ((p->gen & GEN_IDEAL) ? gen_degs[p->gen] : gen_degs_y[p->gen - GEN_IDEAL]) * p->exp;
                return result;
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
std::vector<std::vector<alg::Polynomial<FnCmp>>>& Indecomposables(const Groebner<FnCmp>& gb, std::vector<std::vector<alg::Polynomial<FnCmp>>>& vectors, const array& gen_degs, const array& basis_degs)
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
    for (const auto& v : vectors) {
        PolyI rel;
        int d;
        for (size_t i = 0; i < basis_degs.size(); ++i)
            if (v[i])
                rel += PolyI::Sort((v[i] * Mon{{GEN_IDEAL + i, 1}}).data);
        for (size_t i = 0; i < basis_degs.size(); ++i) {
            if (v[i]) {
                degs.push_back(v[i].GetDeg(gen_degs) + basis_degs[i]);
                break;
            }
        }
        rels.push_back(std::move(rel));
    }
    array indices = ut::range((int)vectors.size());
    std::sort(indices.begin(), indices.end(), [&degs](int i, int j) { return degs[i] < degs[j]; });

    /* Add relations ordered by degree to gb1 */
    GbI gbI = gb.ExtendMO<CmpIdeal<FnCmp>>();
    GbBuffer buffer;
    int deg_max = degs[indices.back()];
    for (int i : indices) {
        PolyI rel = gbI.Reduce(rels[i]);
        if (!rel)
            vectors[i].clear();
        detail::AddRelsIdeal(gbI, buffer, {rel}, gen_degs, basis_degs, degs[i], deg_max);
    }

    /* Keep only the indecomposables in `vectors` */
    ut::RemoveEmptyElements(vectors);
    return vectors;
}

/**
 * Compute the generating set of linear relations among `polys`.
 *
 * The result is truncated by `deg<=deg_max`.
 */
template <typename FnCmp>
std::vector<std::vector<alg::Polynomial<FnCmp>>> ann_seq(const Groebner<FnCmp>& gb, const std::vector<alg::Polynomial<FnCmp>>& polys, const array& gen_degs, int deg_max)
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
                auto mid = std::lower_bound(m.begin(), m.end(), GEN_IDEAL, [](const GenPow& p, unsigned int g) { return p.gen < g; });
                Mon m1(m.begin(), mid), m2(mid, m.end());
                ann[m2.front().gen - GEN_IDEAL] += gb.Reduce(subs(
                                                                 {div(m2, {{m2.front().gen, 1}})}, [&polys](int i) { return polys[i - GEN_IDEAL]; }, gb)
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

//#define TEMPLATE_EXAMPLES
#ifdef TEMPLATE_EXAMPLES
namespace template_examples {
    using FnCmp = alg::CmpLex;
    using FnPred = bool (*)(alg::Mon, alg::Mon);
    using FnDeg = int (*)(alg::Mon);

    void AddRels_(Groebner<FnCmp>& gb, GbBuffer& buffer, const std::vector<Polynomial<FnCmp>>& rels, FnPred pred, FnDeg _get_deg, int deg, int deg_max)
    {
        return alg::AddRels(gb, buffer, rels, pred, _get_deg, deg, deg_max);
    }
}  // namespace template_examples
#endif

} /* namespace alg */

#endif /* GROEBNER_H */
