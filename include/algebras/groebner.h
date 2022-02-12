/** \file groebner.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
#include "linalg.h"
#include "myio.h"  ////
#include <execution>
#include <iostream>  ////
#include <map>
#include <memory>
#include <unordered_map>

/* extension of namespace alg in algebras.h */
namespace alg {

namespace detail {
    void mul(const Mon& mon1, const Mon& mon2, MonOnStack& result);
    void mul(const Mon& mon1, const MonOnStack& mon2, MonOnStack& result);
    void div(const Mon& mon1, const Mon& mon2, MonOnStack& result);
    void div(const Mon& mon1, const MonOnStack& mon2, MonOnStack& result);

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
    inline bool HasGCD(const Mon& mon1, const Mon& mon2, MonTrace t1, MonTrace t2)
    {
        return (t1 & t2) && HasGCD(mon1, mon2);
    }

    template <typename FnGenDeg>
    int DegLCM(const Mon& mon1, const Mon& mon2, FnGenDeg _gen_deg)
    {
        int result = 0;
        MonInd k = mon1.begin(), l = mon2.begin();
        while (k != mon1.end() && l != mon2.end()) {
            if (k->gen < l->gen) {
                result += _gen_deg(k->gen) * k->exp;
                ++k;
            }
            else if (k->gen > l->gen) {
                result += _gen_deg(l->gen) * l->exp;
                ++l;
            }
            else {
                if (k->exp < l->exp)
                    result += _gen_deg(l->gen) * l->exp;
                else
                    result += _gen_deg(k->gen) * k->exp;
                ++k;
                ++l;
            }
        }
        for (; k != mon1.end(); ++k)
            result += _gen_deg(k->gen) * k->exp;
        for (; l != mon2.end(); ++l)
            result += _gen_deg(l->gen) * l->exp;
        return result;
    }

}  // namespace detail

template <typename FnCmp>
class Groebner;

struct CriticalPair
{
    int i1, i2 = -1;
    Mon m1, m2;
    MonTrace trace_m2; /* = Trace(m2) */
    bool trivialSyz;   /* Trivial Syzygy */

    /* Compute the pair for two leading monomials. */
    CriticalPair() = default;
    static void SetFromLM(CriticalPair& result, const Mon& lead1, const Mon& lead2, int i, int j, bool trivial)
    {
        if (trivial) {
            result.m1 = lead2;
            result.m2 = lead1;
        }
        else {
            MonInd k = lead1.begin(), l = lead2.begin();
            while (k != lead1.end() && l != lead2.end()) {
                if (k->gen < l->gen)
                    result.m2.push_back(*k++);
                else if (k->gen > l->gen)
                    result.m1.push_back(*l++);
                else {
                    if (k->exp < l->exp)
                        result.m1.emplace_back(k->gen, l->exp - k->exp);
                    else if (k->exp > l->exp)
                        result.m2.emplace_back(k->gen, k->exp - l->exp);
                    k++;
                    l++;
                }
            }
            if (k != lead1.end())
                result.m2.insert(result.m2.end(), k, lead1.end());
            else
                result.m1.insert(result.m1.end(), l, lead2.end());
        }
        result.i1 = i;
        result.i2 = j;
        result.trace_m2 = Trace(result.m2);
        result.trivialSyz = trivial;
    }

    /* Return `m1 * gb[i1] + m2 * gb[i2]` */
    template <typename FnCmp>
    Polynomial<FnCmp> Sij(const Groebner<FnCmp>& gb) const
    {
        Polynomial<FnCmp> result;
        detail::MonOnStack prod1, prod2;
        result.data.reserve(gb[i1].data.size() + gb[i2].data.size() - 2);
        auto p1 = gb[i1].data.begin() + 1, p2 = gb[i2].data.begin() + 1;
        while (p1 != gb[i1].data.end() && p2 != gb[i2].data.end()) {
            mul(*p1, m1, prod1);
            mul(*p2, m2, prod2);
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
            mul(*p1++, m1, prod1);
            result.data.push_back(Mon(prod1));
        }
        while (p2 != gb[i2].data.end()) {
            mul(*p2++, m2, prod2);
            result.data.push_back(Mon(prod2));
        }
        return result;
    }
};
using CriticalPair1d = std::vector<CriticalPair>;
using CriticalPair2d = std::vector<CriticalPair1d>;
using BufferCriticalPairs = std::map<int, CriticalPair1d>;

/* Groebner basis of critical pairs */
class GbCriPairs
{
public:
    static constexpr uint32_t MODE_ON = 0x0001; /* gb_pairs_ is in use */
    static constexpr uint32_t MODE_GB = 0x0002; /* Groebner basis of critical pairs is recorded */

private:
    int deg_trunc_;                                   /* Truncation degree */
    CriticalPair2d pairs_;                            /* `pairs_[j]` is the set of pairs (i, j) with given j */
    std::map<int, pairint1d> buffer_redundent_pairs_; /* To be removed from generating set */
    BufferCriticalPairs buffer_min_pairs_;            /* To generate the minimal generating set */
    uint32_t mode_;                                   /* This determines how much information will be recorded during the computation */

public:
    GbCriPairs(int d_trunc, int mode) : deg_trunc_(d_trunc), mode_(mode) {}
    int deg_trunc() const
    {
        return deg_trunc_;
    }
    void set_deg_trunc(int d_trunc)
    {
        if (d_trunc <= deg_trunc_) {
            deg_trunc_ = d_trunc;  ////
        }
        else {
            throw MyException(0x540f459dU, "NotImplemented");  ////
        }
    }
    uint32_t mode() const
    {
        return mode_;
    }
    void set_mode(uint32_t mode)
    {
        mode_ = mode;
    }
    const CriticalPair1d& min_pairs(int d) const
    {
        return buffer_min_pairs_.at(d);
    }
    size_t size_min_pairs(int d) const
    {
        auto p = buffer_min_pairs_.find(d);
        return p != buffer_min_pairs_.end() ? p->second.size() : 0;
    }
    bool empty_min_pairs() const
    {
        return buffer_min_pairs_.empty();
    }
    void erase_min_pairs(int d)
    {
        buffer_min_pairs_.erase(d);
    }

    /* Minimize `buffer_min_pairs_[d]` and maintain `pairs_` */
    void AddAndMinimize(const Mon1d& leads, int d)
    {
        /* Add to the Groebner basis of critical pairs */
        if (mode_ & MODE_GB) {
            for (const auto& pair : buffer_min_pairs_[d]) {
                if (pairs_.size() <= pair.i2)
                    pairs_.resize(size_t(pair.i2 + 1));
                pairs_[pair.i2].push_back(pair);
            }
        }

        /* Minimize `buffer_min_pairs_` */
        if (mode_ & MODE_GB) {
            for (const auto [i, j] : buffer_redundent_pairs_[d]) {
                int k = i, l = j;
                while (l != -1) {
                    Mon gcd = GCD(leads[k], leads[l]);
                    Mon m2 = div(leads[k], gcd);
                    auto p = std::find_if(buffer_min_pairs_[d].begin(), buffer_min_pairs_[d].end(), [&m2](const CriticalPair& c) { return c.m2 == m2; });
                    if (p != buffer_min_pairs_[d].end()) { /* Remove it from `buffer_min_pairs_` */
                        p->i2 = -1;
                        l = -1;
                    }
                    else { /* Reduce (k, l) */
                        MonTrace t_m2 = Trace(m2);
#ifndef NDEBUG
                        int old_l = l;
#endif
                        for (const CriticalPair& c : pairs_[l]) {
                            if (divisible(c.m2, m2)) {
                                Mon m1 = div(leads[l], gcd);
                                if (detail::HasGCD(c.m1, m1)) {
                                    l = -1;
                                    break;
                                }
                                else {
                                    l = c.i1;
                                    if (k > l)
                                        std::swap(k, l);
                                    break;
                                }
                            }
                        }
#ifndef NDEBUG
                        if (l == old_l)
                            throw MyException(0x6ddee1b3U, "BUG: l is expected to change.");
#endif
                    }
                }
            }
        }
        else {
            for (const pairint ij : buffer_redundent_pairs_[d]) {
                auto p = std::find_if(buffer_min_pairs_[d].begin(), buffer_min_pairs_[d].end(), [ij](const CriticalPair& c) { return c.i1 == ij.first && c.i2 == ij.second; });
                if (p != buffer_min_pairs_[d].end()) /* Remove it from `buffer_min_pairs_` */
                    p->i2 = -1;
            }
        }
        ut::RemoveIf(buffer_min_pairs_[d], [](const CriticalPair& c) { return c.i2 == -1; });

        /* Delete `bufbuffer_redundent_pairs_[d]` */
        buffer_redundent_pairs_.erase(d);
    }

    /* Remove theses pairs (i, j) from `buffer_min_pairs_` for which Sij = 0 */
    void RemoveTrivialSyzygies(int d)
    {
        ut::RemoveIf(buffer_min_pairs_[d], [](const CriticalPair& c) { return c.trivialSyz; });
    }

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    template <typename FnPred, typename FnDeg>
    void AddToBuffers(const Mon1d& leads, const array& d_leads, const MonTrace1d& traces, const Mon& mon, int d_mon, FnPred pred, FnDeg _gen_deg)
    {
        MonTrace t = Trace(mon);
        std::vector<std::pair<int, CriticalPair>> new_pairs(leads.size());
        ut::Range range(0, (int)leads.size());

        for (size_t i = 0; i < leads.size(); ++i) {
            if (pred(leads[i], mon)) {
                int d_pair = detail::DegLCM(leads[i], mon, _gen_deg);
                if (d_pair <= deg_trunc_) {
                    new_pairs[i].first = d_pair;
                    CriticalPair::SetFromLM(new_pairs[i].second, leads[i], mon, (int)i, (int)leads.size(), d_leads[i] + d_mon == d_pair);
                }
            }
            else
                new_pairs[i].first = -1;
        }

        /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
        for (size_t j = 1; j < new_pairs.size(); ++j) {
            if (new_pairs[j].first != -1) {
                for (size_t i = 0; i < j; ++i) {
                    if (new_pairs[i].first != -1) {
                        if (divisible(new_pairs[i].second.m2, new_pairs[j].second.m2, new_pairs[i].second.trace_m2, new_pairs[j].second.trace_m2)) {
                            new_pairs[j].first = -1;
                            break;
                        }
                        else if (divisible(new_pairs[j].second.m2, new_pairs[i].second.m2, new_pairs[j].second.trace_m2, new_pairs[i].second.trace_m2))
                            new_pairs[i].first = -1;
                        else if (!detail::HasGCD(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                            int dij = detail::DegLCM(leads[i], leads[j], _gen_deg);
                            if (dij <= deg_trunc_) {
                                auto redundent_pair = std::make_pair((int)i, (int)j);
                                if (std::find(buffer_redundent_pairs_[dij].begin(), buffer_redundent_pairs_[dij].end(), redundent_pair) == buffer_redundent_pairs_[dij].end())
                                    buffer_redundent_pairs_[dij].push_back(redundent_pair);
                            }
                        }
                    }
                }
            }
        }
        for (size_t i = 0; i < new_pairs.size(); ++i) {
            if (new_pairs[i].first != -1)
                buffer_min_pairs_[new_pairs[i].first].push_back(std::move(new_pairs[i].second));
        }
    }
};

template <typename FnCmp>
class Groebner  // TODO: add gen_degs
{
private:
    using TypeIndexKey = int;
    using TypeIndex = std::unordered_map<TypeIndexKey, alg::array>;

public:
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;

private:
    Poly1d data_;
    GbCriPairs gb_pairs_; /* Groebner basis of critical pairs */

    /* Caches */
    array data_degs_;   /* Degrees of data_. */
    Mon1d leads_;       /* Leading monomials */
    MonTrace1d traces_; /* Cache for fast divisibility test */
    TypeIndex index_;   /* Cache for fast divisibility test */

public:
    Groebner(int deg_trunc, uint32_t mode = GbCriPairs::MODE_ON | GbCriPairs::MODE_GB) : gb_pairs_(deg_trunc, mode) {}

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    Groebner(int deg_trunc, Poly1d polys) : data_(std::move(polys)), gb_pairs_(deg_trunc, 0)  // TODO: Cache degrees
    {
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            traces_.push_back(Trace(data_[i].GetLead()));
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

    /* Warning: This does not change `gb_pairs_` */
    void InitFrom(Poly1d polys)  // TODO: Cache degrees
    {
        data_ = std::move(polys);
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            traces_.push_back(Trace(data_[i].GetLead()));
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

    /* Initialize from `polys` and `gb_pairs` where `polys` already forms a Groebner basis. */
    Groebner(Poly1d polys, GbCriPairs gb_pairs) : data_(std::move(polys)), gb_pairs_(std::move(gb_pairs))  // TODO: Cache degrees
    {
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            traces_.push_back(Trace(data_[i].GetLead()));
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

    /* Initialize from `polys` and `gb_pairs` where `polys` already forms a Groebner basis. */
    Groebner(Poly1d polys, array polys_degs, GbCriPairs gb_pairs) : data_(std::move(polys)), gb_pairs_(std::move(gb_pairs)), data_degs_(std::move(polys_degs))  // TODO: Cache degrees
    {
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            traces_.push_back(Trace(data_[i].GetLead()));
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

    /* Extend the monomial ordering */
    template <typename FnCmp1>
    Groebner<FnCmp1> ExtendMO() const
    {
        Polynomial1d<FnCmp1> data1;
        for (const auto& p : data_)
            data1.push_back(Polynomial<FnCmp1>{p.data});
        return Groebner<FnCmp1>(std::move(data1), data_degs_, gb_pairs_);
    }

private:
    static TypeIndexKey Key(const Mon& lead)
    {
        return TypeIndexKey{lead.size() == 1 ? lead.back().gen : lead.back().gen + (lead[lead.size() - 2].gen << 16)};
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(const Mon& mon) const
    {
        auto t = Trace(mon);
        for (int i = 0; i < (int)mon.size(); ++i) {
            for (int j = -1; j < i; ++j) {
                auto key = TypeIndexKey{j == -1 ? mon[i].gen : mon[i].gen + (mon[j].gen << 16)};
                auto p = index_.find(key);
                if (p != index_.end()) {
                    for (int k : p->second) {
                        if (divisible(leads_[k], mon, traces_[k], t))
                            return k;
                    }
                }
            }
        }
        return -1;
    }

public: /* Getters and Setters */
    const GbCriPairs& gb_pairs() const
    {
        return gb_pairs_;
    }
    int deg_trunc() const
    {
        return gb_pairs_.deg_trunc();
    }
    void set_deg_trunc(int d_trunc)
    {
        gb_pairs_.set_deg_trunc(d_trunc);
    }
    uint32_t mode() const
    {
        return gb_pairs_.mode();
    }
    void set_mode(uint32_t mode)
    {
        gb_pairs_.set_mode(mode);
    }
    void erase_min_pairs(int d)
    {
        gb_pairs_.erase_min_pairs(d);
    }

    auto size() const
    {
        return data_.size();
    }
    auto& operator[](size_t index) const
    {
        return data_[index];
    }
    template <typename FnPred, typename FnDeg>
    void push_back(Poly g, int deg, FnPred pred, FnDeg _gen_deg)
    {
        const Mon& m = g.GetLead();
        gb_pairs_.AddToBuffers(leads_, data_degs_, traces_, m, deg, pred, _gen_deg);

        data_degs_.push_back(deg);
        leads_.push_back(m);
        traces_.push_back(Trace(m));
        index_[Key(m)].push_back((int)data_.size());
        data_.push_back(std::move(g));
    }
    void AddPairsAndMinimize(int deg)
    {
        gb_pairs_.AddAndMinimize(leads_, deg);
    }

    void RemoveTrivialSyzygies(int deg)
    {
        gb_pairs_.RemoveTrivialSyzygies(deg);
    }
    bool operator==(const Groebner<FnCmp>& rhs) const
    {
        return data_ == rhs.data_;
    }

    const Poly1d& data() const
    {
        return data_;
    }

    /* Return trace of LM(gb[i]) */
    MonTrace trace_LM(size_t i) const
    {
        return traces_[i];
    }

public:
    /* Leadings[i] is the monomials that end with generator i.
     * The result is used to generate basis of $P/I$.
     */
    Mon2d GetLeadings(size_t gens_size) const
    {
        alg::Mon2d result;
        result.resize(gens_size);
        for (size_t i = 0; i < data_.size(); ++i)
            result[leads_[i].back().gen].push_back(leads_[i]);
        return result;
    }

    Poly Reduce(Poly poly) const
    {
        size_t index = 0;
        while (index < poly.data.size()) {
            int gb_index = IndexOfDivisibleLeading(poly.data[index]);
            if (gb_index != -1) {
                alg::detail::Reduce(poly, data_[gb_index], index);
            }
            else
                ++index;
        }
        return poly;
    }
};

using GroebnerLex = Groebner<CmpLex>;
using GroebnerRevlex = Groebner<CmpRevlex>;

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
Polynomial<FnCmp> SubsMGb(const Mon1d& data, const Groebner<FnCmp>& gb, FnMap map)
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
 * @see SubsMGb(const Poly&, const GbType&, FnType)
 */
template <typename FnCmp>
Polynomial<FnCmp> SubsMGb(const Mon1d& data, const Groebner<FnCmp>& gb, const std::vector<Polynomial<FnCmp>>& map)
{
    return SubsMGb(data, gb, [&map](int i) { return map[i]; });
}

/**
 * Comsume relations from 'rels` and `gb.gb_pairs_` in degree `<= deg`
 * while adding new relations back to `buffer` in degree `<= deg_max`.
 */
template <typename FnCmp, typename FnPred, typename FnDeg>
void TplAddRels(Groebner<FnCmp>& gb, const Polynomial1d<FnCmp>& rels, int deg, FnPred pred, FnDeg _gen_deg)
{
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using PPoly1d = std::vector<const Poly*>;

    if (!(gb.mode() & GbCriPairs::MODE_ON))
        throw MyException(0x49fadb32U, "gb is in const mode.");
    int deg_max = gb.deg_trunc();
    if (deg > deg_max)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, PPoly1d> rels_graded;
    for (const auto& rel : rels) {
        if (rel) {
            int d = TplGetDeg(rel.GetLead(), _gen_deg);
            if (d <= deg)
                rels_graded[d].push_back(&rel);
        }
    }
    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg && (d <= deg_max_rels || !gb.gb_pairs().empty_min_pairs()); ++d) {
        // std::cout << "d=" << d << '\n';  ////
        int size_min_pairs_d = (int)gb.gb_pairs().size_min_pairs(d);
        if (size_min_pairs_d) {
            gb.AddPairsAndMinimize(d);
            gb.RemoveTrivialSyzygies(d);
            size_min_pairs_d = (int)gb.gb_pairs().size_min_pairs(d);
            if (!size_min_pairs_d)
                gb.erase_min_pairs(d);
        }
        auto p_rels_d = rels_graded.find(d);
        int size_rels_d = p_rels_d != rels_graded.end() ? (int)p_rels_d->second.size() : 0;
        if (size_min_pairs_d + size_rels_d == 0)
            continue;
        Poly1d rels_tmp(size_min_pairs_d + size_rels_d);

        /* Reduce relations in degree d */
        if (size_min_pairs_d) {
            auto& min_pairs = gb.gb_pairs().min_pairs(d);
            ut::Range range(0, size_min_pairs_d);
            std::for_each(std::execution::par_unseq, range.begin(), range.end(), [&gb, &min_pairs, &rels_tmp](int i) { rels_tmp[i] = gb.Reduce(min_pairs[i].Sij(gb)); });
            gb.erase_min_pairs(d);
        }
        if (size_rels_d) {
            auto& rels_d = p_rels_d->second;
            ut::Range range(0, size_rels_d);
            std::for_each(range.begin(), range.end(), [&gb, &rels_d, &rels_tmp, size_min_pairs_d](int i) { rels_tmp[size_min_pairs_d + i] = gb.Reduce(*rels_d[i]); });
        }

        /* Triangulate these relations */
        Poly1d rels_d;
        for (auto& rel : rels_tmp) {
            for (const Poly& rel1 : rels_d)
                if (std::binary_search(rel.data.begin(), rel.data.end(), rel1.GetLead(), FnCmp::template cmp<Mon, Mon>))
                    rel += rel1;
            if (rel)
                rels_d.push_back(std::move(rel));
        }

        /* Add these relations */
        for (auto& rel : rels_d)
            gb.push_back(std::move(rel), d, pred, _gen_deg);
    }
}

template <typename FnCmp>
void AddRels(Groebner<FnCmp>& gb, const Polynomial1d<FnCmp>& rels, int deg, const array& gen_degs)
{
    TplAddRels(
        gb, rels, deg, [](const Mon&, const Mon&) { return true; }, [&gen_degs](int i) { return gen_degs[i]; });
}

/**
 * Generate buffer in degree `t_min <= t <= t_max`
 */
// template <typename FnCmp>
// BufferCriticalPairs GenerateBuffer(const alg::Groebner<FnCmp>& gb, const array& gen_degs, int d_min, int deg_trunc)
//{
//    BufferCriticalPairs buffer;
//    Mon gcd;
//    int deg_gcd;
//    for (size_t i = 0; i < gb.data.size(); ++i) {
//        int deg1 = gb[i].GetDeg(gen_degs);
//        for (size_t j = i + 1; j < gb.data.size(); ++j) {
//            int deg2 = gb[j].GetDeg(gen_degs);
//            if (d_min < deg1 + deg2 && std::max(deg1, deg2) <= deg_trunc) {
//                if (detail::GCD(gb[i].GetLead(), gb[j].GetLead(), gcd, gen_degs, deg_gcd, std::max(deg1 + deg2 - deg_trunc, 1), deg1 + deg2 - d_min)) {
//                    buffer[deg_gcd].push_back(CriticalPair{gcd, (int)i, (int)j});
//                }
//            }
//        }
//    }
//    return buffer;
//}

/**********************************************************
 * Algorithms that use Groebner basis
 **********************************************************/

inline constexpr int GEN_IDEAL = 0x4000000;
inline constexpr int GEN_SUB = 0x2000000;

namespace detail {
    /* For monomials in the homological groebner basis get the degree of the first part consisting of original generators and set `p` to be
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
}  // namespace detail

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
    void AddRelsIdeal(Groebner<FnCmp>& gb, const std::vector<Polynomial<FnCmp>>& rels, int deg, const array& gen_degs, const array& gen_degs_y)
    {
        TplAddRels(
            gb, rels, deg, [](const Mon& m1, const Mon& m2) { return !(m1.back().gen & GEN_IDEAL) && !(m2.back().gen & GEN_IDEAL); }, [&gen_degs, &gen_degs_y](int i) { return ((i & GEN_IDEAL) ? gen_degs_y[i - GEN_IDEAL] : gen_degs[i]); });
    }
    template <typename FnCmp>
    void AddRelsModule(Groebner<FnCmp>& gb, const std::vector<Polynomial<FnCmp>>& rels, int deg, const array& gen_degs, const array& gen_degs_y)
    {
        TplAddRels(
            gb, rels, deg, [](const Mon& m1, const Mon& m2) { return !((m1.back().gen & GEN_IDEAL) && (m2.back().gen & GEN_IDEAL) && (m1.back().gen != m2.back().gen)); },
            [&gen_degs, &gen_degs_y](int i) { return ((i & GEN_IDEAL) ? gen_degs_y[i - GEN_IDEAL] : gen_degs[i]); });
    }
    template <typename FnCmp>
    void AddRelsHomology(Groebner<FnCmp>& gb, const std::vector<Polynomial<FnCmp>>& rels, int deg, const MayDeg1d& gen_degs, const MayDeg1d& gen_degs_x, const MayDeg1d& gen_degs_b)
    {
        TplAddRels(
            gb, rels, deg, [](const Mon& m1, const Mon& m2) { return true; },
            [&gen_degs, &gen_degs_x, &gen_degs_b](int i) {
                if (i & GEN_IDEAL)
                    return gen_degs_b[i - GEN_IDEAL - GEN_SUB].t;
                else if (i & GEN_SUB)
                    return gen_degs_x[i - GEN_SUB].t;
                else
                    return gen_degs[i].t;
            });
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
        int deg_max = degs[indices.back()];
        GbI gbI = gb.ExtendMO<CmpIdeal<FnCmp>>();
        gbI.set_deg_trunc(deg_max);
        gbI.set_mode(GbCriPairs::MODE_ON);
        for (int i : indices) {
            detail::AddRelsModule(gbI, {}, degs[i], gen_degs, basis_degs);
            PolyI rel = gbI.Reduce(rels[i]);
            if (rel)
                detail::AddRelsModule(gbI, {std::move(rel)}, degs[i], gen_degs, basis_degs);
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
    gbI.set_deg_trunc(deg_max);
    gbI.set_mode(GbCriPairs::MODE_ON);
    detail::AddRelsIdeal(gbI, rels, deg_max, gen_degs, gen_degs_y);

    /* Extract linear relations from gb1 */
    for (const PolyI& g : gbI.data()) {
        if (g.GetLead().back().gen & GEN_IDEAL) {
            Poly1d ann;
            ann.resize(n);
            for (const Mon& m : g.data) {
                auto mid = std::lower_bound(m.begin(), m.end(), GEN_IDEAL, [](const GenPow& p, int g) { return p.gen < g; });
                Mon m1(m.begin(), mid), m2(mid, m.end());  // TODO: improve
                ann[m2.front().gen - GEN_IDEAL] += gb.Reduce(SubsMGb({div(m2, GenPow::Mon(m2.front().gen))}, gb, [&polys](int i) { return polys[i - GEN_IDEAL]; }) * m1);
            }
            result.push_back(std::move(ann));
        }
    }

    /* Add commutators to linear relations */
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (gen_degs_y[i] + gen_degs_y[j] <= deg_max) {
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

/* Propergate `basis` to degree t=`t` */
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
              std::vector<std::string>& gen_names_h, std::vector<Polynomial<FnCmp>>& gen_repr_h, int t_trunc)
{
    using Poly = Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using PolyH = Polynomial<CmpHomology<FnCmp>>;
    using PolyH1d = std::vector<PolyH>;
    using Gb = Groebner<FnCmp>;
    using GbH = Groebner<CmpHomology<FnCmp>>;

    if (gb.deg_trunc() < t_trunc)
        throw MyException(0x12ddc7ccU, "The truncation degree of the homology should be smaller");

    Gb gb_AoAZ = gb;
    gb_AoAZ.set_deg_trunc(t_trunc);
    gb_AoAZ.set_mode(GbCriPairs::MODE_ON);
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
    gb_H.set_deg_trunc(t_trunc);
    gb_H.set_mode(GbCriPairs::MODE_ON);
    MayDeg1d gen_degs_b;
    Poly1d gen_reprs_b; /* dm_i = b_i */

    int deg_t_allX = 0;
    for (auto& deg : gen_degs)
        deg_t_allX += deg.t;
    // std::cout << "deg_t_allX=" << deg_t_allX << '\n';

    for (int t = 1; t <= t_trunc; ++t) {
        // std::cout << "t=" << t << '\n';  //
        PolyH1d rels_h_t;
        Poly1d rels_AoAZ_t;
        if (t <= deg_t_allX) {
            /* Find m_i, b_i and x_i */
            for (size_t i = gb_AoAZ_old_size; i < gb_AoAZ.data().size(); ++i)
                leadings_AoAZ[gb_AoAZ[i].GetLead().back().gen].push_back(gb_AoAZ[i].GetLead());  // TODO: leads_
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
            AddRels(gb_AoAZ, rels_AoAZ_t, t, CmpHomology<FnCmp>::gen_degs);
            rels_AoAZ_t.clear();
        }

        /* Find  y_i */
        size_t i_start_t = gb_H.size();
        detail::AddRelsHomology(gb_H, rels_h_t, t, gen_degs, gen_degs_h, gen_degs_b);
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
                    Poly repr_y = SubsMGb(y.data, gb, [&gen_repr_h, &diff_basis_AoAZ_flat](int i) { return ((i & GEN_IDEAL) ? diff_basis_AoAZ_flat[i - GEN_IDEAL - GEN_SUB] : ((i & GEN_SUB) ? gen_repr_h[i - GEN_SUB] : Poly::Gen(i))); });
                    detail::AddRelsHomology(gb_H, {y + PolyH::Gen(GEN_SUB + (int)gen_degs_h.size())}, t, gen_degs, gen_degs_h, gen_degs_b);
                    gen_names_h.push_back("y_{" + std::to_string(gen_degs_h.size()) + "}");
                    gen_degs_h.push_back(repr_y.GetMayDeg(gen_degs));
                    gen_repr_h.push_back(std::move(repr_y));
                    rels_AoAZ_t.push_back(gen_repr_h.back());
                }
            }
        }
        AddRels(gb_AoAZ, rels_AoAZ_t, t, CmpHomology<FnCmp>::gen_degs);
        rels_AoAZ_t.clear();
    }
    /* Prepare relations for homology */
    PolyRevlex1d polys_h;
    for (auto& p : gb_H.data()) {
        if (p && (p.GetLead().front().gen & GEN_SUB)) {
            auto p_m = std::lower_bound(p.data.begin(), p.data.end(), 0, [](const Mon& m, int) { return !(m.back().gen & GEN_IDEAL); });
            if (p_m != p.data.begin()) {
                Mon1d p1(p.data.begin(), p_m);
                polys_h.push_back(alg::subs<alg::CmpRevlex>(p1, [](int i) { return PolyRevlex::Gen(i - GEN_SUB); }));
            }
        }
    }
    gb_h.InitFrom(polys_h);
}

#define TEMPLATE_EXAMPLES
#ifdef TEMPLATE_EXAMPLES
namespace template_examples {
    using FnCmp = alg::CmpLex;
    using FnPred = bool (*)(alg::Mon, alg::Mon);
    using FnDeg = int (*)(int);

    inline void TplAddRels_(Groebner<FnCmp>& gb, const Polynomial1d<FnCmp>& rels, int deg, FnPred pred, FnDeg _gen_deg)
    {
        return alg::TplAddRels(gb, rels, deg, pred, _gen_deg);
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
