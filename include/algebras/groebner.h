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
#include <unordered_set>

/* extension of namespace alg in algebras.h */
namespace alg {

struct MMod
{
    uint32_t v;
    Mon m;
};

using MMod1d = std::vector<MMod>;

struct Mod
{
};

namespace detail {
    /*
     * Return if `mon1` and `mon2` have a nontrivial common factor.
     */
    bool HasGCD(const Mon& mon1, const Mon& mon2);
    inline bool HasGCD(const Mon& mon1, const Mon& mon2, MonTrace t1, MonTrace t2)
    {
        return (t1 & t2) && HasGCD(mon1, mon2);
    }

    int DegLCM(const Mon& mon1, const Mon& mon2, const int1d& gen_degs);
}  // namespace detail

class Groebner;

struct CriPair
{
    int i1 = -1, i2 = -1;
    Mon m1, m2;

    MonTrace trace_m2 = 0; /* = Trace(m2) */

    /* Compute the pair for two leading monomials. */
    CriPair() = default;
    static void SetFromLM(CriPair& result, const Mon& lead1, const Mon& lead2, int i, int j);

    /* Return `m1 * gb[i1] + m2 * gb[i2]` */
    void SijP(const Groebner& gb, Poly& result, Poly& tmp1, Poly& tmp2) const;
};
using CriPair1d = std::vector<CriPair>;
using CriPair2d = std::vector<CriPair1d>;

/* Groebner basis of critical pairs */
class GbCriPairs
{
private:
    int deg_trunc_;                                                      /* Truncation degree */
    CriPair2d gb_;                                                       /* `pairs_[j]` is the set of pairs (i, j) with given j */
    int2d min_pairs_;                                                    /* Minimal generating set of `pairs_` */
    std::map<int, CriPair2d> buffer_min_pairs_;                          /* To generate `buffer_min_pairs_` and for computing Sij */
    std::map<int, std::unordered_set<uint64_t>> buffer_redundent_pairs_; /* Used to minimize `buffer_min_pairs_` */

public:
    GbCriPairs(int d_trunc) : deg_trunc_(d_trunc) {}
    int deg_trunc() const
    {
        return deg_trunc_;
    }
    CriPair1d Criticals(int d)
    {
        CriPair1d result;
        if (!buffer_min_pairs_.empty() && buffer_min_pairs_.begin()->first == d) {
            auto& b_min_pairs_d = buffer_min_pairs_.begin()->second;
            for (size_t j = 0; j < b_min_pairs_d.size(); ++j)
                for (auto& pair : b_min_pairs_d[j])
                    if (pair.i2 != -1)
                        result.push_back(std::move(pair));
            buffer_min_pairs_.erase(buffer_min_pairs_.begin());
        }
        return result;
    }
    bool empty() const
    {
        return buffer_min_pairs_.empty();
    }

    /* Minimize `buffer_min_pairs_[d]` and maintain `pairs_` */
    void Minimize(const Mon1d& leads, int d);

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const Mon1d& leads, const MonTrace1d& traces, const Mon& mon, const int1d& gen_degs);

    void init(const Mon1d& leads, const MonTrace1d& traces, const int1d& gen_degs, int t_min_buffer);
};

class Groebner  // TODO: add gen_degs
{
private:
    using TypeIndexKey = uint32_t;
    using TypeIndices = std::unordered_map<TypeIndexKey, int1d>;

private:
    GbCriPairs criticals_; /* Groebner basis of critical pairs */

    Poly1d data_;
    int1d data_degs_;     /* Degrees of data_. */
    Mon1d leads_;         /* Leading monomials */
    MonTrace1d traces_;   /* Cache for fast divisibility test */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int1d gen_degs_; /* degree of generators */

public:
    Groebner(int deg_trunc, int1d gen_degs) : criticals_(deg_trunc), gen_degs_(std::move(gen_degs)) {}

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    Groebner(int deg_trunc, int1d gen_degs, Poly1d polys, bool bDynamic = false);

private:
    static TypeIndexKey Key(const Mon& lead)
    {
        return TypeIndexKey{lead.back().g() + (lead.size() == 1 ? 0 : ((lead[lead.size() - 2].g() + 1) << 16))};
    }

public: /* Getters and Setters */
    const GbCriPairs& gb_pairs() const
    {
        return criticals_;
    }
    int deg_trunc() const
    {
        return criticals_.deg_trunc();
    }
    /* This function will erase `gb_pairs_.buffer_min_pairs[d]` */
    CriPair1d Criticals(int d)
    {
        return criticals_.Criticals(d);
    }

    auto size() const
    {
        return data_.size();
    }
    const Poly& operator[](size_t index) const
    {
        return data_[index];
    }
    void push_back(Poly g)
    {
        const Mon& m = g.GetLead();
        criticals_.AddToBuffers(leads_, traces_, m, gen_degs_);

        data_degs_.push_back(GetDeg(m, gen_degs_));
        leads_.push_back(m);
        traces_.push_back(m.Trace());
        indices_[Key(m)].push_back((int)data_.size());
        data_.push_back(std::move(g));
    }
    bool operator==(const Groebner& rhs) const
    {
        return data_ == rhs.data_;
    }

    const auto& data() const
    {
        return data_;
    }

    /* Return trace of LM(gb[i]) */
    MonTrace GetTrace(size_t i) const
    {
        return traces_[i];
    }

public:
    /* Leadings[i] is the monomials that end with generator i.
     * The result is used to generate basis of $P/I$.
     */
    Mon2d GetLeadings(size_t gens_size) const
    {
        Mon2d result;
        result.resize(gens_size);
        for (size_t i = 0; i < data_.size(); ++i)
            result[leads_[i].back().g()].push_back(leads_[i]);
        return result;
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(const Mon& mon) const;

    Poly Reduce(Poly poly) const;

    /**
     * Comsume relations from 'rels` and `gb.criticals_` in degree `<= deg`
     */
    void AddRels(const Poly1d& rels, int deg);
};

/**
 * A fast algorithm that computes
 * `poly ** n` modulo `gb`
 */
void powP(const Poly& poly, int n, const Groebner& gb, Poly& result, Poly& tmp);
inline Poly pow(const Poly& poly, int n, const Groebner& gb)
{
    Poly result, tmp;
    powP(poly, n, gb, result, tmp);
    return result;
}

/**
 * Replace the generators in `poly` with elements given in `map`
 * and evaluate modulo `gb`.
 * @param poly The polynomial to be substituted.
 * @param map `map(i)` is the polynomial that substitutes the generator of id `i`.
 * @param gb A Groebner basis.
 */
template <typename FnMap>
Poly subsMGbTpl(const Poly& poly, const Groebner& gb, const FnMap& map)
{
    Poly result, tmp_prod, tmp;
    for (const Mon& m : poly.data) {
        Poly fm = Poly::Unit();
        for (auto p = m.begin(); p != m.end(); ++p) {
            powP(map(p->g()), p->e(), gb, tmp_prod, tmp);
            fm.imulP(tmp_prod, tmp);
        }
        result.iaddP(fm, tmp);
    }
    return result;
}

inline Poly subsMGb(const Poly& poly, const Groebner& gb, const std::vector<Poly>& map)
{
    return subsMGbTpl(poly, gb, [&map](size_t i) { return map[i]; });
}

} /* namespace alg */

#endif /* GROEBNER_H */
