/** \file groebner.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
#include <map>
#include <unordered_set>

/* extension of namespace alg in algebras.h */
namespace alg2 {

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

    void MutualQuotient(Mon& m1, Mon& m2, const Mon& lead1, const Mon& lead2);
}  // namespace detail

class Groebner;
class GroebnerMod;

inline constexpr uint32_t NULL_INDEX32 = ~uint32_t(0);
inline constexpr uint32_t FLAG_INDEX_X = uint32_t(1) << 31;

struct CriPair
{
    uint32_t i1 = NULL_INDEX32, i2 = NULL_INDEX32;
    Mon m1, m2;

    MonTrace trace_m2 = 0; /* = Trace(m2) */

    /* Compute the pair for two leading monomials. */
    CriPair() = default;
    static void SetFromLM(CriPair& result, const Mon& lead1, const Mon& lead2, int i, int j);

    /* Return `m1 * gb[i1] + m2 * gb[i2]` */
    void SijP(const Groebner& gb, Poly& result, Poly& tmp1, Poly& tmp2) const;
    void SijMod(const Groebner& gb, const GroebnerMod& gbm, Mod& result, Mod& tmp1, Mod& tmp2) const;
};
using CriPair1d = std::vector<CriPair>;
using CriPair2d = std::vector<CriPair1d>;

/* Groebner basis of critical pairs */
class GbCriPairs
{
private:
    int deg_trunc_;                                                      /* Truncation degree */
    CriPair2d gb_;                                                       /* `pairs_[j]` is the set of pairs (i, j) with given j */
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
                    if (pair.i2 != NULL_INDEX32)
                        result.push_back(std::move(pair));
            buffer_min_pairs_.erase(buffer_min_pairs_.begin());
        }
        else if (!buffer_min_pairs_.empty() && buffer_min_pairs_.begin()->first < d)
            throw MyException(0x3321786aU, "buffer_min_pairs_ contains degree < d");
        return result;
    }
    bool empty() const
    {
        return buffer_min_pairs_.empty();
    }

    /* Minimize `buffer_min_pairs_[d]` and maintain `pairs_` */
    void Minimize(const Mon1d& leads, int d);
    void Minimize(const Mon1d& leadsx, const MMod1d& leads, int d);

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const Mon1d& leads, const MonTrace1d& traces, const Mon& mon, const int1d& gen_degs);
    void AddToBuffers(const Mon1d& leadsx, const MonTrace1d& tracesx, const MMod1d& leads, const MonTrace1d& traces, const MMod& mon, const int1d& gen_degs, const int1d& v_degs);
    void init(const Mon1d& leads, const MonTrace1d& traces, const int1d& gen_degs, int t_min_buffer);
    void init(const Mon1d& leadsx, const MonTrace1d& tracesx, const MMod1d& leads, const MonTrace1d& traces, const int1d& gen_degs, const int1d& v_degs, int t_min_buffer);
};

class Groebner
{
private:
    using TypeIndexKey = uint32_t;
    using TypeIndices = std::unordered_map<TypeIndexKey, int1d>;
    friend class GroebnerMod;

private:
    GbCriPairs criticals_; /* Groebner basis of critical pairs */

    Poly1d data_;
    Mon1d leads_;         /* Leading monomials */
    MonTrace1d traces_;   /* Cache for fast divisibility test */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int1d gen_degs_; /* degree of generators */

public:
    Groebner() : criticals_(DEG_MAX) {}
    Groebner(int deg_trunc, int1d gen_degs) : criticals_(deg_trunc), gen_degs_(std::move(gen_degs)) {}

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    Groebner(int deg_trunc, int1d gen_degs, Poly1d polys, bool bDynamic = false);

private:
    static TypeIndexKey Key(const Mon& lead)
    {
        return TypeIndexKey{lead.back().g() + (lead.size() == 1 ? 0 : ((lead[size_t(lead.size() - 2)].g() + 1) << 16))};
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
        criticals_.Minimize(leads_, d);
        return criticals_.Criticals(d);
    }

    auto size() const
    {
        return data_.size();
    }
    const auto& operator[](size_t index) const
    {
        return data_[index];
    }
    void push_back(Poly g)
    {
        const Mon& m = g.GetLead();
        criticals_.AddToBuffers(leads_, traces_, m, gen_degs_);

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

    const auto& gen_degs() const
    {
        return gen_degs_;
    }

    /* Return trace of LM(gb[i]) */
    MonTrace GetTrace(size_t i) const
    {
        return traces_[i];
    }

public:
    /* Leadings[i] is the set of monomials that end with generator i.
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
        fm = gb.Reduce(std::move(fm));
        result.iaddP(fm, tmp);
    }
    return result;
}

inline Poly subsMGb(const Poly& poly, const Groebner& gb, const std::vector<Poly>& map)
{
    return subsMGbTpl(poly, gb, [&map](size_t i) { return map[i]; });
}

/********************************* Modules ****************************************/

class GroebnerMod
{
private:
    using TypeIndexKey = uint32_t;
    using TypeIndices = std::unordered_map<TypeIndexKey, int1d>;

private:
    Groebner* pGb_;
    GbCriPairs criticals_; /* Groebner basis of critical pairs */

    Mod1d data_;
    MMod1d leads_;        /* Leading monomials */
    MonTrace1d traces_;   /* Cache for fast divisibility test */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int1d v_degs_; /* degree of generators of modules */

public:
    GroebnerMod() : pGb_(nullptr), criticals_(DEG_MAX) {}
    GroebnerMod(Groebner* pGb, int deg_trunc, int1d v_degs) : pGb_(pGb), criticals_(deg_trunc), v_degs_(std::move(v_degs)) {}

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    GroebnerMod(Groebner* pGb, int deg_trunc, int1d v_degs, Mod1d polys, bool bDynamic = false);

private:
    static TypeIndexKey Key(const MMod& lead)
    {
        return TypeIndexKey{lead.v + (lead.m.size() > 0 ? ((lead.m.back().g() + 1) << 16) : 0)};
    }

public:  ////
    /**
     * Transform to a submodule.
     *
     * `rels` should be ordered by degree.
     */
    void ToSubMod(const Mod1d& rels, int deg, int1d& index_ind);

public: /* Getters and Setters */
    const auto& gb_pairs() const
    {
        return criticals_;
    }
    int deg_trunc() const
    {
        return criticals_.deg_trunc();
    }
    CriPair1d Criticals(int d)
    {
        criticals_.Minimize(pGb_->leads_, leads_, d);
        return criticals_.Criticals(d);
    }

    auto size() const
    {
        return data_.size();
    }
    const auto& operator[](size_t index) const
    {
        return data_[index];
    }
    void push_back(Mod g)
    {
        const MMod& m = g.GetLead();
        criticals_.AddToBuffers(pGb_->leads_, pGb_->traces_, leads_, traces_, m, pGb_->gen_degs(), v_degs_);

        leads_.push_back(m);
        traces_.push_back(m.m.Trace());
        indices_[Key(m)].push_back((int)data_.size());
        data_.push_back(std::move(g));
    }
    bool operator==(const GroebnerMod& rhs) const
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
    /* Leadings[i] is the set of monomials that end with v_i.
     * The result is used to generate a basis of the module.
     */
    MMod2d GetLeadings(size_t gens_size) const
    {
        MMod2d result;
        result.resize(gens_size);
        for (size_t i = 0; i < data_.size(); ++i)
            result[leads_[i].v].push_back(leads_[i]);
        return result;
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(const MMod& mon) const;

    Mod Reduce(Mod poly) const;

    /**
     * Comsume relations from 'rels` and `gb.criticals_` in degree `<= deg`
     */
    void AddRels(const Mod1d& rels, int deg);
};

}  // namespace alg2

#endif /* GROEBNER_H */
