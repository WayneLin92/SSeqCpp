/** \file groebnerZ.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNERZ_H
#define GROEBNERZ_H

#include "algebras/algebrasZ.h"
#include "algebras/groebner.h"
#include <unordered_map>

/* extension of namespace alg in algebras.h */
namespace algZ {

class Groebner;
class GroebnerMod;

inline constexpr uint32_t NULL_INDEX32 = ~uint32_t(0);
inline constexpr uint32_t FLAG_INDEX_X = uint32_t(1) << 31;

struct CriPair
{
    uint32_t i1 = NULL_INDEX32, i2 = NULL_INDEX32;
    Mon m1, m2;
    int O = -1;
    MonTrace trace_m2 = 0; /* = Trace(m2) */

    CriPair() = default;
    /* Compute the pair for two leading monomials. */
    static void SetFromLM(CriPair& result, const Mon& lead1, const Mon& lead2, int O1, int O2, int i, int j, const AdamsDeg1d& gen_degs);

    /* Return `m1 * gb[i1] + m2 * gb[i2]` */
    void SijP(const Groebner& gb, Poly& result, Poly& tmp) const;
    void SijMod(const Groebner& gb, const GroebnerMod& gbm, Mod& result, Mod& tmp) const;
};
using CriPair1d = std::vector<CriPair>;
using CriPair2d = std::vector<CriPair1d>;

/* Groebner basis of critical pairs */
class GbCriPairs
{
public:
    int t_trunc_;                                                        /* Truncation degree */
    std::map<int, std::unordered_map<int, CriPair1d>> buffer_min_pairs_; /* To generate `buffer_min_pairs_` and for computing Sij */
    std::vector<std::pair<int, CriPair>> new_pairs__;                    /* tmp variable to be used in functions */

public:
    GbCriPairs(int t_trunc) : t_trunc_(t_trunc) {}
    int deg_trunc() const
    {
        return t_trunc_;
    }
    CriPair1d Criticals(int d)
    {
        CriPair1d result;
        if (!buffer_min_pairs_.empty() && buffer_min_pairs_.begin()->first == d) {
            auto& b_min_pairs_d = buffer_min_pairs_.begin()->second;
            for (auto& [_, pairs] : b_min_pairs_d)
                for (auto& pair : pairs)
                    if (pair.i2 != NULL_INDEX32)
                        result.push_back(std::move(pair));
            buffer_min_pairs_.erase(buffer_min_pairs_.begin());
        }
        else if (!buffer_min_pairs_.empty() && buffer_min_pairs_.begin()->first < d)
            throw ErrorIdMsg(0x3321786aU, "buffer_min_pairs_ contains degree < d");
        return result;
    }
    int NextD() const
    {
        return buffer_min_pairs_.empty() ? -1 : buffer_min_pairs_.begin()->first;
    }
    bool empty() const
    {
        return buffer_min_pairs_.empty();
    }

    void Reset()
    {
        buffer_min_pairs_.clear();
    }

    void ClearBuffer()
    {
        buffer_min_pairs_.clear();
    }

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const Mon1d& leads, const MonTrace1d& traces, const int1d& leads_O, const Mon& mon, int O, const AdamsDeg1d& gen_degs);
    void AddToBuffers(const Mon1d& leadsx, const int1d& leadsx_O, const MMod1d& leads, const int1d& leads_O, const MMod& mon, int O, const AdamsDeg1d& gen_degs, const AdamsDeg1d& v_degs);
    void AddToBuffersX(const Mon1d& leadsx, const int1d& leadsx_O, const MMod1d& leads, const int1d& leads_O, const AdamsDeg1d& gen_degs, const AdamsDeg1d& v_degs, size_t i_start);
};

int NextO(const ut::map_seq2d<int, 0>& possEinf, int t_max, int stem, int O_min);

/* 2 is considered as the first generator */
class Groebner
{
private:
    using TypeIndexKey = uint32_t;
    friend class GroebnerMod;

private:
    GbCriPairs criticals_; /* Groebner basis of critical pairs */

    Poly1d data_;
    Mon1d leads_;                                                /* Leading monomials */
    MonTrace1d traces_;                                          /* Cache for fast divisibility test */
    int1d data_O_;                                               /* Cache for certainty of data_ */
    std::unordered_map<TypeIndexKey, int1d> leads_group_by_key_; /* Cache for fast divisibility test */
    int2d leads_group_by_last_gen_;                              /* Cache for generating a basis */
    std::map<AdamsDeg, int1d> leads_group_by_deg_;               /* Cache for generating a basis */
    std::map<int, int1d> leads_group_by_t_;                      /* Cache for iteration */

    AdamsDeg1d gen_degs_; /* degree of generators */
    int2d nodes_gen_2tor_degs_; /* 2 torsion degree of generators */
    std::vector<size_t> nodes_gen_size_;
    std::vector<size_t> nodes_data_size_;

public:
    Groebner() : criticals_(DEG_MAX), gen_degs_({AdamsDeg(1, 1)}), nodes_gen_2tor_degs_({{FIL_MAX + 1}}) {}
    Groebner(int t_trunc, AdamsDeg1d gen_degs) : criticals_(t_trunc), gen_degs_(std::move(gen_degs)), nodes_gen_2tor_degs_({})
    {
        for (AdamsDeg deg : gen_degs_) {
            if (deg == AdamsDeg(1, 1))
                nodes_gen_2tor_degs_.back().push_back(FIL_MAX + 1);
            else
                nodes_gen_2tor_degs_.back().push_back((deg.stem() + 5) / 2);
        }
    }

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    Groebner(int t_trunc, AdamsDeg1d gen_degs, Poly1d polys);

private:
    static TypeIndexKey Key(const Mon& lead)
    {
        return TypeIndexKey{lead.backg() + (lead.backg2p1() << 16)};
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
        // criticals_.Minimize(leads_, data_O_, d, gen_degs_);
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

    void push_back_data(Poly g, AdamsDeg deg)
    {
        const Mon& m = g.GetLead();

        leads_.push_back(m);
        traces_.push_back(m.Trace());
        data_O_.push_back(g.UnknownFil());
        int index = (int)data_.size();
        leads_group_by_key_[Key(m)].push_back(index);
        leads_group_by_deg_[deg].push_back(index);
        leads_group_by_t_[deg.t].push_back(index);
        uint32_t backg = m.backg();
        ut::get(leads_group_by_last_gen_, backg).push_back(index);

        if (g.data.size() == 1) {
            if (m.frontg() == backg && backg > 0 && m.m().begin()->e_masked() == 1) {
                size_t i = backg;
                nodes_gen_2tor_degs_.back()[i] = std::min(m.c(), nodes_gen_2tor_degs_.back()[i]);
            }
        }

        data_.push_back(std::move(g));
    }

    void push_back_data_buffer(Poly g, AdamsDeg deg)
    {
        criticals_.AddToBuffers(leads_, traces_, data_O_, g.GetLead(), g.UnknownFil(), gen_degs_);
        push_back_data(std::move(g), deg);
    }

    void ResetRels()
    {
        criticals_.Reset();
        leads_.clear();
        traces_.clear();
        data_O_.clear();
        leads_group_by_key_.clear();
        leads_group_by_deg_.clear();
        leads_group_by_t_.clear();
        leads_group_by_last_gen_.clear();
        data_.clear();

        for (size_t i = 1; i < gen_degs_.size(); ++i)
            nodes_gen_2tor_degs_.back()[i] = (gen_degs_[i].stem() + 5) / 2;  //// TODO: modify
    }

    /* Restore the algebra to a previous status */
    void AddNode();
    void PopNode();
    void debug_print() const;

    bool operator==(const Groebner& rhs) const
    {
        return data_ == rhs.data_;
    }

    const auto& data() const
    {
        return data_;
    }

    const auto& leads() const
    {
        return leads_;
    }

    const auto& gen_degs() const
    {
        return gen_degs_;
    }

    const auto& gen_2tor_degs() const
    {
        return nodes_gen_2tor_degs_.back();
    }

    const auto& leads_group_by_deg() const
    {
        return leads_group_by_deg_;
    }

    const auto& leads_group_by_t() const
    {
        return leads_group_by_t_;
    }

    auto OutputForDatabase() const -> const std::map<AdamsDeg, Poly1d>
    {
        std::map<AdamsDeg, Poly1d> result;
        for (auto& [deg, indices] : leads_group_by_deg_) {
            for (int i : indices)
                result[deg].push_back(data_[i]);
        }
        return result;
    }

    /* Return trace of LM(gb[i]) */
    MonTrace GetTrace(size_t i) const
    {
        return traces_[i];
    }

    Mon Gen(uint32_t gen_id, uint32_t exp = 1) const
    {
        return Mon::Gen(gen_id, exp, gen_degs_[gen_id].s * exp, gen_degs_[gen_id].stem() % 2 == 0);
    }

public:
    /* Compute the basis in `deg` assuming all basis have been calculated in lower total degrees
     * The relations in `deg` are ignored
     */
    bool IsNewBaseByLastGen(const Mon& mon, uint32_t last_gen) const;

    /* Return relations with leadings with the given deg */
    Poly1d RelsLF(AdamsDeg deg) const;

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(const Mon& mon, int eff_min) const;
    int IndexOfDivisibleLeadingV2(const Mon& mon) const;

    /* This function does not decrease the certainty of poly */
    Poly Reduce(Poly poly) const;
    Poly ReduceForGbRel(Poly poly) const;
    /* This reduce poly by gb with all certainties until the lead is in the basis */
    Poly ReduceV2(Poly poly) const;

    void AddGen(AdamsDeg deg)
    {
        gen_degs_.push_back(deg);
        if (deg == AdamsDeg(1, 1))
            nodes_gen_2tor_degs_.back().push_back(FIL_MAX + 1);
        else
            nodes_gen_2tor_degs_.back().push_back((deg.stem() + 5) / 2); ////TODO: modify
    }

    /**
     * Comsume relations from 'rels` and `gb.criticals_` in degree `<= deg`
     */
    void AddRels(Poly1d rels, int deg, const ut::map_seq2d<int, 0>& possEinf);

    void SimplifyRels();
    /* Apply differential monommial orderings to old database */
    void SimplifyRelsReorder(const ut::map_seq2d<int, 0>& possEinf);
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
        Poly fm = Poly::twoTo(m.c());
        for (auto p = m.m().begin(); p != m.m().end(); ++p) {
            powP(map(p->g()), p->e_masked(), gb, tmp_prod, tmp);
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

private:
    const Groebner* pGb_;
    size_t old_pGb_size_;
    GbCriPairs criticals_; /* Groebner basis of critical pairs */

    Mod1d data_;
    MMod1d leads_;                                               /* Leading monomials */
    MonTrace1d traces_;                                          /* Cache for fast divisibility test */
    int1d data_O_;                                               /* Cache for certainty of data_ */
    std::unordered_map<TypeIndexKey, int1d> leads_group_by_key_; /* Cache for fast divisibility test */
    int2d leads_group_by_v_;                                     /* Cache for generating a basis */
    std::map<AdamsDeg, int1d> leads_group_by_deg_;               /* Cache for iteration */
    std::map<int, int1d> leads_group_by_t_;                      /* Cache for iteration */

    AdamsDeg1d v_degs_; /* degree of generators of modules */

    std::vector<size_t> nodes_gen_size_;
    std::vector<size_t> nodes_data_size_;

public:
    GroebnerMod() : pGb_(nullptr), criticals_(DEG_MAX), old_pGb_size_(0) {}
    GroebnerMod(const Groebner* pGb, int deg_trunc, AdamsDeg1d v_degs) : pGb_(pGb), criticals_(deg_trunc), v_degs_(std::move(v_degs)), old_pGb_size_(pGb->leads_.size()) {}

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    GroebnerMod(const Groebner* pGb, int deg_trunc, AdamsDeg1d v_degs, Mod1d polys);

private:
    static TypeIndexKey Key(const MMod& lead)
    {
        return TypeIndexKey{lead.v + (lead.m ? ((lead.m.backg() + 1) << 16) : 0)};
    }

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
    void push_back_data(Mod g, AdamsDeg deg)
    {
        const MMod& m = g.GetLead();

        leads_.push_back(m);
        traces_.push_back(m.m.Trace());
        data_O_.push_back(g.UnknownFil());
        int index = (int)data_.size();
        leads_group_by_key_[Key(m)].push_back(index);
        leads_group_by_deg_[deg].push_back(index);
        leads_group_by_t_[deg.t].push_back(index);
        ut::get(leads_group_by_v_, m.v).push_back(index);
        data_.push_back(std::move(g));
    }
    /* This is used for initialization */
    void push_back_data_buffer(Mod g, AdamsDeg deg)
    {
        criticals_.AddToBuffers(pGb_->leads_, pGb_->data_O_, leads_, data_O_, g.GetLead(), g.UnknownFil(), pGb_->gen_degs(), v_degs_);
        push_back_data(std::move(g), deg);
    }

    void ResetRels()
    {
        old_pGb_size_ = pGb_->leads_.size();
        criticals_.Reset();
        leads_.clear();
        traces_.clear();
        data_O_.clear();
        leads_group_by_key_.clear();
        leads_group_by_deg_.clear();
        leads_group_by_t_.clear();
        leads_group_by_v_.clear();
        data_.clear();
    }

    void AddNode();
    void PopNode();
    void debug_print() const;  ////

    bool operator==(const GroebnerMod& rhs) const
    {
        return data_ == rhs.data_;
    }

    auto& data() const
    {
        return data_;
    }

    const auto& leads() const
    {
        return leads_;
    }

    auto& v_degs() const
    {
        return v_degs_;
    }

    const auto& gen_2tor_degs() const
    {
        return pGb_->gen_2tor_degs();
    }

    const auto& leads_group_by_deg() const
    {
        return leads_group_by_deg_;
    }

    const auto& leads_group_by_t() const
    {
        return leads_group_by_t_;
    }

    const std::map<AdamsDeg, Mod1d> OutputForDatabase() const
    {
        std::map<AdamsDeg, Mod1d> result;
        for (auto& [deg, indices] : leads_group_by_deg_) {
            for (int i : indices)
                result[deg].push_back(data_[i]);
        }
        return result;
    }

    /* Return trace of LM(gb[i]) */
    MonTrace GetTrace(size_t i) const
    {
        return traces_[i];
    }

    MMod Gen(uint32_t v_id) const
    {
        return MMod(Mon(), v_id, v_degs_[v_id].s);
    }

public:
    /* Compute the basis in `deg` based on the basis of pGb_ */
    bool IsNewBaseByV(const Mon& mon, uint32_t v) const;

    /* Return relations with leadings with the given deg */
    Mod1d RelsLF(AdamsDeg deg) const;

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(const MMod& mon, int eff_min) const;
    /* Return -1 if not found
     * Return the index of gb with the biggest effective number
     */
    int IndexOfDivisibleLeadingV2(const MMod& mon) const;

    /* This function does not decrease the certainty of poly */
    Mod Reduce(Mod poly) const;
    Mod ReduceForGbRel(Mod poly) const;
    /* This reduce poly by gb with all certainties until the lead is in the basis */
    Mod ReduceV2(Mod poly) const;

    void AddGen(AdamsDeg deg)
    {
        v_degs_.push_back(deg);
    }

    /**
     * Comsume relations from 'rels` and `gb.criticals_` in degree `<= deg`
     */
    void AddRels(Mod1d rels, int deg, const ut::map_seq2d<int, 0>& possEinf);

    void SimplifyRels();
    void SimplifyRelsReorder(const ut::map_seq2d<int, 0>& possEinf);
};

inline bool IsValidRel(const Poly& poly)
{
    return poly && !poly.GetLead().IsUnKnown();
}

inline bool IsValidRel(const Mod& x)
{
    return x && !x.GetLead().IsUnKnown();
}

}  // namespace algZ

#endif /* GROEBNER_H */
