#ifndef GROEBNERZV2_H
#define GROEBNERZV2_H

#include "algebras/algebrasZ.h"
#include "algebras/groebner.h"
#include <unordered_map>

namespace algZ {

struct PiBase
{
    algZ::Mon1d pi_basis;
    int2d Einf;
};

using PiBasis = std::map<AdamsDeg, PiBase>;
using PiBasis1d = std::vector<PiBasis>;

/* 2 is considered as the first generator */
class PiGroebner
{
private:
    using TypeIndexKey = uint32_t;
    friend class PiGroebnerMod;

private:
    int t_trunc_ = DEG_MAX;
    Poly1d data_;
    Mon1d leads_;                                                /* Leading monomials */
    MonTrace1d traces_;                                          /* Cache for fast divisibility test */
    int1d leads_O_;                                              /* Cache for certainty of data_ */
    std::unordered_map<TypeIndexKey, int1d> leads_group_by_key_; /* Cache for fast divisibility test */
    int2d leads_group_by_last_gen_;                              /* Cache for generating a basis */
    std::map<AdamsDeg, int1d> leads_group_by_deg_;               /* Cache for generating a basis */
    std::map<int, int1d> leads_group_by_t_;                      /* Cache for iteration */
    std::vector<size_t> nodes_data_size_;

    AdamsDeg1d gen_degs_ = {AdamsDeg(1, 1)};      /* degree of generators */
    alg2::Poly1d gen_Einf = {alg2::Poly::Gen(0)}; /* Projection onto the E_infty page */
    int2d nodes_gen_2tor_degs_ = {{FIL_MAX + 1}}; /* 2 torsion degree of generators */
    std::vector<size_t> nodes_gen_degs_size_;

    PiBasis1d nodes_basis_ = {{{AdamsDeg(0, 0), {{algZ::Mon()}, {{0}}}}}};

public:
    PiGroebner() {}
    PiGroebner(int t_trunc, AdamsDeg1d gen_degs) : t_trunc_(t_trunc), gen_degs_(std::move(gen_degs)), nodes_gen_2tor_degs_({{}})
    {
        for (AdamsDeg deg : gen_degs_) {
            if (deg == AdamsDeg(1, 1))
                nodes_gen_2tor_degs_[0].push_back(FIL_MAX + 1);
            else
                nodes_gen_2tor_degs_[0].push_back((deg.stem() + 5) / 2);
        }
    }

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    PiGroebner(int t_trunc, AdamsDeg1d gen_degs, Poly1d polys);

private:
    static TypeIndexKey Key(const Mon& lead)
    {
        return TypeIndexKey{lead.backg() + (lead.backg2p1() << 16)};
    }

public: /* Getters and Setters */
    int deg_trunc() const
    {
        return t_trunc_;
    }
    auto size() const
    {
        return data_.size();
    }
    const auto& operator[](size_t index) const
    {
        return data_[index];
    }

    bool operator==(const PiGroebner& rhs) const
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

    const PiBase* GetRecentPiBasis(AdamsDeg deg) const
    {
        for (auto p = nodes_basis_.rbegin(); p != nodes_basis_.rend(); ++p)
            if (p->find(deg) != p->end())
                return &p->at(deg);
        return nullptr;
    }

    void push_back_data(Poly g, AdamsDeg deg)
    {
        const Mon& m = g.GetLead();

        leads_.push_back(m);
        traces_.push_back(m.Trace());
        leads_O_.push_back(g.UnknownFil());
        int index = (int)data_.size();
        leads_group_by_key_[Key(m)].push_back(index);
        leads_group_by_deg_[deg].push_back(index);
        leads_group_by_t_[deg.t].push_back(index);
        uint32_t backg = m.backg();
        ut::push_back(leads_group_by_last_gen_, (size_t)backg, index);

        if (g.data.size() == 1) {
            if (m.frontg() == backg && backg > 0 && m.m().begin()->e_masked() == 1) {
                size_t index = backg;
                nodes_gen_2tor_degs_.back()[index] = std::min(m.c(), nodes_gen_2tor_degs_.back()[index]);
            }
        }

        data_.push_back(std::move(g));
    }

    void ResetRels()
    {
        leads_.clear();
        traces_.clear();
        leads_O_.clear();
        leads_group_by_key_.clear();
        leads_group_by_deg_.clear();
        leads_group_by_t_.clear();
        leads_group_by_last_gen_.clear();
        data_.clear();

        for (size_t i = 1; i < gen_degs_.size(); ++i)
            nodes_gen_2tor_degs_.back()[i] = (gen_degs_[i].stem() + 5) / 2;  //// TODO: modify
    }

    /* Restore the algebra to a previous status */
    void Pop(size_t gen_size, size_t rel_size);
    void debug_print() const;

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
            nodes_gen_2tor_degs_.back().push_back((deg.stem() + 5) / 2);
    }

    /**
     * Comsume relations from 'rels` and `gb.criticals_` in degree `<= deg`
     */
    void AddRel(Poly rel, int deg, const ut::map_seq2d<int, 0>& possEinf);
    void AddRels(Poly1d rels, int deg, const ut::map_seq2d<int, 0>& possEinf);

    void SimplifyRels(const ut::map_seq2d<int, 0>& possEinf);
    void SimplifyRelsReorder(const ut::map_seq2d<int, 0>& possEinf);

    void Sync();
};

/**
 * A fast algorithm that computes
 * `poly ** n` modulo `gb`
 */
void powP(const Poly& poly, int n, const PiGroebner& gb, Poly& result, Poly& tmp);
inline Poly pow(const Poly& poly, int n, const PiGroebner& gb)
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
Poly subsMGbTpl(const Poly& poly, const PiGroebner& gb, const FnMap& map)
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

inline Poly subsMGb(const Poly& poly, const PiGroebner& gb, const std::vector<Poly>& map)
{
    return subsMGbTpl(poly, gb, [&map](size_t i) { return map[i]; });
}

/********************************* Modules ****************************************/

class PiGroebnerMod
{
private:
    using TypeIndexKey = uint32_t;

private:
    int t_trunc_;
    PiGroebner* pGb_;
    size_t old_pGb_size_;

    Mod1d data_;
    MMod1d leads_;                                               /* Leading monomials */
    MonTrace1d traces_;                                          /* Cache for fast divisibility test */
    int1d leads_O_;                                              /* Cache for certainty of data_ */
    std::unordered_map<TypeIndexKey, int1d> leads_group_by_key_; /* Cache for fast divisibility test */
    int2d leads_group_by_v_;                                     /* Cache for generating a basis */
    std::map<AdamsDeg, int1d> leads_group_by_deg_;               /* Cache for iteration */
    std::map<int, int1d> leads_group_by_t_;                      /* Cache for iteration */

    AdamsDeg1d v_degs_; /* degree of generators of modules */

public:
    PiGroebnerMod() : pGb_(nullptr), t_trunc_(DEG_MAX), old_pGb_size_(0) {}
    PiGroebnerMod(PiGroebner* pGb, int deg_trunc, AdamsDeg1d v_degs) : pGb_(pGb), t_trunc_(deg_trunc), v_degs_(std::move(v_degs)), old_pGb_size_(pGb->leads_.size()) {}

    /* Initialize from `polys` which already forms a Groebner basis. The instance will be in const mode. */
    PiGroebnerMod(PiGroebner* pGb, int deg_trunc, AdamsDeg1d v_degs, Mod1d polys, bool bDynamic = false);

private:
    static TypeIndexKey Key(const MMod& lead)
    {
        return TypeIndexKey{lead.v + (lead.m ? ((lead.m.backg() + 1) << 16) : 0)};
    }

public:
    /**
     * Transform to a submodule.
     *
     * `rels` should be ordered by degree.
     */
    // void ToSubMod(const Mod1d& rels, int deg, int1d& index_ind);

public: /* Getters and Setters */
    int deg_trunc() const
    {
        return t_trunc_;
    }
    auto size() const
    {
        return data_.size();
    }
    const auto& operator[](size_t index) const
    {
        return data_[index];
    }
    void push_back_data_init(Mod g, AdamsDeg deg)
    {
        const MMod& m = g.GetLead();

        leads_.push_back(m);
        traces_.push_back(m.m.Trace());
        leads_O_.push_back(g.UnknownFil());
        int index = (int)data_.size();
        leads_group_by_key_[Key(m)].push_back(index);
        leads_group_by_deg_[deg].push_back(index);
        leads_group_by_t_[deg.t].push_back(index);
        ut::push_back(leads_group_by_v_, (size_t)m.v, index);
        data_.push_back(std::move(g));
    }
    /* This is used for initialization */
    void push_back_data(Mod g, AdamsDeg deg)
    {
        push_back_data_init(std::move(g), deg);
    }

    void ResetRels()
    {
        old_pGb_size_ = pGb_->leads_.size();
        leads_.clear();
        traces_.clear();
        leads_O_.clear();
        leads_group_by_key_.clear();
        leads_group_by_deg_.clear();
        leads_group_by_t_.clear();
        leads_group_by_v_.clear();
        data_.clear();
    }

    void Pop(size_t gen_size, size_t rel_size);
    void debug_print() const;  ////

    bool operator==(const PiGroebnerMod& rhs) const
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
    /*
     * Return the index of gb with the biggest effective number
     * Return -1 if not found
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

    void SimplifyRels(const ut::map_seq2d<int, 0>& possEinf);
    void SimplifyRelsReorder(const ut::map_seq2d<int, 0>& possEinf);
};

/* Return if x is nonzero and the leading term is known */
template<typename T>
bool IsValidRel(const T& x)
{
    return x && !x.GetLead().IsUnKnown();
}

}  // namespace algZ

#endif /* GROEBNER_H */