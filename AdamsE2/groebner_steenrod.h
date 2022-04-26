/** \file groebner_steenrod.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_STEENROD_H
#define GROEBNER_STEENROD_H

#include "algebras/steenrod.h"
#include <execution>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace steenrod {

struct CriPairMRes
{
    int i1 = -1, i2 = -1;
    MMilnor m1, m2;

    /* Compute the pair for two leading monomials. */
    CriPairMRes() = default;
    static void SetFromLM(CriPairMRes& result, MMilnor lead1, MMilnor lead2, int i, int j)
    {
        MMilnor gcd = gcdLF(lead1, lead2);
        result.m1 = divLF(lead2, gcd);
        result.m2 = divLF(lead1, gcd);
        result.i1 = i;
        result.i2 = j;
    }
    static CriPairMRes Single(MMilnor m2, int j)
    {
        CriPairMRes result;
        result.m2 = m2;
        result.i1 = -1;
        result.i2 = j;
        return result;
    }
};
using CriPairMRes1d = std::vector<CriPairMRes>;
using CriPairMRes2d = std::vector<CriPairMRes1d>;
using CriPairMRes3d = std::vector<CriPairMRes2d>;
using PCriPairMRes1d = std::vector<CriPairMRes*>;
using PCriPairMRes2d = std::vector<PCriPairMRes1d>;

struct AdamsDeg
{
    int s, t;
    bool operator<(AdamsDeg rhs) const
    {
        if (t < rhs.t)
            return true;
        if (t > rhs.t)
            return false;
        if (s > rhs.s)
            return true;
        return false;
    }
};

/* Groebner basis of critical pairs */
class GbCriPairsMRes
{
    using TypeRedSing = std::vector<std::vector<std::unordered_set<uint64_t>>>;

private:
    int deg_trunc_;                                                           /* Truncation degree */
    CriPairMRes3d pairs_;                                                     /* `pairs_[s][j]` is the set of pairs (i, j) with given j in degree s */
    TypeRedSing redundent_singles_;                                           /* `redundent_singles_[s][i]` is the set of generators that should not be multiplied by leads[s][i] in reduction */
    std::map<AdamsDeg, CriPairMRes2d> buffer_min_pairs_;                      /* `buffer_min_pairs_[st]` To generate minimal pairs to compute Sij */
    std::map<AdamsDeg, std::unordered_set<uint64_t>> buffer_redundent_pairs_; /* Used to minimize `buffer_min_pairs_` */
    std::map<AdamsDeg, CriPairMRes1d> buffer_singles_;                        /* For computing Sj. `buffer_singles_` stores indices of singles_ */

public:
    GbCriPairsMRes(int d_trunc) : deg_trunc_(d_trunc) {}

    int deg_trunc() const
    {
        return deg_trunc_;
    }
    AdamsDeg next_st() const
    {
        AdamsDeg st = buffer_min_pairs_.empty() ? AdamsDeg{0, 1024} : buffer_min_pairs_.begin()->first;
        if (!buffer_singles_.empty())
            st = std::min(st, buffer_singles_.begin()->first);
        if (st.t == 1024)
            st = AdamsDeg{1, 1};
        return st;
    }
    bool empty_pairs_for_gb() const
    {
        return buffer_min_pairs_.empty() && buffer_singles_.empty();
    }
    bool empty_min_pairs_for_gb(AdamsDeg st) const
    {
        return buffer_min_pairs_.find(st) == buffer_min_pairs_.end();
    }
    void resize_pairs(size_t s)
    {
        pairs_.resize(s);
    }
    CriPairMRes1d pairs_for_gb(AdamsDeg st)
    {
        CriPairMRes1d result;
        if (buffer_singles_.find(st) != buffer_singles_.end()) {
            std::swap(result, buffer_singles_.at(st));
            buffer_singles_.erase(st);
        }
        if (buffer_min_pairs_.find(st) != buffer_min_pairs_.end()) {
            for (int j = 0; j < (int)buffer_min_pairs_.at(st).size(); ++j)
                for (auto& pair : buffer_min_pairs_.at(st)[j])
                    if (pair.i2 != -1)
                        result.push_back(std::move(pair));
            buffer_min_pairs_.erase(st);
        }

        return result;
    }

    /* Minimize `buffer_min_pairs_[t]` and maintain `pairs_` */
    void Minimize(const MMod1d& leads, AdamsDeg st);

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const MMod1d& leads, MMod mon, int d_mon_base, int s);

    void init(const MMod2d& leads, const array2d& basis_degrees);
};

struct DataMRes
{
    Mod x1, x2;
    DataMRes() {}
    DataMRes(Mod x1_, Mod x2_) : x1(std::move(x1_)), x2(std::move(x2_)) {}
    DataMRes operator+(const DataMRes& rhs) const
    {
        return DataMRes(x1 + rhs.x1, x2 + rhs.x2);
    }
    DataMRes& operator+=(const DataMRes& rhs)
    {
        x1 += rhs.x1;
        x2 += rhs.x2;
        return *this;
    }
};

using DataMRes1d = std::vector<DataMRes>;
using DataMRes2d = std::vector<DataMRes1d>;

class GroebnerMRes
{
private:
    using TypeIndices = std::vector<std::unordered_map<uint32_t, array>>;

private:
    GbCriPairsMRes gb_pairs_; /* Groebner basis of critical pairs */

    DataMRes2d data_;
    MMod2d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    array2d basis_degrees_; /* `basis_degrees[s][i]` is the degree of v_{s,i} */

public:
    GroebnerMRes(int deg_trunc, array2d basis_degrees) : gb_pairs_(deg_trunc), basis_degrees_(std::move(basis_degrees)) {}

    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    GroebnerMRes(int deg_trunc, DataMRes2d data, array2d basis_degrees) : gb_pairs_(deg_trunc), data_(std::move(data)), basis_degrees_(std::move(basis_degrees))
    {
        if (basis_degrees_.empty())
            basis_degrees_.push_back({});
        if (basis_degrees_[0].empty())
            basis_degrees_[0].push_back(0);

        leads_.resize(data_.size());
        indices_.resize(data_.size());
        gb_pairs_.resize_pairs(data_.size());

        for (size_t s = 0; s < data_.size(); ++s) {
            for (int j = 0; j < (int)data_[s].size(); ++j) {
                leads_[s].push_back(data_[s][j].x1.GetLead());
                indices_[s][Key(data_[s][j].x1.GetLead())].push_back(j);
            }
        }
        gb_pairs_.init(leads_, basis_degrees_);
    }

private:
    static uint32_t Key(MMod lead)
    {
        return uint32_t(lead.v());
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(MMod mon, int s) const
    {
        auto key = uint32_t(mon.v());
        auto p = indices_[s].find(key);
        if (p != indices_[s].end())
            for (int k : p->second)
                if (divisibleLF(leads_[s][k], mon))
                    return k;
        return -1;
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeadingX2(MMod mon, int s) const
    {
        for (size_t k = 0; k < leads_[s].size(); ++k)
            if (data_[s][k].x2.GetLead().v_raw() == mon.v_raw() && divisibleLF(data_[s][k].x2.GetLead(), mon))
                return (int)k;
        return -1;
    }

public: /* Getters and Setters */
    const GbCriPairsMRes& gb_pairs() const
    {
        return gb_pairs_;
    }
    int deg_trunc() const
    {
        return gb_pairs_.deg_trunc();
    }
    /* This function will erase `gb_pairs_.buffer_min_pairs[t]` */
    CriPairMRes1d pairs(AdamsDeg st)
    {
        return gb_pairs_.pairs_for_gb(st);
    }
    const array2d& basis_degs() const
    {
        return basis_degrees_;
    }
    auto size() const
    {
        return data_.size();
    }
    auto& operator[](size_t index) const
    {
        return data_[index];
    }
    void resize_data(int s)
    {
        if (data_.size() <= (size_t)s) {
            size_t sp1 = size_t(s + 1); /* s plus 1 */
            data_.resize(sp1);
            leads_.resize(sp1);
            indices_.resize(sp1);
            gb_pairs_.resize_pairs(sp1);
        }
    }
    void push_back(DataMRes g, int s)
    {
        MMod mv = g.x1.GetLead();
        gb_pairs_.AddToBuffers(leads_[s], mv, basis_degrees_[s][mv.v()], s);  // TODO: modify the counterpart in algebras/groebner.h

        leads_[s].push_back(mv);
        indices_[s][Key(mv)].push_back((int)data_[s].size());

        data_[s].push_back(std::move(g));
    }

    /* Add x2 + v_{s+1,i} */
    void push_back_kernel(Mod x2, DataMRes1d& rels, int s, int t)
    {
        size_t sp1 = size_t(s + 1); /* s plus 1 */
        if (basis_degrees_.size() <= sp1) {
            size_t sp2 = size_t(s + 2);
            basis_degrees_.resize(sp2);
        }
        basis_degrees_[sp1].push_back(t);

        Mod x3 = MMod(MMilnor(), basis_degrees_[sp1].size() - 1);

        DataMRes g = DataMRes(std::move(x2), std::move(x3));
        rels.push_back(g);
        push_back(std::move(g), s);
    }
    void MinimizePairs(AdamsDeg st)
    {
        if (!gb_pairs_.empty_min_pairs_for_gb(st))
            gb_pairs_.Minimize(leads_[st.s], st);
    }

    const auto& data() const
    {
        return data_;
    }

public:
    DataMRes Reduce(CriPairMRes& p, int s) const
    {
        DataMRes result;
        size_t sp1 = size_t(s + 1);

        Milnor tmp_a;
        Mod tmp_m1;
        Mod tmp_m2;
        tmp_a.data.reserve(50);
        tmp_m1.data.reserve(100);
        tmp_m2.data.reserve(100);

        if (p.i1 >= 0) {
            result.x1.iaddmul(p.m1, data_[s][p.i1].x1, tmp_a, tmp_m1, tmp_m2).iaddmul(p.m2, data_[s][p.i2].x1, tmp_a, tmp_m1, tmp_m2);
            result.x2.iaddmul(p.m1, data_[s][p.i1].x2, tmp_a, tmp_m1, tmp_m2).iaddmul(p.m2, data_[s][p.i2].x2, tmp_a, tmp_m1, tmp_m2);
        }
        else {
            result.x1.iaddmul(p.m2, data_[s][p.i2].x1, tmp_a, tmp_m1, tmp_m2);
            result.x2.iaddmul(p.m2, data_[s][p.i2].x2, tmp_a, tmp_m1, tmp_m2);
        }

        size_t index;
        index = 0;
        while (index < result.x1.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x1.data[index], s);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x1.data[index], data_[s][gb_index].x1.data[0]);
                result.x1.iaddmul(m, data_[s][gb_index].x1, tmp_a, tmp_m1, tmp_m2);
                result.x2.iaddmul(m, data_[s][gb_index].x2, tmp_a, tmp_m1, tmp_m2);
            }
            else
                ++index;
        }
        index = 0;
        while (index < result.x2.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x2.data[index], s + 1);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x2.data[index], data_[sp1][gb_index].x1.data[0]);
                result.x2.iaddmul(m, data_[sp1][gb_index].x1, tmp_a, tmp_m1, tmp_m2);
            }
            else
                ++index;
        }

        return result;
    }

public:
    static GroebnerMRes load(const std::string& filename, int t_trunc);
};

/**
 * Comsume relations from 'rels` and `gb.gb_pairs_` in degree `<= deg`
 * `min_gb` stores the minimal generating set of gb.
 * return the dimension of the calculated range for debugging.
 */
size_t AddRelsMRes(GroebnerMRes& gb, const Mod1d& rels, int deg);

}  // namespace steenrod

#endif