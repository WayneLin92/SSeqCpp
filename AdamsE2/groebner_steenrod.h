/** \file groebner_steenrod.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_STEENROD_H
#define GROEBNER_STEENROD_H

#include "algebras/steenrod.h"
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace steenrod {

struct CPMilnor
{
    int i1 = -1, i2 = -1;
    MMilnor m1, m2;

    /* Compute the pair for two leading monomials. */
    CPMilnor() = default;
    static void SetFromLM(CPMilnor& result, MMilnor lead1, MMilnor lead2, int i, int j)
    {
        MMilnor gcd = gcdLF(lead1, lead2);
        result.m1 = divLF(lead2, gcd);
        result.m2 = divLF(lead1, gcd);
        result.i1 = i;
        result.i2 = j;
    }
    static CPMilnor Single(MMilnor m2, int j)
    {
        CPMilnor result;
        result.m2 = m2;
        result.i1 = -1;
        result.i2 = j;
        return result;
    }
};
using CPMilnor1d = std::vector<CPMilnor>;
using CPMilnor2d = std::vector<CPMilnor1d>;
using CPMilnor3d = std::vector<CPMilnor2d>;
using PtCPMilnor1d = std::vector<CPMilnor*>;
using PtCPMilnor2d = std::vector<PtCPMilnor1d>;

/* Groebner basis of critical pairs */
class CPMilnors
{
private:
    int t_trunc_;                                                        /* Truncation degree */
    CPMilnor2d gb_;                                                      /* `pairs_[j]` is the set of pairs (i, j) with given j */
    std::map<int, CPMilnor2d> buffer_min_pairs_;                         /* `buffer_min_pairs_[t]` To generate minimal pairs to compute Sij */
    std::map<int, std::unordered_set<uint64_t>> buffer_redundent_pairs_; /* Used to minimize `buffer_min_pairs_` */
    std::map<int, CPMilnor1d> buffer_singles_;                           /* For computing Sj. `buffer_singles_` stores indices of singles_ */

public:
    CPMilnors(int t_trunc) : t_trunc_(t_trunc) {}

    int t_trunc() const
    {
        return t_trunc_;
    }
    bool empty_min_pairs_for_gb(int t) const
    {
        return buffer_min_pairs_.find(t) == buffer_min_pairs_.end();
    }
    CPMilnor1d pairs_for_gb(int t)
    {
        CPMilnor1d result;
        if (!buffer_singles_.empty() && buffer_singles_.begin()->first == t) {
            std::swap(result, buffer_singles_.begin()->second);
            buffer_singles_.erase(t);
        }
        if (!buffer_min_pairs_.empty() && buffer_min_pairs_.begin()->first == t) {
            auto& b_min_pairs_t = buffer_min_pairs_.begin()->second;
            for (size_t j = 0; j < b_min_pairs_t.size(); ++j)
                for (auto& pair : b_min_pairs_t[j])
                    if (pair.i2 != -1)
                        result.push_back(std::move(pair));
            buffer_min_pairs_.erase(buffer_min_pairs_.begin());
        }

        return result;
    }

    /* Minimize `buffer_min_pairs_[t]` and maintain `pairs_` */
    void Minimize(const MMod1d& leads, int t);

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const MMod1d& leads, MMod mon, int t_v);

    void init(const MMod1d& leads, const array& basis_degrees);
};
using CPMilnors1d = std::vector<CPMilnors>;

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

// struct DataMRes
//{
//     Mod x1, x2, x3;
//     DataMRes() {}
//     DataMRes(Mod x1_, Mod x2_, Mod x3_) : x1(std::move(x1_)), x2(std::move(x2_)), x3(std::move(x3_)) {}
//     DataMRes operator+(const DataMRes& rhs) const
//     {
//         return DataMRes(x1 + rhs.x1, x2 + rhs.x2, x3 + rhs.x3);
//     }
//     DataMRes& operator+=(const DataMRes& rhs)
//     {
//         x1 += rhs.x1;
//         x2 += rhs.x2;
//         x3 += rhs.x3;
//         return *this;
//     }
// };

using DataMRes1d = std::vector<DataMRes>;
using DataMRes2d = std::vector<DataMRes1d>;

class GroebnerMRes
{
private:
    using TypeIndices = std::vector<std::unordered_map<uint32_t, array>>;

private:
    int t_trunc_;

    DataMRes2d gb_;
    MMod2d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    CPMilnors1d pairs_; /* Groebner basis of critical pairs */

    array2d basis_degrees_; /* `basis_degrees[s][i]` is the degree of v_{s,i} */

public:
    GroebnerMRes(int t_trunc, array2d basis_degrees) : t_trunc_(t_trunc), basis_degrees_(std::move(basis_degrees)) {}

    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    GroebnerMRes(int t_trunc, DataMRes2d data, array2d basis_degrees) : t_trunc_(t_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees))
    {
        if (basis_degrees_.empty())
            basis_degrees_.push_back({0});
        if (basis_degrees_[0].empty())
            basis_degrees_[0].push_back(0);

        leads_.resize(gb_.size());
        indices_.resize(gb_.size());

        for (size_t s = 0; s < gb_.size(); ++s) {
            for (int j = 0; j < (int)gb_[s].size(); ++j) {
                leads_[s].push_back(gb_[s][j].x1.GetLead());
                indices_[s][Key(gb_[s][j].x1.GetLead())].push_back(j);
            }
        }

        for (size_t s = 0; s < gb_.size(); ++s) {
            pairs_.push_back(CPMilnors(t_trunc_));
            pairs_.back().init(leads_[s], basis_degrees_[s]);
        }
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
            if (gb_[s][k].x2.GetLead().v_raw() == mon.v_raw() && divisibleLF(gb_[s][k].x2.GetLead(), mon))
                return (int)k;
        return -1;
    }

public: /* Getters and Setters */
    int t_trunc() const
    {
        return t_trunc_;
    }
    int t_begin() const
    {
        int t = 0;
        for (size_t s = 0; s < gb_.size(); ++s) {
            if (!gb_[s].empty()) {
                auto lead = gb_[s].back().x1.GetLead();
                int t1 = lead.deg_m() + basis_degrees_[s][lead.v()];
                if (t1 > t)
                    t = t1;
            }
        }
        return t;
    }
    /* This function will erase `gb_pairs_.buffer_min_pairs[t]` */
    CPMilnor1d pairs(size_t s, int t)
    {
        return pairs_[s].pairs_for_gb(t);
    }
    const array2d& basis_degs() const
    {
        return basis_degrees_;
    }
    auto gb_size() const
    {
        return gb_.size();
    }
    auto& operator[](size_t index) const
    {
        return gb_[index];
    }
    void resize_gb(size_t s)
    {
        if (gb_.size() < s) {
            gb_.resize(s);
            leads_.resize(s);
            indices_.resize(s);
            while (pairs_.size() < s)
                pairs_.push_back(CPMilnors(t_trunc_));
        }
    }
    void push_back(DataMRes g, size_t s)
    {
        MMod mv = g.x1.GetLead();
        pairs_[s].AddToBuffers(leads_[s], mv, basis_degrees_[s][mv.v()]);

        leads_[s].push_back(mv);
        indices_[s][Key(mv)].push_back((int)gb_[s].size());
        gb_[s].push_back(std::move(g));
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
    void MinimizePairs(size_t s, int t)
    {
        pairs_[s].Minimize(leads_[s], t);
    }

    const auto& data() const
    {
        return gb_;
    }

public:
    DataMRes Reduce(CPMilnor& p, int s) const
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
            result.x1.iaddmul(p.m1, gb_[s][p.i1].x1, tmp_a, tmp_m1, tmp_m2).iaddmul(p.m2, gb_[s][p.i2].x1, tmp_a, tmp_m1, tmp_m2);
            result.x2.iaddmul(p.m1, gb_[s][p.i1].x2, tmp_a, tmp_m1, tmp_m2).iaddmul(p.m2, gb_[s][p.i2].x2, tmp_a, tmp_m1, tmp_m2);
        }
        else {
            result.x1.iaddmul(p.m2, gb_[s][p.i2].x1, tmp_a, tmp_m1, tmp_m2);
            result.x2.iaddmul(p.m2, gb_[s][p.i2].x2, tmp_a, tmp_m1, tmp_m2);
        }

        size_t index;
        index = 0;
        while (index < result.x1.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x1.data[index], s);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x1.data[index], gb_[s][gb_index].x1.data[0]);
                result.x1.iaddmul(m, gb_[s][gb_index].x1, tmp_a, tmp_m1, tmp_m2);
                result.x2.iaddmul(m, gb_[s][gb_index].x2, tmp_a, tmp_m1, tmp_m2);
            }
            else
                ++index;
        }
        index = 0;
        while (index < result.x2.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x2.data[index], s + 1);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x2.data[index], gb_[sp1][gb_index].x1.data[0]);
                result.x2.iaddmul(m, gb_[sp1][gb_index].x1, tmp_a, tmp_m1, tmp_m2);
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
 * Comsume relations from 'rels` and `gb.pairs_` in degree `<= deg`.
 *
 * return the dimension of the calculated range for debugging.
 */
size_t AddRelsMRes(GroebnerMRes& gb, const Mod1d& rels, int deg);

}  // namespace steenrod

#endif