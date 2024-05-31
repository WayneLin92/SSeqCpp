/** \file groebner_steenrod.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_STEENROD_H
#define GROEBNER_STEENROD_H

#include "benchmark.h"
#include "steenrod.h"
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace steenrod {

class Groebner;

struct CriMilnor
{
    int i1 = -1, i2 = -1;
    MMilnor m1, m2;

    /* Compute the pair for two leading monomials. */
    CriMilnor() = default;
    static void SetFromLM(CriMilnor& result, MMilnor lead1, MMilnor lead2, int i, int j)
    {
        MMilnor gcd = gcdLF(lead1, lead2);
        result.m1 = divLF(lead2, gcd);
        result.m2 = divLF(lead1, gcd);
        result.i1 = i;
        result.i2 = j;
    }
    static CriMilnor Single(MMilnor m2, int j)
    {
        CriMilnor result;
        result.m2 = m2;
        result.i1 = -1;
        result.i2 = j;
        return result;
    }
    void Sij(const Groebner& gb, Mod& result, Mod& tmp1, Mod& tmp2) const;
};
using CriMilnor1d = std::vector<CriMilnor>;
using CriMilnor2d = std::vector<CriMilnor1d>;
using CriMilnor3d = std::vector<CriMilnor2d>;
using PtCriMilnor1d = std::vector<CriMilnor*>;
using PtCriMilnor2d = std::vector<PtCriMilnor1d>;

/* Groebner basis of critical pairs */
class CriMilnors
{
private:
    int t_trunc_;                                                        /* Truncation degree */
    CriMilnor2d gb_;                                                     /* `pairs_[j]` is the set of pairs (i, j) with given j */
    std::map<int, CriMilnor2d> buffer_min_pairs_;                        /* `buffer_min_pairs_[t]` To generate minimal pairs to compute Sij */
    std::map<int, std::unordered_set<uint64_t>> buffer_redundent_pairs_; /* Used to minimize `buffer_min_pairs_` */
    std::map<int, CriMilnor1d> buffer_singles_;                          /* For computing Sj. `buffer_singles_` stores indices of singles_ */

public:
    CriMilnors(int t_trunc) : t_trunc_(t_trunc) {}

    int t_trunc() const
    {
        return t_trunc_;
    }
    /* Return both critical pairs and critical singles */
    CriMilnor1d Criticals(int t);
    bool empty() const
    {
        return buffer_min_pairs_.empty() && buffer_singles_.empty();
    }

    void Reset()
    {
        gb_.clear();
        buffer_min_pairs_.clear();
        buffer_redundent_pairs_.clear();
        buffer_singles_.clear();
    }

    /* Minimize `buffer_min_pairs_[t]` and maintain `pairs_` */
    void Minimize(const MMod1d& leads, int t);

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const MMod1d& leads, MMod mon, int t_v);

    void init(const MMod1d& leads, const int1d& basis_degrees, int t_min_buffer);
};
using CriMilnors1d = std::vector<CriMilnors>;

/********************************************************
 *                    class Groebner
 ********************************************************/

class Groebner
{
    using TypeIndices = std::unordered_map<uint64_t, int1d>;

private:
    int d_trunc_;

    CriMilnors criticals_; /* Groebner basis of critical pairs */

    Mod1d gb_;
    MMod1d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int1d v_degs_; /* `basis_degrees[i]` is the degree of v_{s,i} */

public:
    Groebner(int d_trunc, Mod1d data, int1d v_degs);

public:
    int d_trunc() const
    {
        return d_trunc_;
    }

    auto& v_degs() const
    {
        return v_degs_;
    }

    void set_v_degs(int1d v_degs)
    {
        v_degs_ = std::move(v_degs);
    }

    Mod new_gen(int t)
    {
        v_degs_.push_back(t);
        return MMod(MMilnor(), v_degs_.size() - 1);
    }

    void push_back(Mod g)
    {
        MMod m = g.GetLead();
        criticals_.AddToBuffers(leads_, m, v_degs_[m.v()]);

        leads_.push_back(m);
        indices_[m.v_raw()].push_back((int)gb_.size());
        gb_.push_back(std::move(g));
    }

    void ResetRels()
    {
        criticals_.Reset();
        leads_.clear();
        indices_.clear();
        gb_.clear();
        v_degs_.clear();
    }

    auto& data() const
    {
        return gb_;
    }

    bool IsBasis(MMod m) const;

    CriMilnor1d Criticals(int t);
    Mod Reduce(const CriMilnor& cp) const;
    Mod Reduce(Mod x) const;

    void AddRels(const Mod1d& rels, int deg_max, int1d& min_rels);
    void AddRels(const Mod1d& rels, int deg_max);

    /* Assume that the v_degs is in increasing order.
     * Simplify generators and relations.
     * output: cells[i] is the new presentation of the old v_i.
     * min_rels: set of indices of indecomposable relations
     */
    void MinimizeOrderedGensRels(Mod1d& cells, int1d& min_rels);
};

}  // namespace steenrod

#endif