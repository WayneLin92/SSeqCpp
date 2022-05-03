/** \file groebner_steenrod.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_STEENROD_H
#define GROEBNER_STEENROD_H

#include "algebras/benchmark.h"
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
    friend class GroebnerMRes;

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
    CPMilnor1d cpairs_for_gb(int t)
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
                        result.push_back(pair);
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

struct Filtr
{
    uint64_t data;
    Filtr() : data(~0) {}
    Filtr(MMod m) : data(m.v_raw() | m.w_may()) {}
    Filtr(uint64_t data_) : data(data_) {}
    bool operator<(Filtr rhs) const
    {
        return data < rhs.data;
    }
    bool operator==(Filtr rhs) const
    {
        return data == rhs.data;
    }
    Filtr operator+(uint64_t w_may) const
    {
        return Filtr(data + w_may);
    }
};

struct DataMRes
{
    Mod x1, x2, x2m;
    Filtr fil;
    DataMRes() {}
    DataMRes(Mod x1_, Mod x2_, Mod x2m_) : x1(std::move(x1_)), x2(std::move(x2_)), x2m(std::move(x2m_))
    {
        if (x1)
            fil = Filtr(x1.GetLead());
    }
    DataMRes& operator+=(const DataMRes& rhs)
    {
#ifndef NDEBUG
        if (!rhs.valid_x2m())
            throw MyException(0x2ae8baa3U, "Add only when rhs.x2m is valid");
#endif
        if (valid_x2m() && fil == rhs.fil)
            x2m += rhs.x2m;
        x1 += rhs.x1;
        x2 += rhs.x2;
        return *this;
    }
    bool valid_x2m() const
    {
        if (x1)
            return fil == Filtr(x1.GetLead());
        if (x2m)
            return false;
        return true;
    }
};

using DataMRes1d = std::vector<DataMRes>;
using DataMRes2d = std::vector<DataMRes1d>;

class GroebnerMRes
{
private:
    using TypeIndices = std::vector<std::unordered_map<uint32_t, array>>;

private:
    int t_trunc_;

    CPMilnors1d cpairs_; /* Groebner basis of critical pairs */

    DataMRes2d gb_;
    MMod2d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    array2d basis_degrees_; /* `basis_degrees[s][i]` is the degree of v_{s,i} */

    CPMilnors1d cpairs_x2m_; /* Groebner basis of critical pairs */

    Mod2d gb_x2m_;
    MMod2d leads_x2m_;        /* Leading monomials */
    TypeIndices indices_x2m_; /* Cache for fast divisibility test */

    array2d basis_degrees_x2m_; /* `basis_degrees_x2m_[s][i]` is the degree of w_{s,i} */

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
            cpairs_.push_back(CPMilnors(t_trunc_));
            cpairs_.back().init(leads_[s], basis_degrees_[s]);
        }

        ////x2m
    }

private:
    static uint32_t Key(MMod lead)
    {
        return uint32_t(lead.v());
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(MMod mon, size_t s) const
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
    int IndexOfDivisibleLeadingX2m(MMod mon, size_t s) const
    {
        auto key = uint32_t(mon.v());
        auto p = indices_x2m_[s].find(key);
        if (p != indices_x2m_[s].end())
            for (int k : p->second)
                if (divisibleLF(leads_x2m_[s][k], mon))
                    return k;
        return -1;
    }

    Mod ReduceByGbX2m(Mod x2m, size_t s) const;
    Mod ReduceByGbX2m(const CPMilnor& p, size_t s) const;
    void AddRelsX2m(size_t s, int t);

public:
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

    CPMilnor1d cpairs(size_t s, int t)
    {
        AddRelsX2m(s, t);
        cpairs_[s].Minimize(leads_[s], t);
        CPMilnor1d cps = cpairs_[s].cpairs_for_gb(t);
        std::vector<Filtr> fils(cps.size());
        for (size_t i = 0; i < cps.size(); ++i)
            fils[i] = gb_[s][cps[i].i2].fil + cps[i].m2.w_may();
        auto indices = ut::size_t_range(cps.size());
        std::sort(indices.begin(), indices.end(), [&fils](size_t i, size_t j) { return fils[i] < fils[j]; });
        CPMilnor1d result;
        result.reserve(cps.size());
        for (size_t i = 0; i < cps.size(); ++i)
            result.push_back(cps[indices[i]]);

        for (auto& cp : result) {
            if (!ReduceX2m(cp, s)) {
                bench::Counter(4);
                cp.i2 = -1;
            }
        }
        ut::RemoveIf(result, [](const CPMilnor& cp) { return cp.i2 == -1; });
        return result;
    }

    const array2d& basis_degs() const  //// para s
    {
        return basis_degrees_;
    }
    const array& basis_degs_x2m(size_t s) const
    {
        return basis_degrees_x2m_[s];
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
            while (cpairs_.size() < s)
                cpairs_.push_back(CPMilnors(t_trunc_));
        }
    }
    void resize_gb_x2m(size_t s)
    {
        if (gb_x2m_.size() < s) {
            gb_x2m_.resize(s);
            leads_x2m_.resize(s);
            indices_x2m_.resize(s);
            while (cpairs_x2m_.size() < s)
                cpairs_x2m_.push_back(CPMilnors(t_trunc_));
        }
    }
    Mod new_gen(size_t s, int t)
    {
        if (basis_degrees_.size() < s + 1)  ////
            basis_degrees_.resize(s + 1);
        basis_degrees_[s].push_back(t);
        return MMod(MMilnor(), basis_degrees_[s].size() - 1);
    }
    Mod new_gen_x2m(size_t s, int t)
    {
        if (basis_degrees_x2m_.size() < s + 1)  ////
            basis_degrees_x2m_.resize(s + 1);
        basis_degrees_x2m_[s].push_back(t);
        return MMod(MMilnor(), basis_degrees_x2m_[s].size() - 1);
    }
    void push_back(DataMRes g, size_t s)
    {
        MMod m = g.x1.GetLead();
        cpairs_[s].AddToBuffers(leads_[s], m, basis_degrees_[s][m.v()]);

        leads_[s].push_back(m);
        indices_[s][Key(m)].push_back((int)gb_[s].size());
        gb_[s].push_back(std::move(g));
    }
    void push_back_x2m(Mod g, size_t s)
    {
        MMod m = g.GetLead();
        cpairs_x2m_[s].AddToBuffers(leads_x2m_[s], m, basis_degrees_x2m_[s][m.v()]);

        leads_x2m_[s].push_back(m);
        indices_x2m_[s][Key(m)].push_back((int)gb_x2m_[s].size());
        gb_x2m_[s].push_back(std::move(g));
    }

    /* Add x2 + v_{s+1,i} */
    void push_back_kernel(Mod x2, DataMRes1d& rels, size_t s, int t)
    {
        DataMRes g = DataMRes(std::move(x2), new_gen(s + 1, t), new_gen_x2m(s, t));
        rels.push_back(g);
        push_back(std::move(g), s);
    }

    const auto& data() const
    {
        return gb_;
    }
    const auto& data_x2m() const
    {
        return gb_x2m_;
    }

public:
    DataMRes Reduce(const CPMilnor& cp, size_t s) const;
    Mod ReduceX2(const CPMilnor& cp, int s) const;  ////
    Mod ReduceX2m(const CPMilnor& cp, size_t s) const;

public:
    static GroebnerMRes load(const std::string& filename, int t_trunc);
};

/**
 * Comsume relations from 'rels` and `gb.cpairs_` in degree `<= deg`.
 *
 * return the dimension of the calculated range for debugging.
 */
size_t AddRelsMRes(GroebnerMRes& gb, const Mod1d& rels, int deg);

}  // namespace steenrod

#endif