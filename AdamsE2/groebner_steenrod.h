/** \file groebner_steenrod.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_STEENROD_H
#define GROEBNER_STEENROD_H

#include "algebras/steenrod.h"
#include <execution>
#include <map>
#include <unordered_map >
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

inline void print(const ModCpt& x)  // For debugging
{
    if (x.data.empty()) {
        std::cout << '0';
    }
    for (auto pm = x.data.begin(); pm != x.data.end(); ++pm) {
        if (pm != x.data.begin())
            std::cout << '+';
        for (int i : pm->m())
            std::cout << "P_{" << MMILNOR_GEN_I[i] << MMILNOR_GEN_J[i] << '}';
        std::cout << "v_{" << pm->v() << '}';
    }
}

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
    void set_deg_trunc(int d_trunc)
    {
        if (d_trunc <= deg_trunc_) {
            deg_trunc_ = d_trunc;  ////
        }
        else {
            throw MyException(0x6e1832b6U, "NotImplemented");  ////
        }
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
                    result.push_back(std::move(pair));
            buffer_min_pairs_.erase(st);
        }

        return result;
    }

    /* Minimize `buffer_min_pairs_[t]` and maintain `pairs_` */
    void Minimize(const MModCpt1d& leads, AdamsDeg st);

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const MModCpt1d& leads, MModCpt mon, int d_mon_base, int s);

    void init(const MModCpt2d& leads, const array2d& basis_degrees);
};

struct DataMRes
{
    ModCpt x1, x2, x2LF;
    DataMRes operator+(const DataMRes& rhs) const
    {
        return DataMRes{x1 + rhs.x1, x2 + rhs.x2};
    }
    DataMRes& operator+=(const DataMRes& rhs)
    {
        x1 += rhs.x1;
        x2 += rhs.x2;
        return *this;
    }
};
inline DataMRes operator*(const Milnor& a, const DataMRes& x)
{
    return DataMRes{a * x.x1, a * x.x2};
}
using DataMRes1d = std::vector<DataMRes>;
using DataMRes2d = std::vector<DataMRes1d>;

class GroebnerMRes
{
private:
    using TypeIndices = std::vector<std::unordered_map<uint32_t, array>>;

private:
    GbCriPairsMRes gb_pairs_; /* Groebner basis of critical pairs */
    DataMRes2d data_;
    array2d basis_degrees_; /* `basis_degrees[s][i]` is the degree of v_{s,i} */
    ModCpt2d kernel_;       /* kernel[s] is the generators of Kernel of $F_{s}\to F_{s-1}$ */

    /* Caches */
    MModCpt2d leads_;     /* Leading monomials */
    ModCpt3d dataX2LF_;   /* dataLF_[s][t] is the vector space of leading forms of x2 */
    TypeIndices indices_; /* Cache for fast divisibility test */

public:
    GroebnerMRes(int deg_trunc, array2d basis_degrees) : gb_pairs_(deg_trunc), basis_degrees_(std::move(basis_degrees)) {}

    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    GroebnerMRes(int deg_trunc, DataMRes2d data, array2d basis_degrees) : gb_pairs_(deg_trunc), data_(std::move(data)), basis_degrees_(std::move(basis_degrees))
    {
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
    static uint32_t Key(MModCpt lead)
    {
        return uint32_t(lead.v());
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(MModCpt mon, int s) const
    {
        auto key = uint32_t(mon.v());
        auto p = indices_[s].find(key);
        if (p != indices_[s].end())
            for (int k : p->second)
                if (divisibleLF(leads_[s][k].m(), mon.m()))
                    return k;
        return -1;
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeadingX2(MModCpt mon, int s) const
    {
        /*auto key = uint32_t(mon.v());
        auto p = indices_[s].find(key);
        if (p != indices_[s].end())
            for (int k : p->second)
                if (divisibleLF(leads_[s][k].m(), mon.m()))
                    return k;*/
        for (size_t k = 0; k < leads_[s].size(); ++k)
            if (data_[s][k].x2.GetLead().v() == mon.v() && divisibleLF(data_[s][k].x2.GetLead().m(), mon.m()))
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
    void set_deg_trunc(int d_trunc)
    {
        gb_pairs_.set_deg_trunc(d_trunc);
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
            data_.resize(size_t(s + 1));
            dataX2LF_.resize(size_t(s + 1));
            leads_.resize(size_t(s + 1));
            indices_.resize(size_t(s + 1));
            gb_pairs_.resize_pairs(size_t(s + 1));
        }
    }
    void resize_dataX2LF(int s, int t)
    {
        dataX2LF_.resize(size_t(t + 1));
    }
    void push_back(DataMRes g, int s, int t)
    {
        MModCpt mv = g.x1.GetLead();
        gb_pairs_.AddToBuffers(leads_[s], mv, basis_degrees_[s][mv.v()], s);  // TODO: modify the counterpart in algebras/groebner.h

        leads_[s].push_back(mv);
        indices_[s][Key(mv)].push_back((int)data_[s].size());

        dataX2LF_[s].resize(size_t(t + 1));
        dataX2LF_[s][t].push_back(g.x2LF);
        data_[s].push_back(std::move(g));
    }
    /* Add x2 + v_{s+1,k} */
    void push_back_kernel(ModCpt x2, DataMRes1d& rels, int s, int t)  // TODO: s, t
    {
        if (kernel_.size() <= s)
            kernel_.resize(size_t(s + 1));
        if (basis_degrees_.size() <= size_t(s + 1))
            basis_degrees_.resize(size_t(s + 2));
        kernel_[s].push_back(x2);
        basis_degrees_[size_t(s + 1)].push_back(t);

        ModCpt x3 = MModCpt(MMilnor(0), (int)basis_degrees_[size_t(s + 1)].size() - 1);
        DataMRes g = DataMRes{std::move(x2), x3, x3};
        rels.push_back(g);
        push_back(std::move(g), s, t);
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
    const auto& kernel() const
    {
        return kernel_;
    }

public:
    DataMRes Reduce(CriPairMRes& p, int s) const
    {
        DataMRes result;

        if (p.i1 >= 0)
            result = Milnor(p.m1) * data_[s][p.i1] + Milnor(p.m2) * data_[s][p.i2];
        else
            result = Milnor(p.m2) * data_[s][p.i2];

        size_t index;

        index = 0;
        while (index < result.x1.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x1.data[index], s);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x1.data[index].m(), data_[s][gb_index].x1.data[0].m());
                result += Milnor(m) * data_[s][gb_index];
            }
            else
                ++index;
        }
        index = 0;
        while (index < result.x2.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x2.data[index], s + 1);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x2.data[index].m(), data_[size_t(s + 1)][gb_index].x1.data[0].m());
                result.x2 += Milnor(m) * data_[size_t(s + 1)][gb_index].x1;
            }
            else
                ++index;
        }

        return result;
    }

    /* Reduce x2 by data_[s] */
    ModCpt ReduceX2(ModCpt x2, int s) const
    {
        size_t index = 0;
        while (index < x2.data.size()) {
            int gb_index = IndexOfDivisibleLeading(x2.data[index], s);
            if (gb_index != -1) {
                MMilnor m = divLF(x2.data[index].m(), data_[s][gb_index].x1.data[0].m());
                x2 += Milnor(m) * data_[s][gb_index].x1;
            }
            else
                ++index;
        }
        return x2;
    }

    /* Reduce x2 by data_[s] */
    ModCpt ReduceX2LF(CriPairMRes& p, int s, int t) const
    {
        ModCpt result;

        if (p.i1 >= 0)
            result = mulLF(p.m1, data_[s][p.i1].x2LF) + mulLF(p.m2, data_[s][p.i2].x2LF);
        else
            result = mulLF(p.m2, data_[s][p.i2].x2LF);

        /*for (auto& x : data_[size_t(s + 1)])
            if (std::binary_search(result.data.begin(), result.data.end(), x.x1.GetLead()))
                result += x.x1;*/
        return result;
    }

public:
    static GroebnerMRes load(const std::string& filename, int t_trunc);
};

/**
 * Comsume relations from 'rels` and `gb.gb_pairs_` in degree `<= deg`
 * `min_gb` stores the minimal generating set of gb.
 */
void AddRelsMRes(GroebnerMRes& gb, const ModCpt1d& rels, int deg);
void generate_basis(const std::string& filename, int t_trunc);

}  // namespace steenrod

#endif