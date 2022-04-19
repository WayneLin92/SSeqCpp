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

/*
 * This class allows exponents up to 255
*/
struct MMilnorE
{
    std::array<unsigned char, MMILNOR_INDEX_NUM + 1> data = {};
    MMilnorE mulLF(MMilnor m) const
    {
        MMilnorE result = *this;
        result.data[0] += m.weight();
        for (int i : m)
            ++result.data[size_t(i + 1)];
        return result;
    }
    MMilnorE& imul(MMilnor m)
    {
        data[0] += m.weight();
        for (int i : m)
            ++data[size_t(i + 1)];
        return *this;
    }
    bool operator<(const MMilnorE& rhs) const
    {
        return std::memcmp(data.data(), rhs.data.data(), sizeof(data)) < 0;
    }
    bool operator==(const MMilnorE& rhs) const
    {
        return std::memcmp(data.data(), rhs.data.data(), sizeof(data)) == 0;
    }
    std::string Str() const;
};
using MMilnorE1d = std::vector<MMilnorE>;
using MMilnorE2d = std::vector<MMilnorE1d>;

inline std::ostream& operator<<(std::ostream& sout, const MMilnorE& x)
{
    return sout << x.Str();
}

inline MMilnorE mulE(MMilnor m1, MMilnor m2)
{
    return MMilnorE().imul(m1).imul(m2);
}

struct CmpMMod
{
    const MMilnorE1d& mDiffLF_;
    const array& ind_indices_;
    CmpMMod(const MMilnorE1d& mDiffLF, const array& ind_indices) : mDiffLF_(mDiffLF), ind_indices_(ind_indices) {}
    bool operator()(MMod mv1, MMod mv2)  // TODO: if v1 == v2 then no need to compute m1, m2
    {
        int v1 = mv1.v(), v2 = mv2.v();
        MMilnorE m1 = mDiffLF_[v1].mulLF(mv1.m());
        MMilnorE m2 = mDiffLF_[v2].mulLF(mv2.m());
        if (m1 < m2)
            return true;
        if (m2 < m1)
            return false;
        if (ind_indices_[v2] < ind_indices_[v1])
            return true;
        return false;
    }
};

struct DataMRes
{
    Mod x1, x2, x1LF, x2LF;
    DataMRes() {}
    DataMRes(Mod x1_, Mod x2_, Mod x1LF_, Mod x2LF_) : x1(std::move(x1_)), x2(std::move(x2_)), x1LF(std::move(x1LF_)), x2LF(std::move(x2LF_)) {}
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

    //Mod2d diff_;              /* `diff_[s][i]` equals d(v_{s,i}) */
    array2d basis_degrees_;   /* `basis_degrees[s][i]` is the degree of v_{s,i} */
    MMilnorE2d mDiffLF_;      /* `mDiffLF_[s][i]` is the iterated leading monomial of v_{s,i} */
    array2d vDiffLF_;         /* `vDiffLF_[s][i]` is the index of d(v_{s,i}) */
    array2d ordered_indices_; /* `ordered_indices_[s]` is the set of sorted indices depending on `vDiffLF_` */
    array2d ind_indices_;     /* `ind_indices_[s][i]` is the position of i in `ordered_indices_[s]` */

public:
    GroebnerMRes(int deg_trunc, array2d basis_degrees) : gb_pairs_(deg_trunc), basis_degrees_(std::move(basis_degrees)) {}

    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    GroebnerMRes(int deg_trunc, DataMRes2d data, array2d basis_degrees) : gb_pairs_(deg_trunc), data_(std::move(data)), basis_degrees_(std::move(basis_degrees))
    {
        if (basis_degrees_.empty()) {
            basis_degrees_.push_back({});
            mDiffLF_.push_back({});
            vDiffLF_.push_back({});
            ordered_indices_.push_back({});
            ind_indices_.push_back({});
        }
        if (basis_degrees_[0].empty()) {
            basis_degrees_[0].push_back(0);
            mDiffLF_[0].push_back(MMilnorE());
            vDiffLF_[0].push_back(0);
            ordered_indices_[0].push_back(0);
            ind_indices_[0].push_back(0);
        }

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
                if (divisibleLF(leads_[s][k].m(), mon.m()))
                    return k;
        return -1;
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeadingX2(MMod mon, int s) const
    {
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
            size_t sp1 = size_t(s + 1); /* s plus 1 */
            data_.resize(sp1);
            leads_.resize(sp1); 
            indices_.resize(sp1);
            gb_pairs_.resize_pairs(sp1);
        }
    }
    void push_back(DataMRes g, int s, int t)
    {
        g.x1LF = GetModLF(g.x1, s);
        g.x2LF = GetModLF(g.x2, s + 1);
        MMod mv = g.x1.GetLead();
        gb_pairs_.AddToBuffers(leads_[s], mv, basis_degrees_[s][mv.v()], s);  // TODO: modify the counterpart in algebras/groebner.h

        leads_[s].push_back(mv);
        indices_[s][Key(mv)].push_back((int)data_[s].size());

        data_[s].push_back(std::move(g));
    }

    /* Add x2 + v_{s+1,i} */
    void push_back_kernel(Mod x2, DataMRes1d& rels, int s, int t)
    {
        //if (diff_.size() <= s)
        //    diff_.resize(size_t(s + 1));
        //diff_[s].push_back(x2);
        size_t sp1 = size_t(s + 1); /* s plus 1 */
        if (basis_degrees_.size() <= sp1) {
            size_t sp2 = size_t(s + 2);
            basis_degrees_.resize(sp2);
            mDiffLF_.resize(sp2);
            vDiffLF_.resize(sp2);
            ordered_indices_.resize(sp2);
            ind_indices_.resize(sp2);
        }
        basis_degrees_[sp1].push_back(t);
        int v = x2.GetLead().v();
        mDiffLF_[sp1].push_back(mDiffLF_[s][v].mulLF(x2.GetLead().m()));
        vDiffLF_[sp1].push_back(v);

        int n = (int)ordered_indices_[sp1].size();
        auto p = std::upper_bound(ordered_indices_[sp1].begin(), ordered_indices_[sp1].end(), n, [&](int i, int j) { return ind_indices_[s][vDiffLF_[sp1][i]] < ind_indices_[s][vDiffLF_[sp1][j]]; });
        size_t index = p - ordered_indices_[sp1].begin();
        ordered_indices_[sp1].insert(p, n);

        ind_indices_[sp1].resize(ordered_indices_[sp1].size());
        for (size_t i = index; i < ordered_indices_[sp1].size(); ++i)
            ind_indices_[sp1][ordered_indices_[sp1][i]] = (int)i;


        Mod x3 = MMod(MMilnor(0), (int)basis_degrees_[sp1].size() - 1);

        DataMRes g = DataMRes{std::move(x2), std::move(x3), Mod(), Mod()};
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
    CmpMMod GetCmpMMod(int s) const
    {
        return CmpMMod(mDiffLF_[s], ind_indices_[s]);
    }
    /* Return -1 for <, 1 for > and 0 for == */
    int IntCmpDiff(MMilnor m1, MMod x1, MMilnor m2, MMod x2, int s) const
    {
        int v1 = x1.v(), v2 = x2.v();
        MMilnorE me1 = mDiffLF_[s][v1].mulLF(x1.m()).mulLF(m1);
        MMilnorE me2 = mDiffLF_[s][v2].mulLF(x2.m()).mulLF(m2);
        if (me1 < me2)
            return -1;
        if (me2 < me1)
            return 1;
        size_t sm1 = size_t(s - 1);
        if (ind_indices_[sm1][vDiffLF_[s][v1]] > ind_indices_[sm1][vDiffLF_[s][v2]])
            return -1;
        if (ind_indices_[sm1][vDiffLF_[s][v1]] < ind_indices_[sm1][vDiffLF_[s][v2]])
            return 1;
        return 0;
    }
    Mod GetModLF(const Mod& x, int s) const
    {
        Mod result(x.GetLead());
        int v = result.GetLead().v();
        int dv = vDiffLF_[s][v];
        MMilnorE m_lead = mDiffLF_[s][v].mulLF(result.GetLead().m());
        for (size_t i = 1; i < x.data.size(); ++i) {
            int vi = x.data[i].v();
            if (vDiffLF_[s][vi] == dv && mDiffLF_[s][vi].mulLF(x.data[i].m()) == m_lead)
                result.data.push_back(x.data[i]);
            else
                break;
        }
        return result;
    }

public:
    DataMRes Reduce(CriPairMRes& p, int s) const
    {
        DataMRes result;
        size_t sp1 = size_t(s + 1);
        auto cmp_s = CmpMMod(mDiffLF_[s], ind_indices_[s]);
        auto cmp_sp1 = CmpMMod(mDiffLF_[sp1], ind_indices_[sp1]);

        if (p.i1 >= 0) {
            result.x1 = mulMod(p.m1, data_[s][p.i1].x1, cmp_s).add(mulMod(p.m2, data_[s][p.i2].x1, cmp_s), cmp_s);
            result.x2 = mulMod(p.m1, data_[s][p.i1].x2, cmp_sp1).add(mulMod(p.m2, data_[s][p.i2].x2, cmp_sp1), cmp_sp1);
        }
        else {
            result.x1 = mulMod(p.m2, data_[s][p.i2].x1, cmp_s);
            result.x2 = mulMod(p.m2, data_[s][p.i2].x2, cmp_sp1);
        }

        size_t index;
        index = 0;
        while (index < result.x1.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x1.data[index], s);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x1.data[index].m(), data_[s][gb_index].x1.data[0].m());
                result.x1.iadd(mulMod(m, data_[s][gb_index].x1, cmp_s), cmp_s);
                result.x2.iadd(mulMod(m, data_[s][gb_index].x2, cmp_sp1), cmp_sp1);
            }
            else
                ++index;
        }
        index = 0;
        while (index < result.x2.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.x2.data[index], s + 1);
            if (gb_index != -1) {
                MMilnor m = divLF(result.x2.data[index].m(), data_[sp1][gb_index].x1.data[0].m());
                result.x2.iadd(mulMod(m, data_[sp1][gb_index].x1, cmp_sp1), cmp_sp1);
            }
            else
                ++index;
        }

        return result;
    }

    /* Reduce x2LF of `p` by `data_[s + 1]` */
    Mod ReduceX2LF(CriPairMRes& p, int s) const
    {
        Mod result;
        size_t sp1 = size_t(s + 1);

        if (p.i1 >= 0) {
            int iCmp = IntCmpDiff(p.m1, data_[s][p.i1].x2LF.GetLead(), p.m2, data_[s][p.i2].x2LF.GetLead(), s + 1);

            if (iCmp < 0)
                result = mulLF(p.m1, data_[s][p.i1].x2LF);
            else if (iCmp > 0)
                result = mulLF(p.m2, data_[s][p.i2].x2LF);
            else
                result = mulLF(p.m1, data_[s][p.i1].x2LF).addLF(mulLF(p.m2, data_[s][p.i2].x2LF));
        }
        else
            result = mulLF(p.m2, data_[s][p.i2].x2LF);

        size_t index = 0;
        while (index < result.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.data[index], s + 1);
            if (gb_index != -1) {
                MMilnor m = divLF(result.data[index].m(), data_[sp1][gb_index].x1LF.GetLead().m());
                result.iaddLF(mulLF(m, data_[sp1][gb_index].x1LF));
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
 */
void AddRelsMRes(GroebnerMRes& gb, const Mod1d& rels, int deg);

}  // namespace steenrod

#endif