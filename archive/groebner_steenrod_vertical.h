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

namespace detail {
    /* Assume that f[i] is divisable by g[0].
     * Set m = LF(f[i]/g[0])
     * Replace f with f + m * g.
     */
    inline void Reduce(Mod& f, const Mod& g, const size_t index, MMilnor& m)
    {
        m = divLF(f.data[index].m, g.data[0].m);
        f += Milnor(m) * g;
    }

    inline void Reduce(Mod& f, const Mod& g, const size_t index)
    {
        MMilnor m = divLF(f.data[index].m, g.data[0].m);
        f += Milnor(m) * g;
    }
}  // namespace detail

struct CriticalPair
{
    int i1 = -1, i2 = -1;
    MMilnor m1, m2;
    Mod x;

    /* Compute the pair for two leading monomials. */
    CriticalPair() = default;
    static void SetFromLM(CriticalPair& result, MMilnor lead1, MMilnor lead2, int i, int j)
    {
        MMilnor gcd = gcdLF(lead1, lead2);
        result.m1 = divLF(lead2, gcd);
        result.m2 = divLF(lead1, gcd);
        result.i1 = i;
        result.i2 = j;
    }
    static CriticalPair Single(MMilnor m2, int j)
    {
        CriticalPair result;
        result.i1 = -1;
        result.i2 = j;
        result.m2 = m2;
        return result;
    }
};
using CriticalPair1d = std::vector<CriticalPair>;
using PCriticalPair1d = std::vector<CriticalPair*>;
using CriticalPair2d = std::vector<CriticalPair1d>;

/* Groebner basis of critical pairs */
class GbCriPairs
{
private:
    int deg_trunc_;                                                      /* Truncation degree */
    CriticalPair2d pairs_;                                               /* `pairs_[j]` is the set of pairs (i, j) with given j */
    CriticalPair2d singles_;                                             /* `singles[d]` is the set of gen *  */
    array3d min_pairs_;                                                  /* Minimal generating set of `pairs_`. `min_pairs_` stores indices of `pairs` */
    std::map<int, CriticalPair2d> buffer_min_pairs_;                     /* To generate `min_pairs_` and for computing Sij */
    std::map<int, array2d> buffer_singles_;                              /* For computing Sj. `buffer_singles_` store indices of singles_ */
    std::map<int, std::unordered_set<uint64_t>> buffer_redundent_pairs_; /* Used to minimize `buffer_min_pairs_` */

public:
    GbCriPairs(int d_trunc) : deg_trunc_(d_trunc) {}
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
    bool empty_pairs_for_gb_d(int d) const
    {
        return buffer_min_pairs_.find(d) == buffer_min_pairs_.end() && buffer_singles_.find(d) == buffer_singles_.end();
    }
    bool empty_pairs_for_gb() const
    {
        return buffer_min_pairs_.empty() && buffer_singles_.empty();
    }
    PCriticalPair1d pairs_for_gb(int d)  // Warning: return pointers
    {
        PCriticalPair1d result;
        if ((size_t)d < min_pairs_.size())
            for (size_t j = 0; j < min_pairs_[d].size(); ++j)
                for (auto i : min_pairs_[d][j])
                    result.push_back(&pairs_[j][i]);
        if (buffer_singles_.find(d) != buffer_singles_.end())
            for (size_t j = 0; j < buffer_singles_[d].size(); ++j)
                for (auto i : buffer_singles_[d][j])
                    result.push_back(&singles_[j][i]);

        buffer_min_pairs_.erase(d);
        buffer_singles_.erase(d);
        return result;
    }
    Mod1d MinSyzOfGb() const
    {
        Mod1d result;
        for (auto& min_pairs_d : min_pairs_)
            for (size_t j = 0; j < min_pairs_d.size(); ++j)
                for (auto i : min_pairs_d[j])
                    result.push_back(pairs_[j][i].x);
        for (auto& singles_d : singles_)
            for (auto& s : singles_d)
                result.push_back(s.x);
        return result;
    }

    /* Minimize `buffer_min_pairs_[d]` and maintain `pairs_` */
    void AddAndMinimize(const MMod1d& leads, int d)
    {
        /* Add to the Groebner basis of critical pairs */
        if (buffer_min_pairs_.find(d) == buffer_min_pairs_.end())
            return;
        auto& b_min_pairs_d = buffer_min_pairs_.at(d);
        if (pairs_.size() < b_min_pairs_d.size())
            pairs_.resize(b_min_pairs_d.size());
        array old_sizes(b_min_pairs_d.size());
        for (int j = 0; j < (int)b_min_pairs_d.size(); ++j) {
            old_sizes[j] = (int)pairs_[j].size();
            for (const auto& pair : b_min_pairs_d[j])
                pairs_[j].push_back(pair);
        }

        /* Minimize `buffer_min_pairs_` */
        for (uint64_t ij : buffer_redundent_pairs_[d]) {
            int i, j;
            ut::GetPair(ij, i, j);
            while (j != -1) {
                MMilnor gcd = gcdLF(leads[i].m, leads[j].m);
                MMilnor m2 = divLF(leads[i].m, gcd);
                if (j < (int)buffer_min_pairs_[d].size()) {
                    auto p = std::find_if(buffer_min_pairs_[d][j].begin(), buffer_min_pairs_[d][j].end(), [&m2](const CriticalPair& c) { return c.m2 == m2; });
                    /* Remove it from `buffer_min_pairs_` */
                    if (p != buffer_min_pairs_[d][j].end()) {
                        p->i2 = -1;
                        break;
                    }
                }

                /* Reduce (i, j) */
                auto c = pairs_[j].begin();
                auto end = pairs_[j].end();
                for (; c < end; ++c) {
                    if (divisibleLF(c->m2, m2)) {
                        MMilnor m1 = divLF(leads[j].m, gcd);
                        if (gcdLF(c->m1, m1)) {
                            j = -1;
                            break;
                        }
                        else {
                            j = c->i1;
                            if (i > j)
                                std::swap(i, j);
                            break;
                        }
                    }
                }
#ifndef NDEBUG
                if (c == end)
                    throw MyException(0xfa5db14U, "Should not happen because pairs_ is groebner");
#endif
            }
        }

        if ((size_t)d >= min_pairs_.size())
            min_pairs_.resize((size_t)d + 1);
        for (int j = 0; j < (int)b_min_pairs_d.size(); ++j)
            for (int i = 0; i < (int)b_min_pairs_d[j].size(); ++i)
                if (b_min_pairs_d[j][i].i2 != -1) {
                    if (min_pairs_[d].size() <= (size_t)j)
                        min_pairs_[d].resize((size_t)j + 1);
                    min_pairs_[d][j].push_back(old_sizes[j] + i);
                }

        /* Delete `bufbuffer_redundent_pairs_[d]` */
        buffer_redundent_pairs_.erase(d);
    }

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const MMod1d& leads, MMod mon, const array& basis_degrees)
    {
        std::vector<std::pair<int, CriticalPair>> new_pairs(leads.size());
        size_t s = leads.size();
        ut::Range range(0, (int)leads.size());

        std::for_each(std::execution::seq, range.begin(), range.end(), [&](int i) {
            if (leads[i].v == mon.v) {
                int d_pair = lcmLF(leads[i].m, mon.m).deg() + basis_degrees[mon.v];
                if (d_pair <= deg_trunc_) {
                    new_pairs[i].first = d_pair;
                    CriticalPair::SetFromLM(new_pairs[i].second, leads[i].m, mon.m, (int)i, (int)s);
                }
                else
                    new_pairs[i].first = -1;
            }
            else
                new_pairs[i].first = -1;
        });

        /* Compute sigma_j */
        int d_mon = mon.m.deg() + basis_degrees[mon.v];
        singles_.resize(s + 1);
        for (int i : mon.m) {
            MMilnor m = MMilnor::FromIndex(i);
            int d = d_mon + m.deg();
            if (d <= deg_trunc_) {
                buffer_singles_[d].resize(s + 1);
                buffer_singles_[d][s].push_back((int)singles_[s].size());
                singles_[s].push_back(CriticalPair::Single(m, (int)s));
            }
        }

        /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
        for (size_t j = 1; j < new_pairs.size(); ++j) {
            if (new_pairs[j].first != -1) {
                for (size_t i = 0; i < j; ++i) {
                    if (new_pairs[i].first != -1) {
                        if (divisibleLF(new_pairs[i].second.m2, new_pairs[j].second.m2)) {
                            new_pairs[j].first = -1;
                            break;
                        }
                        else if (divisibleLF(new_pairs[j].second.m2, new_pairs[i].second.m2))
                            new_pairs[i].first = -1;
                        else if (!gcdLF(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                            int dij = lcmLF(leads[i].m, leads[j].m).deg() + basis_degrees[mon.v];
                            if (dij <= deg_trunc_)
                                buffer_redundent_pairs_[dij].insert(ut::BindPair((uint32_t)i, (uint32_t)j));
                        }
                    }
                }
            }
        }
        for (size_t i = 0; i < new_pairs.size(); ++i) {
            if (new_pairs[i].first != -1) {
                buffer_min_pairs_[new_pairs[i].first].resize(s + 1);
                buffer_min_pairs_[new_pairs[i].first][s].push_back(new_pairs[i].second);
            }
        }
    }
};

class Groebner
{
private:
    using TypeIndexKey = int;
    using TypeIndex = std::unordered_map<TypeIndexKey, array>;

private:
    GbCriPairs gb_pairs_; /* Groebner basis of critical pairs */
    array basis_degrees_;
    Mod1d data_;

    /* Caches */
    MMod1d leads_;    /* Leading monomials */
    TypeIndex index_; /* Cache for fast divisibility test */

public:
    Groebner(int deg_trunc, array basis_degrees) : gb_pairs_(deg_trunc), basis_degrees_(basis_degrees) {}

    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    Groebner(int deg_trunc, array basis_degrees, Mod1d polys) : gb_pairs_(deg_trunc), basis_degrees_(basis_degrees), data_(std::move(polys))
    {
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

    /* Initialize from `polys` and `gb_pairs` where `polys` already forms a Groebner basis. */
    Groebner(Mod1d polys, array basis_degrees, GbCriPairs gb_pairs) : gb_pairs_(std::move(gb_pairs)), basis_degrees_(std::move(basis_degrees)), data_(std::move(polys))
    {
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

private:
    static TypeIndexKey Key(MMod lead)
    {
        auto p = lead.m.begin();
        while (p != lead.m.end())
            ++p;
        return TypeIndexKey{lead.v + (*p << 20)}; /* biggest i such that (MMILNOR_ONE >> (i - 1)) & lead.m is nonzero */
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(MMod mon) const
    {
        if (!mon.m) {
            auto key = TypeIndexKey{mon.v};
            auto p = index_.find(key);
            if (p != index_.end()) {
                for (int k : p->second) {
                    if (!leads_[k].m)
                        return k;
                }
            }
        }
        for (int i : mon.m) {
            auto key = TypeIndexKey{mon.v + ((i + 1) << 20)};
            auto p = index_.find(key);
            if (p != index_.end()) {
                for (int k : p->second) {
                    if (divisibleLF(leads_[k].m, mon.m))
                        return k;
                }
            }
        }
        return -1;
    }

public: /* Getters and Setters */
    const GbCriPairs& gb_pairs() const
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
    const array& basis_degs() const
    {
        return basis_degrees_;
    }
    /* This function will erase `gb_pairs_.buffer_min_pairs[d]` */
    PCriticalPair1d pairs(int d)
    {
        return gb_pairs_.pairs_for_gb(d);
    }
    auto size() const
    {
        return data_.size();
    }
    auto& operator[](size_t index) const
    {
        return data_[index];
    }
    void push_back(Mod g)
    {
        MMod mv = g.GetLead();
        gb_pairs_.AddToBuffers(leads_, mv, basis_degrees_);

        leads_.push_back(mv);
        index_[Key(mv)].push_back((int)data_.size());
        data_.push_back(std::move(g));
    }
    void AddPairsAndMinimize(int deg)
    {
        gb_pairs_.AddAndMinimize(leads_, deg);
    }
    Mod1d MinSyzOfGb() const
    {
        return gb_pairs_.MinSyzOfGb();
    }

    const Mod1d& data() const
    {
        return data_;
    }

public:
    Mod Reduce(CriticalPair& p) const
    {
        Mod result;
        if (p.i1 == -1) {
            result = Milnor(p.m2) * data_[p.i2];
            p.x = Mod(MMod{p.m2, p.i2});
        }
        else {
            result = Milnor(p.m1) * data_[p.i1] + Milnor(p.m2) * data_[p.i2];
            p.x = Mod(MMod{p.m1, p.i1}) + Mod(MMod{p.m2, p.i2});
        }

        size_t index = 0;
        while (index < result.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.data[index]);
            if (gb_index != -1) {
                MMilnor m;
                detail::Reduce(result, data_[gb_index], index, m);
                p.x += Mod(MMod{m, gb_index});
            }
            else
                ++index;
        }
        return result;
    }

    Mod Reduce(Mod x) const
    {
        size_t index = 0;
        while (index < x.data.size()) {
            int gb_index = IndexOfDivisibleLeading(x.data[index]);
            if (gb_index != -1)
                detail::Reduce(x, data_[gb_index], index);
            else
                ++index;
        }
        return x;
    }
};

/**
 * Comsume relations from 'rels` and `gb.gb_pairs_` in degree `<= deg`
 * `min_gb` stores the minimal generating set of gb.
 */
void AddRels(Groebner& gb, const Mod1d& rels, int deg, array& min_gb, Mod1d& redundent_gb);

/*************************************************************
 ********************** Min Res ******************************
 **************************************************************/

inline constexpr int S_MAX = 1024;

/* MMod.v = bind(S_MAX - s, index) */
inline int GetS(int v)
{
    return S_MAX - (v >> 16);
}
inline int GetV(int s, int index)
{
    return index + ((S_MAX - s) << 16);
}
inline int GetBaseDeg(const array2d& basis_degrees, int v)
{
    return basis_degrees[GetS(v)][v & 0xffff];
}

struct CriticalPairMinRes
{
    int i1 = -1, i2 = -1;
    int s = -1;
    MMilnor m1, m2;

    /* Compute the pair for two leading monomials. */
    CriticalPairMinRes() = default;
    static void SetFromLM(CriticalPairMinRes& result, MMilnor lead1, MMilnor lead2, int i, int j, int s)
    {
        MMilnor gcd = gcdLF(lead1, lead2);
        result.m1 = divLF(lead2, gcd);
        result.m2 = divLF(lead1, gcd);
        result.i1 = i;
        result.i2 = j;
        result.s = s;
    }
    static CriticalPairMinRes Single(MMilnor m2, int j, int s)
    {
        CriticalPairMinRes result;
        result.m2 = m2;
        result.i1 = -1;
        result.i2 = j;
        result.s = s;
        return result;
    }
};
using CriticalPairMinRes1d = std::vector<CriticalPairMinRes>;
using CriticalPairMinRes2d = std::vector<CriticalPairMinRes1d>;
using PCriticalPairMinRes1d = std::vector<CriticalPairMinRes*>;
using PCriticalPairMinRes2d = std::vector<PCriticalPairMinRes1d>;

inline void print(const Mod& x)  ////
{
    if (x.data.empty()) {
        std::cout << '0';
    }
    for (auto pm = x.data.begin(); pm != x.data.end(); ++pm) {
        if (pm != x.data.begin())
            std::cout << '+';
        for (int i : pm->m)
            std::cout << "P_{" << MMILNOR_GEN_I[i] << MMILNOR_GEN_J[i] << '}';
        std::cout << "v_{" << GetS(pm->v) << ',' << (pm->v & 0xffff) << '}';
    }
}

/* Groebner basis of critical pairs */
class GbCriPairsMinRes
{
private:
    int deg_trunc_;                                                      /* Truncation degree */
    CriticalPairMinRes2d pairs_;                                         /* `pairs_[j]` is the set of pairs (i, j) with given j */
    CriticalPairMinRes2d singles_;                                       /* `singles[d]` is the set of gen *  */
    array3d min_pairs_;                                                  /* Minimal generating set of `pairs_`. `min_pairs_` stores indices of `pairs` */
    std::map<int, CriticalPairMinRes2d> buffer_min_pairs_;               /* To generate `min_pairs_` and for computing Sij */
    std::map<int, std::unordered_set<uint64_t>> buffer_redundent_pairs_; /* Used to minimize `buffer_min_pairs_` */
    std::map<int, array2d> buffer_singles_;                              /* For computing Sj. `buffer_singles_` stores indices of singles_ */

public:
    GbCriPairsMinRes(int d_trunc) : deg_trunc_(d_trunc) {}
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
    bool empty_pairs_for_gb_d(int d) const
    {
        return buffer_min_pairs_.find(d) == buffer_min_pairs_.end() && buffer_singles_.find(d) == buffer_singles_.end();
    }
    bool empty_pairs_for_gb() const
    {
        return buffer_min_pairs_.empty() && buffer_singles_.empty();
    }
    PCriticalPairMinRes2d pairs_for_gb(int d)  // Warning: return pointers
    {
        PCriticalPairMinRes2d result;
        if ((size_t)d < min_pairs_.size())
            for (size_t j = 0; j < min_pairs_[d].size(); ++j)
                for (auto i : min_pairs_[d][j]) {
                    if (result.size() <= pairs_[j][i].s)
                        result.resize(size_t(pairs_[j][i].s + 1));
                    result[pairs_[j][i].s].push_back(&pairs_[j][i]);
                }
        if (buffer_singles_.find(d) != buffer_singles_.end())
            for (size_t j = 0; j < buffer_singles_[d].size(); ++j)
                for (auto i : buffer_singles_[d][j]) {
                    if (result.size() <= singles_[j][i].s)
                        result.resize(size_t(singles_[j][i].s + 1));
                    result[singles_[j][i].s].push_back(&singles_[j][i]);
                }

        buffer_min_pairs_.erase(d);
        buffer_singles_.erase(d);
        return result;
    }

    /* Minimize `buffer_min_pairs_[d]` and maintain `pairs_` */
    void AddAndMinimize(const MMod1d& leads, int d)
    {
        /* Add to the Groebner basis of critical pairs */
        if (buffer_min_pairs_.find(d) == buffer_min_pairs_.end())
            return;
        auto& b_min_pairs_d = buffer_min_pairs_.at(d);
        if (pairs_.size() < b_min_pairs_d.size())
            pairs_.resize(b_min_pairs_d.size());
        array old_sizes(b_min_pairs_d.size());
        for (int j = 0; j < (int)b_min_pairs_d.size(); ++j) {
            old_sizes[j] = (int)pairs_[j].size();
            for (const auto& pair : b_min_pairs_d[j])
                pairs_[j].push_back(pair);
        }

        /* Minimize `buffer_min_pairs_` */
        for (uint64_t ij : buffer_redundent_pairs_[d]) {
            int i, j;
            ut::GetPair(ij, i, j);
            while (j != -1) {
                MMilnor gcd = gcdLF(leads[i].m, leads[j].m);
                MMilnor m2 = divLF(leads[i].m, gcd);
                if (j < (int)buffer_min_pairs_[d].size()) {
                    auto p = std::find_if(buffer_min_pairs_[d][j].begin(), buffer_min_pairs_[d][j].end(), [&m2](const CriticalPairMinRes& c) { return c.m2 == m2; });
                    /* Remove it from `buffer_min_pairs_` */
                    if (p != buffer_min_pairs_[d][j].end()) {
                        p->i2 = -1;
                        break;
                    }
                }

                /* Reduce (i, j) */
                auto c = pairs_[j].begin();
                auto end = pairs_[j].end();
                for (; c < end; ++c) {
                    if (divisibleLF(c->m2, m2)) {
                        MMilnor m1 = divLF(leads[j].m, gcd);
                        if (gcdLF(c->m1, m1)) {
                            j = -1;
                            break;
                        }
                        else {
                            j = c->i1;
                            if (i > j)
                                std::swap(i, j);
                            break;
                        }
                    }
                }
#ifndef NDEBUG
                if (c == end) {
                    std::cout << "i=" << i << '\n';
                    std::cout << "j=" << j << '\n';
                    std::cout << "leads[i]=";
                    print(Mod(leads[i]));
                    std::cout << "\nleads[j]=";
                    print(Mod(leads[j]));
                    std::cout << '\n';
                    std::cout << "m2=" << Milnor(m2) << '\n';
                    for (size_t k = 0; k < pairs_[j].size(); ++k)
                        std::cout << "pairs_[j][k].m2" << Milnor(pairs_[j][k].m2) << '\n';
                    throw MyException(0xfa5db14U, "Should not happen because pairs_ is groebner");
                }
#endif
            }
        }

        if ((size_t)d >= min_pairs_.size())
            min_pairs_.resize((size_t)d + 1);
        for (int j = 0; j < (int)b_min_pairs_d.size(); ++j)
            for (int i = 0; i < (int)b_min_pairs_d[j].size(); ++i)
                if (b_min_pairs_d[j][i].i2 != -1) {
                    if (min_pairs_[d].size() <= (size_t)j)
                        min_pairs_[d].resize((size_t)j + 1);
                    min_pairs_[d][j].push_back(old_sizes[j] + i);
                }

        /* Delete `bufbuffer_redundent_pairs_[d]` */
        buffer_redundent_pairs_.erase(d);
    }

    /* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
    ** `buffer_min_pairs_` will become a Groebner basis at this stage.
    */
    void AddToBuffers(const MMod1d& leads, MMod mon, int d_mon_base)
    {
        std::vector<std::pair<int, CriticalPairMinRes>> new_pairs(leads.size());
        size_t s = leads.size();
        ut::Range range(0, (int)leads.size());

        std::for_each(std::execution::seq, range.begin(), range.end(), [&](int i) {
            if (leads[i].v == mon.v) {
                int d_pair = lcmLF(leads[i].m, mon.m).deg() + d_mon_base;
                if (d_pair <= deg_trunc_) {
                    new_pairs[i].first = d_pair;
                    CriticalPairMinRes::SetFromLM(new_pairs[i].second, leads[i].m, mon.m, (int)i, (int)s, GetS(mon.v));
                }
                else
                    new_pairs[i].first = -1;
            }
            else
                new_pairs[i].first = -1;
        });

        /* Compute sigma_j */
        int d_mon = mon.m.deg() + d_mon_base;
        singles_.resize(s + 1);
        for (int i : mon.m) {
            MMilnor m = MMilnor::FromIndex(i);
            int d = d_mon + m.deg();
            if (d <= deg_trunc_) {
                buffer_singles_[d].resize(s + 1);
                buffer_singles_[d][s].push_back((int)singles_[s].size());
                singles_[s].push_back(CriticalPairMinRes::Single(m, (int)s, GetS(mon.v)));
            }
        }

        /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
        for (size_t j = 1; j < new_pairs.size(); ++j) {
            if (new_pairs[j].first != -1) {
                for (size_t i = 0; i < j; ++i) {
                    if (new_pairs[i].first != -1) {
                        if (divisibleLF(new_pairs[i].second.m2, new_pairs[j].second.m2)) {
                            new_pairs[j].first = -1;
                            break;
                        }
                        else if (divisibleLF(new_pairs[j].second.m2, new_pairs[i].second.m2))
                            new_pairs[i].first = -1;
                        else if (!gcdLF(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                            int dij = lcmLF(leads[i].m, leads[j].m).deg() + d_mon_base;
                            if (dij <= deg_trunc_) {
                                /*if (i == 31 && j == 80) {
                                    std::cout << "new_pairs[i].second.m2=" << Milnor(new_pairs[i].second.m2) << '\n';
                                    std::cout << "new_pairs[j].second.m2=" << Milnor(new_pairs[j].second.m2) << '\n';
                                    std::cout << "s=" << s << "  \n";
                                    std::cout << "mon=";
                                    print(Mod(mon));
                                    std::cout << '\n';
                                    for (size_t k = 0; k < leads.size(); ++k) {
                                        std::cout << "leads[" << k << "]=";
                                        print(Mod(leads[k]));
                                        std::cout << '\n';
                                    }
                                    throw MyException(0, "BUG");
                                }*/
                                buffer_redundent_pairs_[dij].insert(ut::BindPair((uint32_t)i, (uint32_t)j));
                            }
                        }
                    }
                }
            }
        }
        for (size_t i = 0; i < new_pairs.size(); ++i) {
            if (new_pairs[i].first != -1) {
                buffer_min_pairs_[new_pairs[i].first].resize(s + 1);
                buffer_min_pairs_[new_pairs[i].first][s].push_back(new_pairs[i].second);
            }
        }
    }
};

class GroebnerMinRes
{
private:
    using TypeIndexKey = uint64_t;
    using TypeIndex = std::unordered_map<TypeIndexKey, array>;

private:
    GbCriPairsMinRes gb_pairs_; /* Groebner basis of critical pairs */
    array2d basis_degrees_;
    Mod1d data_;
    Mod2d kernel_; /* kernel[s] is the generators of Kernel of $F_{s}\to F_{s-1}$ */

    /* Caches */
    MMod1d leads_;    /* Leading monomials */
    TypeIndex index_; /* Cache for fast divisibility test */

public:
    GroebnerMinRes(int deg_trunc, array2d basis_degrees) : gb_pairs_(deg_trunc), basis_degrees_(basis_degrees) {}

    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    GroebnerMinRes(int deg_trunc, array2d basis_degrees, Mod1d polys) : gb_pairs_(deg_trunc), basis_degrees_(basis_degrees), data_(std::move(polys))
    {
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

    /* Initialize from `polys` and `gb_pairs` where `polys` already forms a Groebner basis. */
    GroebnerMinRes(Mod1d polys, array2d basis_degrees, GbCriPairsMinRes gb_pairs) : gb_pairs_(std::move(gb_pairs)), basis_degrees_(std::move(basis_degrees)), data_(std::move(polys))
    {
        for (int i = 0; i < (int)data_.size(); ++i) {
            leads_.push_back(data_[i].GetLead());
            index_[Key(data_[i].GetLead())].push_back(i);
        }
    }

private:
    static TypeIndexKey Key(MMod lead)
    {
        auto p = lead.m.begin();
        while (p != lead.m.end())
            ++p;
        return TypeIndexKey(lead.v) + (TypeIndexKey(*p) << 32); /* biggest i such that (MMILNOR_ONE >> (i - 1)) & lead.m is nonzero */
    }

    /* Return -1 if not found */
    int IndexOfDivisibleLeading(MMod mon) const
    {
        if (!mon.m) {
            auto key = TypeIndexKey(mon.v);
            auto p = index_.find(key);
            if (p != index_.end()) {
                for (int k : p->second) {
                    if (!leads_[k].m)
                        return k;
                }
            }
        }
        for (int i : mon.m) {
            auto key = TypeIndexKey(mon.v) + (TypeIndexKey(i + 1) << 32);
            auto p = index_.find(key);
            if (p != index_.end()) {
                for (int k : p->second) {
                    if (divisibleLF(leads_[k].m, mon.m))
                        return k;
                }
            }
        }
        return -1;
    }

public: /* Getters and Setters */
    const GbCriPairsMinRes& gb_pairs() const
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
    /* This function will erase `gb_pairs_.buffer_min_pairs[d]` */
    PCriticalPairMinRes2d pairs(int d)
    {
        return gb_pairs_.pairs_for_gb(d);
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
    void push_back(Mod g)
    {
        MMod mv = g.GetLead();
        gb_pairs_.AddToBuffers(leads_, mv, GetBaseDeg(basis_degrees_, mv.v));  // TODO: modify the counterpart in algebras/groebner.h

        leads_.push_back(mv);
        index_[Key(mv)].push_back((int)data_.size());
        data_.push_back(std::move(g));
    }
    void push_back_kernel(Mod g, int d, int s)
    {
        if (kernel_.size() <= s) {
            kernel_.resize(size_t(s + 1));
            basis_degrees_.resize(size_t(s + 2));
        }
        kernel_[s].push_back(g);
        basis_degrees_[size_t(s + 1)].push_back(d);

        g.data.push_back(MMod{MMilnor(0), GetV(s + 1, (int)basis_degrees_[size_t(s + 1)].size() - 1)});
        push_back(std::move(g));
    }
    void AddPairsAndMinimize(int deg)
    {
        gb_pairs_.AddAndMinimize(leads_, deg);
    }

    const Mod1d& data() const
    {
        return data_;
    }
    const Mod2d& kernel() const
    {
        return kernel_;
    }

public:
    Mod Reduce(CriticalPairMinRes& p) const
    {
        Mod result;
        if (p.i1 == -1) {
            result = Milnor(p.m2) * data_[p.i2];
        }
        else {
            result = Milnor(p.m1) * data_[p.i1] + Milnor(p.m2) * data_[p.i2];
        }

        size_t index = 0;
        while (index < result.data.size()) {
            int gb_index = IndexOfDivisibleLeading(result.data[index]);
            if (gb_index != -1) {
                MMilnor m;  ////
                detail::Reduce(result, data_[gb_index], index, m);
            }
            else
                ++index;
        }
        return result;
    }

    Mod Reduce(Mod x) const
    {
        size_t index = 0;
        while (index < x.data.size()) {
            int gb_index = IndexOfDivisibleLeading(x.data[index]);
            if (gb_index != -1)
                detail::Reduce(x, data_[gb_index], index);
            else
                ++index;
        }
        return x;
    }

    Mod HeadReduce(CriticalPairMinRes& p) const
    {
        Mod result;
        if (p.i1 == -1) {
            result = Milnor(p.m2) * data_[p.i2];
        }
        else {
            result = Milnor(p.m1) * data_[p.i1] + Milnor(p.m2) * data_[p.i2];
        }

        size_t index = 0;
        while (index < result.data.size() && GetS(result.data[index].v) == p.s) {
            int gb_index = IndexOfDivisibleLeading(result.data[index]);
            if (gb_index != -1)
                detail::Reduce(result, data_[gb_index], index);
            else
                ++index;
        }
        return result;
    }

    Mod HeadReduce(Mod x, int s) const
    {
        size_t index = 0;
        while (index < x.data.size() && GetS(x.data[index].v) == s) {
            int gb_index = IndexOfDivisibleLeading(x.data[index]);
            if (gb_index != -1)
                detail::Reduce(x, data_[gb_index], index);
            else
                ++index;
        }

        return x;
    }
};

/**
 * Comsume relations from 'rels` and `gb.gb_pairs_` in degree `<= deg`
 * `min_gb` stores the minimal generating set of gb.
 */
void AddRelsMinRes(GroebnerMinRes& gb, const Mod1d& rels, int deg);

}  // namespace steenrod

#endif