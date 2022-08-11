#include "groebner_steenrod.h"
#include "algebras/benchmark.h"
#include "algebras/database.h"
#include <atomic>
#include <cstring>
#include <mutex>

namespace steenrod {

/********************************************************
 *                    class CriMilnors
 ********************************************************/

CriMilnor1d CriMilnors::Criticals(int t)
{
    CriMilnor1d result;
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

void CriMilnors::Minimize(const MMod1d& leads, int t)
{
    if (buffer_redundent_pairs_.empty() || buffer_redundent_pairs_.begin()->first != t)
        return;
    if (buffer_min_pairs_.empty() || buffer_min_pairs_.begin()->first != t)
        return;

    constexpr uint64_t NULL_J = ~uint64_t(0);
    auto& b_min_pairs_t = buffer_min_pairs_.begin()->second;
    for (uint64_t ij : buffer_redundent_pairs_.begin()->second) {
        uint64_t i, j;
        ut::UnBind(ij, i, j);
        while (j != NULL_J) {
            MMilnor gcd = gcdLF(leads[i].m_no_weight(), leads[j].m_no_weight());
            MMilnor m2 = divLF(leads[i].m(), gcd);
            if (j < (int)b_min_pairs_t.size()) {
                auto p = std::find_if(b_min_pairs_t[j].begin(), b_min_pairs_t[j].end(), [&m2](const CriMilnor& c) { return c.m2 == m2; });
                /* Mark it to be removed from `buffer_min_pairs_` */
                if (p != b_min_pairs_t[j].end()) {
                    p->i2 = -1;
                    break;
                }
            }

            /* Reduce (i, j) */
            auto c = gb_[j].begin();
            auto end = gb_[j].end();
            for (; c < end; ++c) {
                if (divisibleLF(c->m2, m2)) {
                    MMilnor m1 = divLF(leads[j].m(), gcd);
                    if (gcdLF(c->m1, m1)) {
                        j = NULL_J;
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
                throw MyException(0xfa5db14U, "Should not happen because gb_ is groebner");
#endif
        }
    }

    /* Delete `buffer_redundent_pairs_[t]` */
    buffer_redundent_pairs_.erase(buffer_redundent_pairs_.begin());
}

/**
 * Populate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
 * `buffer_min_pairs_` will be reduced later by `CPMilnors::Minimize()`.
 */
void CriMilnors::AddToBuffers(const MMod1d& leads, MMod mon, int t_v)
{
    size_t lead_size = leads.size();
    if (gb_.size() < lead_size + 1)
        gb_.resize(lead_size + 1);
    std::vector<std::pair<int, CriMilnor>> new_pairs(lead_size);

    /* Populate `new_pairs` */
    for (size_t i = 0; i < lead_size; ++i) {
        new_pairs[i].first = -1;
        if (leads[i].v_raw() == mon.v_raw()) {
            int d_pair = lcmLF(leads[i].m_no_weight(), mon.m_no_weight()).deg() + t_v;
            if (d_pair <= t_trunc_) {
                new_pairs[i].first = d_pair;
                CriMilnor::SetFromLM(new_pairs[i].second, leads[i].m(), mon.m(), (int)i, (int)lead_size);
            }
        }
    }

    /* Remove some new pairs to form Groebner basis and discover redundent pairs */
    for (size_t j = 1; j < new_pairs.size(); ++j) {
        if (new_pairs[j].first != -1) {
            for (size_t i = 0; i < j; ++i) {
                if (new_pairs[i].first != -1) {
                    if (divisibleLF(new_pairs[i].second.m2, new_pairs[j].second.m2)) {
                        new_pairs[j].first = -1;
                        break;
                    }
                    else if (divisibleLF(new_pairs[j].second.m2, new_pairs[i].second.m2)) {
                        new_pairs[i].first = -1;
                    }
                    else if (!gcdLF(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                        int dij = lcmLF(leads[i].m_no_weight(), leads[j].m_no_weight()).deg() + t_v;
                        if (dij <= t_trunc_)
                            buffer_redundent_pairs_[dij].insert(ut::Bind(i, j));
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < new_pairs.size(); ++i) {
        if (new_pairs[i].first != -1) {
            gb_[lead_size].push_back(new_pairs[i].second);
            buffer_min_pairs_[new_pairs[i].first].resize(lead_size + 1);
            buffer_min_pairs_[new_pairs[i].first][lead_size].push_back(new_pairs[i].second);
        }
    }

    /* Populate `buffer_singles_` */
    int t_mon = mon.deg_m() + t_v;
    for (int i : mon.m_no_weight()) {
        MMilnor m = MMilnor::FromIndex(i);
        int t = t_mon + m.deg();
        if (t <= t_trunc_)
            buffer_singles_[t].push_back(CriMilnor::Single(m, (int)lead_size));
    }
}

void CriMilnors::init(const MMod1d& leads, const int1d& basis_degrees, int t_min_buffer)
{
    MMod1d tmp_leads;
    for (auto& mon : leads) {
        int t_v = basis_degrees[mon.v()];
        AddToBuffers(tmp_leads, mon, t_v);
        tmp_leads.push_back(mon);
    }

    for (auto it = buffer_min_pairs_.begin(); it != buffer_min_pairs_.end();) {
        if (it->first < t_min_buffer)
            it = buffer_min_pairs_.erase(it);
        else
            ++it;
    }
    for (auto it = buffer_redundent_pairs_.begin(); it != buffer_redundent_pairs_.end();) {
        if (it->first < t_min_buffer)
            it = buffer_redundent_pairs_.erase(it);
        else
            ++it;
    }
    for (auto it = buffer_singles_.begin(); it != buffer_singles_.end();) {
        if (it->first < t_min_buffer)
            it = buffer_singles_.erase(it);
        else
            ++it;
    }
}

/********************************************************
 *                    class GroebnerX2m
 ********************************************************/

/** Return i such that mon divides leads[i].
 *
 * Return -1 if not found.
 * @param indices indices[v_raw] stores the indices of elements of leads with the given v_raw.
 */
int IndexOfDivisibleLeading(const MMod1d& leads, const std::unordered_map<uint64_t, int1d>& indices, MMod mon)
{
    auto key = mon.v_raw();
    auto p = indices.find(key);
    if (p != indices.end())
        for (int k : p->second)
            if (divisibleLF(leads[k], mon))
                return k;
    return -1;
}

GroebnerX2m::GroebnerX2m(int t_trunc, int stem_trunc, Mod2d data, int2d basis_degrees, std::map<int, int>& latest_st) : t_trunc_(t_trunc), stem_trunc_(stem_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees))
{
    leads_.resize(gb_.size());
    indices_.resize(gb_.size());

    for (size_t s = 0; s < gb_.size(); ++s) {
        for (int j = 0; j < (int)gb_[s].size(); ++j) {
            leads_[s].push_back(gb_[s][j].GetLead());
            indices_[s][gb_[s][j].GetLead().v_raw()].push_back(j);
        }
    }

    for (size_t s = 0; s < gb_.size(); ++s) {
        criticals_.push_back(CriMilnors(std::min(t_trunc_, int(s) + stem_trunc_ + 2)));
        criticals_.back().init(leads_[s], basis_degrees_[s], latest_st[(int)s] + 1);
    }
}

Mod GroebnerX2m::Reduce(Mod x2m, size_t s) const
{
    Mod tmp;
    tmp.data.reserve(50);

    size_t index;
    index = 0;
    while (index < x2m.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], x2m.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(x2m.data[index], gb_[s][gb_index].data[0]);
            x2m.iaddmulMay(m, gb_[s][gb_index], tmp);
        }
        else
            ++index;
    }

    return x2m;
}

Mod GroebnerX2m::Reduce(const CriMilnor& p, size_t s) const
{
    Mod result;

    Mod tmp;
    tmp.data.reserve(50);

    if (p.i1 >= 0)
        result.iaddmulMay(p.m1, gb_[s][p.i1], tmp).iaddmulMay(p.m2, gb_[s][p.i2], tmp);
    else
        result.iaddmulMay(p.m2, gb_[s][p.i2], tmp);

    size_t index;
    index = 0;
    while (index < result.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], result.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.data[index], gb_[s][gb_index].data[0]);
            result.iaddmulMay(m, gb_[s][gb_index], tmp);
        }
        else
            ++index;
    }

    return result;
}

Mod1d GroebnerX2m::AddRels(size_t s, int t)
{
    Mod tmp_Mod;

    /* Populate `rels_tmp` */
    resize(s + 1);
    Mod1d rels_tmp;

    criticals_[s].Minimize(leads_[s], t);
    CriMilnor1d pairs_st = criticals_[s].Criticals(t);
    if (!pairs_st.empty()) {
        rels_tmp.resize(pairs_st.size());
        ut::for_each_seq((int)rels_tmp.size(), [&](size_t i) { rels_tmp[i] = Reduce(pairs_st[i], s); });
    }

    Mod1d rels_st;
    if (rels_tmp.empty())
        return rels_st;

    /* Triangulate these relations */
    for (size_t i = 0; i < rels_tmp.size(); ++i) {
        steenrod::Reduce(rels_tmp[i], rels_st, tmp_Mod);
        if (rels_tmp[i])
            rels_st.push_back(std::move(rels_tmp[i]));
    }

    /* Add these relations */
    for (size_t i = 0; i < rels_st.size(); ++i)
        push_back(rels_st[i], s);
    return rels_st;
}

/********************************************************
 *                    class AdamsRes
 ********************************************************/

AdamsRes::AdamsRes(int t_trunc, int stem_trunc, DataMRes2d data, int2d basis_degrees, Mod2d data_x2m, int2d basis_degrees_x2m, std::map<int, int>& latest_st)
    : t_trunc_(t_trunc), stem_trunc_(stem_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees)), gb_x2m_(t_trunc, stem_trunc, std::move(data_x2m), std::move(basis_degrees_x2m), latest_st)
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
            indices_[s][gb_[s][j].x1.GetLead().v_raw()].push_back(j);
        }
    }

    for (size_t s = 0; s < gb_.size(); ++s) {
        criticals_.push_back(CriMilnors(std::min(t_trunc_, int(s) + stem_trunc_ + 2)));
        criticals_.back().init(leads_[s], basis_degrees_[s], latest_st[(int)s] + 1);
    }
}

CriMilnor1d AdamsRes::Criticals(size_t s, int t, Mod1d& rels_x2m)
{
    rels_x2m = gb_x2m_.AddRels(s, t);
    criticals_[s].Minimize(leads_[s], t);
    CriMilnor1d cris = criticals_[s].Criticals(t);
    std::vector<Filtr> fils(cris.size());
    for (size_t i = 0; i < cris.size(); ++i)
        fils[i] = gb_[s][cris[i].i2].fil + cris[i].m2.w_may();
    auto indices = ut::size_t_range(cris.size());
    std::stable_sort(indices.begin(), indices.end(), [&fils](size_t i, size_t j) { return fils[j] < fils[i]; });
    CriMilnor1d result;
    result.reserve(cris.size());
    for (size_t i = 0; i < cris.size(); ++i)
        result.push_back(cris[indices[i]]);

    Mod1d x2ms;
    Mod tmp_Mod;
    for (auto& cp : result) {
        Mod x2m = ReduceX2m(cp, s);
        steenrod::Reduce(x2m, x2ms, tmp_Mod);
        if (x2m)
            x2ms.push_back(std::move(x2m));
        else
            cp.i2 = -1;
    }
    ut::RemoveIf(result, [](const CriMilnor& cp) { return cp.i2 == -1; });
    return result;
}

DataMRes AdamsRes::Reduce(const CriMilnor& cp, size_t s) const
{
    DataMRes result;

    Milnor tmp_a;
    Mod tmp_x1, tmp_x2;
    tmp_a.data.reserve(32);
    tmp_x1.data.reserve(64);
    tmp_x2.data.reserve(64);

    if (cp.i1 >= 0) {
        result.x1.iaddmulP(cp.m1, gb_[s][cp.i1].x1, tmp_a, tmp_x1, tmp_x2).iaddmulP(cp.m2, gb_[s][cp.i2].x1, tmp_a, tmp_x1, tmp_x2);
        result.x2.iaddmulP(cp.m1, gb_[s][cp.i1].x2, tmp_a, tmp_x1, tmp_x2).iaddmulP(cp.m2, gb_[s][cp.i2].x2, tmp_a, tmp_x1, tmp_x2);
        result.x2m.iaddmulMay(cp.m1, gb_[s][cp.i1].x2m, tmp_x1).iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_x1);
    }
    else {
        result.x1.iaddmulP(cp.m2, gb_[s][cp.i2].x1, tmp_a, tmp_x1, tmp_x2);
        result.x2.iaddmulP(cp.m2, gb_[s][cp.i2].x2, tmp_a, tmp_x1, tmp_x2);
        result.x2m.iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_x1);
    }
    result.fil = gb_[s][cp.i2].fil + cp.m2.w_may();

    size_t index;
    index = 0;
    while (index < result.x1.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], result.x1.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.x1.data[index], gb_[s][gb_index].x1.data[0]);
            if (result.valid_x2m() && result.fil == Filtr(result.x1.data[index]))
                result.x2m.iaddmulMay(m, gb_[s][gb_index].x2m, tmp_x1);
            result.x1.iaddmulP(m, gb_[s][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
            result.x2.iaddmulP(m, gb_[s][gb_index].x2, tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }
    size_t sp1 = s + 1;
    index = 0;
    while (index < result.x2.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[sp1], indices_[sp1], result.x2.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.x2.data[index], gb_[sp1][gb_index].x1.data[0]);
            result.x2.iaddmulP(m, gb_[sp1][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }
    result.x2m = gb_x2m_.Reduce(std::move(result.x2m), s);

    return result;
}

struct IndexMMod
{
    MMod m;
    unsigned i, index;
    bool operator<(IndexMMod rhs)
    {
        return rhs.m < m;
    }
};
using IndexMMod1d = std::vector<IndexMMod>;

/**
 * results are expected to be zero initially
 */
void AdamsRes::ReduceBatch(const CriMilnor1d& cps, DataMRes1d& results, size_t s) const
{
    Mod tmp_x, tmp_x1, tmp_x2, tmp_x3;
    Milnor tmp_a;
    tmp_x.data.reserve(64);
    tmp_x1.data.reserve(64);
    tmp_x2.data.reserve(64);
    tmp_x3.data.reserve(32);
    tmp_a.data.reserve(64);

    for (size_t i = 0; i < cps.size(); ++i) {
        if (cps[i].i1 >= 0) {
            results[i].x1.iaddmulP(cps[i].m1, gb_[s][cps[i].i1].x1, tmp_a, tmp_x1, tmp_x2).iaddmulP(cps[i].m2, gb_[s][cps[i].i2].x1, tmp_a, tmp_x1, tmp_x2);
            results[i].x2.iaddmulP(cps[i].m1, gb_[s][cps[i].i1].x2, tmp_a, tmp_x1, tmp_x2).iaddmulP(cps[i].m2, gb_[s][cps[i].i2].x2, tmp_a, tmp_x1, tmp_x2);
            results[i].x2m.iaddmulMay(cps[i].m1, gb_[s][cps[i].i1].x2m, tmp_x1).iaddmulMay(cps[i].m2, gb_[s][cps[i].i2].x2m, tmp_x1);
        }
        else {
            results[i].x1.iaddmulP(cps[i].m2, gb_[s][cps[i].i2].x1, tmp_a, tmp_x1, tmp_x2);
            results[i].x2.iaddmulP(cps[i].m2, gb_[s][cps[i].i2].x2, tmp_a, tmp_x1, tmp_x2);
            results[i].x2m.iaddmulMay(cps[i].m2, gb_[s][cps[i].i2].x2m, tmp_x1);
        }
        results[i].fil = gb_[s][cps[i].i2].fil + cps[i].m2.w_may();
    }

    IndexMMod1d heap;
    for (size_t i = 0; i < results.size(); ++i)
        if (results[i].x1)
            heap.push_back(IndexMMod{results[i].x1.data[0], (unsigned)i, 0});
    std::make_heap(heap.begin(), heap.end());

    while (!heap.empty()) {
        MMod term = heap.front().m;
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], term);
        if (gb_index != -1) {
            MMilnor m = divLF(term, gb_[s][gb_index].x1.data[0]);
            mulP(m, gb_[s][gb_index].x1, tmp_x1, tmp_a);
            mulP(m, gb_[s][gb_index].x2, tmp_x2, tmp_a);
            MulMayP(m, gb_[s][gb_index].x2m, tmp_x3, tmp_a);

            while (!heap.empty() && heap.front().m == term) {

                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());

                if (results[i].valid_x2m() && results[i].fil == Filtr(results[i].x1.data[index]))
                    results[i].x2m.iaddP(tmp_x3, tmp_x);
                results[i].x1.iaddP(tmp_x1, tmp_x);
                results[i].x2.iaddP(tmp_x2, tmp_x);

                if (index < results[i].x1.data.size()) {
                    heap.back() = IndexMMod{results[i].x1.data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
        else {
            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());
                if (++index < results[i].x1.data.size()) {
                    heap.back() = IndexMMod{results[i].x1.data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
    }

    size_t sp1 = s + 1;
    heap.clear();
    for (size_t i = 0; i < results.size(); ++i)
        if (results[i].x2)
            heap.push_back(IndexMMod{results[i].x2.GetLead(), (unsigned)i, 0});
    std::make_heap(heap.begin(), heap.end());

    while (!heap.empty()) {
        MMod term = heap.front().m;
        int gb_index = IndexOfDivisibleLeading(leads_[sp1], indices_[sp1], term);
        if (gb_index != -1) {
            MMilnor m = divLF(term, gb_[sp1][gb_index].x1.data[0]);
            mulP(m, gb_[sp1][gb_index].x1, tmp_x1, tmp_a);

            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());

                results[i].x2.iaddP(tmp_x1, tmp_x);

                if (index < results[i].x2.data.size()) {
                    heap.back() = IndexMMod{results[i].x2.data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
        else {
            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());
                if (++index < results[i].x2.data.size()) {
                    heap.back() = IndexMMod{results[i].x2.data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
    }

    for (size_t i = 0; i < results.size(); ++i)
        results[i].x2m = gb_x2m_.Reduce(std::move(results[i].x2m), s);
}

Mod AdamsRes::Reduce(Mod x, size_t s) const
{
    size_t index;
    index = 0;
    Milnor tmp_a;
    Mod tmp_x1, tmp_x2;
    while (index < x.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], x.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(x.data[index], gb_[s][gb_index].x1.data[0]);
            x.iaddmulP(m, gb_[s][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }

    return x;
}

Mod AdamsRes::ReduceX2m(const CriMilnor& cp, size_t s) const
{
    Mod result;

    Mod tmp_m2;
    tmp_m2.data.reserve(100);

    if (cp.i1 >= 0)
        result.iaddmulMay(cp.m1, gb_[s][cp.i1].x2m, tmp_m2).iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_m2);
    else
        result.iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_m2);

    return gb_x2m_.Reduce(result, s);
}

class DbSteenrod : public myio::Database
{
    using Statement = myio::Statement;

private:
    std::future<void> f_;

public:
    DbSteenrod() = default;
    explicit DbSteenrod(const std::string& filename) : Database(filename) {}

public:
    void create_generators(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, diff BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_generators_x2m(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_X2m_generators (id INTEGER PRIMARY KEY, s SMALLINT, t SMALLINT);");
    }
    void create_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (id INTEGER PRIMARY KEY, x1 BLOB, x2 BLOB, x2m BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_relations_x2m(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_X2m_relations (id INTEGER PRIMARY KEY, x2m BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_time(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_time (id INTEGER PRIMARY KEY, s SMALLINT, t SMALLINT, time REAL, UNIQUE(s, t));");
    }

    void create_generators_and_delete(const std::string& table_prefix) const
    {
        create_generators(table_prefix);
        delete_from(table_prefix + "_generators");
    }
    void create_relations_and_delete(const std::string& table_prefix) const
    {
        create_relations(table_prefix);
        delete_from(table_prefix + "_relations");
    }
    void create_generators_x2m_and_delete(const std::string& table_prefix) const
    {
        create_generators_x2m(table_prefix);
        delete_from(table_prefix + "_X2m_generators");
    }
    void create_relations_x2m_and_delete(const std::string& table_prefix) const
    {
        create_relations_x2m(table_prefix);
        delete_from(table_prefix + "_X2m_relations");
    }
    void create_time_and_delete(const std::string& table_prefix) const
    {
        create_time(table_prefix);
        delete_from(table_prefix + "_time");
    }

    void save_generators(const std::string& table_prefix, const DataMRes1d& rels, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (diff, s, t) VALUES (?1, ?2, ?3);");

        for (auto& x : rels) {
            if (x.x2.data.size() == 1 && x.x2.GetLead().deg_m() == 0) {
                stmt.bind_blob(1, x.x1.data);
                stmt.bind_int(2, s + 1);
                stmt.bind_int(3, t);
                stmt.step_and_reset();
            }
        }
    }

    /* insert v_{0, 0} */
    void save_fil_0(const std::string& table_prefix) const
    {
        if (get_int("SELECT MIN(t) from " + table_prefix + "_generators;") != 0) {
            Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (id, diff, s, t) VALUES (?1, ?2, ?3, ?4);");
            stmt.bind_int(1, 0);
            stmt.bind_blob(2, Mod().data);
            stmt.bind_int(3, 0);
            stmt.bind_int(4, 0);
            stmt.step_and_reset();
        }
    }

    void save_relations(const std::string& table_prefix, const DataMRes1d& rels, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (x1, x2, x2m, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");

        for (auto& rel : rels) {
            stmt.bind_blob(1, rel.x1.data);
            stmt.bind_blob(2, rel.x2.data);
            stmt.bind_blob(3, rel.x2m.data);
            stmt.bind_int(4, s);
            stmt.bind_int(5, t);
            stmt.step_and_reset();
        }
    }

    void save_generators_x2m(const std::string& table_prefix, int num, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_X2m_generators (s, t) VALUES (?1, ?2);");

        for (int k = 0; k < num; ++k) {
            stmt.bind_int(1, s);
            stmt.bind_int(2, t);
            stmt.step_and_reset();
        }
    }

    void save_relations_x2m(const std::string& table_prefix, const Mod1d& rels, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_X2m_relations (x2m, s, t) VALUES (?1, ?2, ?3);");

        for (auto& rel : rels) {
            stmt.bind_blob(1, rel.data);
            stmt.bind_int(2, s);
            stmt.bind_int(3, t);
            stmt.step_and_reset();
        }
    }

    void save_time(const std::string& table_prefix, int t, double time)
    {
        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_time (s, t, time) VALUES (?1, ?2, ?3);");

        stmt.bind_int(1, -1);
        stmt.bind_int(2, t);
        stmt.bind_double(3, time);
        stmt.step_and_reset();
    }

    int2d load_basis_degrees(const std::string& table_prefix) const
    {
        int2d basis_degrees;
        Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int s = stmt.column_int(0), t = stmt.column_int(1);
            if (basis_degrees.size() <= (size_t)s)
                basis_degrees.resize(size_t(s + 1));
            basis_degrees[s].push_back(t);
        }
        return basis_degrees;
    }

    int2d load_basis_degrees_x2m(const std::string& table_prefix) const
    {
        return load_basis_degrees(table_prefix + "_X2m");
    }

    DataMRes2d load_data(const std::string& table_prefix) const
    {
        DataMRes2d data;
        Statement stmt(*this, "SELECT x1, x2, x2m, s FROM " + table_prefix + "_relations ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            DataMRes x;
            x.x1.data = stmt.column_blob_tpl<MMod>(0);
            x.x2.data = stmt.column_blob_tpl<MMod>(1);
            x.x2m.data = stmt.column_blob_tpl<MMod>(2);

            x.fil = Filtr(x.x1.GetLead());

            size_t s = (size_t)stmt.column_int(3);
            if (data.size() <= s)
                data.resize(s + 1);
            data[s].push_back(std::move(x));
        }
        return data;
    }

    Mod2d load_data_x2m(const std::string& table_prefix) const
    {
        Mod2d data;
        Statement stmt(*this, "SELECT x2m, s FROM " + table_prefix + "_X2m_relations ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            Mod x;
            x.data = stmt.column_blob_tpl<MMod>(0);
            size_t s = (size_t)stmt.column_int(1);
            if (data.size() <= s)
                data.resize(s + 1);
            data[s].push_back(std::move(x));
        }
        return data;
    }

    std::map<int, int> latest_st(const std::string& table_prefix) const
    {
        std::map<int, int> result;
        Statement stmt(*this, "SELECT s, max(t) FROM " + table_prefix + "_relations GROUP BY s ORDER BY s;");
        while (stmt.step() == MYSQLITE_ROW) {
            int s = stmt.column_int(0);
            int t = stmt.column_int(1);
            result[s] = t;
        }
        return result;
    }

    void save(const DataMRes2d& data, const Mod2d& rels_x2m_cri, const Mod2d& rels_x2m, const int1d& num_x2m, double time, int t)
    {
        if (f_.valid())
            f_.wait();
        f_ = std::async(std::launch::async, [this, data, rels_x2m_cri, rels_x2m, num_x2m, time, t]() {
            begin_transaction();
            for (size_t s = data.size(); s-- > 0;) {
                save_generators("SteenrodMRes", data[s], (int)s, t);
                save_relations("SteenrodMRes", data[s], (int)s, t);
            }
            for (size_t s = num_x2m.size(); s-- > 0;)
                save_generators_x2m("SteenrodMRes", num_x2m[s], (int)s, t);
            for (size_t s = rels_x2m_cri.size(); s-- > 0;) {
                save_relations_x2m("SteenrodMRes", rels_x2m_cri[s], (int)s, t);
                save_relations_x2m("SteenrodMRes", rels_x2m[s], (int)s, t);
            }
            save_time("SteenrodMRes", t, time);
            end_transaction();
        });
    }
};

void Resolve(AdamsRes& gb, const Mod1d& rels, int t_max, int stem_max, const std::string& db_filename)
{
    int t_trunc = gb.t_trunc();
    if (t_max > t_trunc)
        throw MyException(0xb2474e19U, "t_max is bigger than the truncation degree.");
    Mod tmp_Mod;

#ifdef MYDEPLOY
    DbSteenrod db(db_filename);
    db.create_generators("SteenrodMRes");
    db.create_relations("SteenrodMRes");
    db.create_generators_x2m("SteenrodMRes");
    db.create_relations_x2m("SteenrodMRes");
    db.create_time("SteenrodMRes");
    db.save_fil_0("SteenrodMRes");

    bench::Timer timer;
    timer.SuppressPrint();
#endif

    /* Group `rels` by degree */
    std::map<int, int1d> rels_graded;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            const auto& lead = rels[i].GetLead();
            int t = lead.deg_m();
            if (t <= t_max && t <= stem_max)
                rels_graded[t].push_back((int)i);
        }
    }

    for (int t = 1; t <= t_max; ++t) {
        size_t tt = (size_t)t;
        size_t s_min = (size_t)std::max(0, t - stem_max - 2);
        gb.resize_gb(t);
        std::vector<unsigned> old_size_x2m(tt);
        for (size_t s = 0; s < tt; ++s)
            old_size_x2m[s] = (unsigned)gb.basis_degrees_x2m(s).size();

        DataMRes1d data_tmp_neg;
        auto p_rels_d = rels_graded.find(t);
        if (p_rels_d != rels_graded.end()) {
            data_tmp_neg.resize(p_rels_d->second.size());
            auto& rels_d = p_rels_d->second;
            ut::for_each_seq(data_tmp_neg.size(), [&](size_t i) { data_tmp_neg[i] = DataMRes(Mod(), gb.Reduce(rels[rels_d[i]], 0), Mod()); });
        }

        CriMilnor2d cris(tt - 1);
        Mod2d rels_x2m_cri(tt - 1);
        DataMRes2d data_tmps(tt - 1);
        std::vector<unsigned> arr_s;
        for (size_t s = s_min; s < tt - 1; ++s) {
            cris[s] = gb.Criticals(s, t, rels_x2m_cri[s]);
            if (!cris[s].empty()) {
                data_tmps[s].resize(cris[s].size());
                arr_s.push_back((unsigned)s);
            }
        }

        std::atomic<int> threadsLeft = (int)arr_s.size();
        std::mutex print_mutex = {};
        ut::for_each_par(arr_s.size(), [&arr_s, &gb, &cris, &data_tmps, &print_mutex, &threadsLeft, t](size_t i) {
            size_t s = arr_s[i];
            gb.ReduceBatch(cris[s], data_tmps[s], s);

            {
                std::scoped_lock lock(print_mutex);
                --threadsLeft;
                std::cout << "t=" << t << " s=" << s << " threadsLeft=" << threadsLeft << std::endl;
            }
        });

        /* Triangulate these relations */
        DataMRes2d data(tt);
        Mod2d rels_x2m(tt - 1);
        int s_min1 = std::max(0, t - stem_max - 2) - 1; /* This can be -1 */
        for (int s = t - 2; s >= s_min1; --s) {
            size_t ss = (size_t)s; /* When ss is -1 it will not be used */
            size_t sp1 = size_t(s + 1);
            size_t sp2 = size_t(s + 2);

            auto& data_tmp_s = s >= 0 ? data_tmps[s] : data_tmp_neg;

            Mod1d x2m_st_tmp;
            Mod1d kernel_sp1_tmp;
            for (size_t i = 0; i < data_tmp_s.size(); ++i) {
                if (s >= 0) {
                    for (size_t j = 0; j < data[ss].size(); ++j)
                        if (std::binary_search(data_tmp_s[i].x1.data.begin(), data_tmp_s[i].x1.data.end(), data[ss][j].x1.GetLead()))
                            data_tmp_s[i].iaddP(data[ss][j], tmp_Mod);
                    Reduce(data_tmp_s[i].x2, data[sp1], tmp_Mod);
                }
                if (data_tmp_s[i].x1) {
                    /* Determine if x2m aligns with x1 */
                    if (!data_tmp_s[i].valid_x2m()) {
                        x2m_st_tmp.push_back(gb.new_gen_x2m(s, t));
                        std::swap(x2m_st_tmp.back(), data_tmp_s[i].x2m);
                        data_tmp_s[i].fil = Filtr(data_tmp_s[i].x1.GetLead());
                    }
                    data[ss].push_back(std::move(data_tmp_s[i]));
                }
                else {
                    if (data_tmp_s[i].x2)
                        kernel_sp1_tmp.push_back(std::move(data_tmp_s[i].x2));
                    if (data_tmp_s[i].x2m)
                        x2m_st_tmp.push_back(std::move(data_tmp_s[i].x2m));
                }
            }
            for (size_t i = 0; i < kernel_sp1_tmp.size(); ++i) {
                Reduce(kernel_sp1_tmp[i], data[sp1], tmp_Mod);
                if (kernel_sp1_tmp[i])
                    data[sp1].push_back(DataMRes(std::move(kernel_sp1_tmp[i]), gb.new_gen(sp2, t), gb.new_gen_x2m(sp1, t)));
            }
            for (size_t i = 0; i < x2m_st_tmp.size(); ++i) {
                Reduce(x2m_st_tmp[i], rels_x2m[ss], tmp_Mod);
                if (x2m_st_tmp[i])
                    rels_x2m[ss].push_back(std::move(x2m_st_tmp[i]));
            }
        }

        /* Save the result */
#ifdef MYDEPLOY
        double time = timer.Elapsed();
        int1d num_x2m;
        for (size_t s = 0; s < tt; ++s)
            num_x2m.push_back((unsigned)gb.basis_degrees_x2m(s).size() - old_size_x2m[s]);
        db.save(data, rels_x2m_cri, rels_x2m, num_x2m, time, t);
        std::cout << "    time=" << time << std::endl;
        timer.Reset();
#endif

        for (size_t s = tt; s-- > s_min;)
            for (size_t i = 0; i < data[s].size(); ++i)
                gb.push_back(data[s][i], s);
        for (size_t s = tt - 1; s-- > s_min;)
            for (size_t i = 0; i < rels_x2m[s].size(); ++i)
                gb.push_back_x2m(rels_x2m[s][i], s);
    }
}

AdamsRes AdamsRes::load(const std::string& db_filename, const std::string& tablename, int t_trunc, int stem_trunc)
{
    DbSteenrod db(db_filename);
    db.create_generators(tablename);
    db.create_relations(tablename);
    db.create_generators_x2m(tablename);
    db.create_relations_x2m(tablename);
    DataMRes2d data = db.load_data(tablename);
    int2d basis_degrees = db.load_basis_degrees(tablename);
    Mod2d data_x2m = db.load_data_x2m(tablename);
    int2d basis_degrees_x2m = db.load_basis_degrees_x2m(tablename);
    auto latest_st = db.latest_st(tablename);
    return AdamsRes(t_trunc, stem_trunc, std::move(data), std::move(basis_degrees), std::move(data_x2m), std::move(basis_degrees_x2m), latest_st);
}

void ResetDb(const std::string& filename, const std::string& tablename)
{
    DbSteenrod db(filename);
    db.create_generators_and_delete(tablename);
    db.create_relations_and_delete(tablename);
    db.create_generators_x2m_and_delete(tablename);
    db.create_relations_x2m_and_delete(tablename);
    db.create_time_and_delete(tablename);
}

}  // namespace steenrod