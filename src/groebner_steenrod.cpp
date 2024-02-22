#include "groebner_steenrod.h"
#include "benchmark.h"
#include "database.h"
#include <atomic>
#include <cstring>
#include <mutex>

/********************************************************
 *                    class CriMilnors
 ********************************************************/

namespace steenrod {

void CriMilnor::Sij(const Groebner& gb, Mod& result, Milnor& tmp_a, Mod& tmp1, Mod& tmp2) const
{
    if (i1 >= 0) {
        result.iaddmulP(m1, gb.data()[i1], tmp_a, tmp1, tmp2).iaddmulP(m2, gb.data()[i2], tmp_a, tmp1, tmp2);
    }
    else {
        result.iaddmulP(m2, gb.data()[i2], tmp_a, tmp1, tmp2);
    }
}

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
 *                    class Groebner
 ********************************************************/
Groebner::Groebner(int d_trunc, Mod1d data, int1d v_degs) : d_trunc_(d_trunc), criticals_(d_trunc), gb_(std::move(data)), v_degs_(std::move(v_degs))
{
    for (int j = 0; j < (int)gb_.size(); ++j) {
        leads_.push_back(gb_[j].GetLead());
        indices_[gb_[j].GetLead().v_raw()].push_back(j);
    }
    criticals_.init(leads_, v_degs_, 0);
}

CriMilnor1d Groebner::Criticals(int t)
{
    criticals_.Minimize(leads_, t);
    return criticals_.Criticals(t);
}

inline int IndexOfDivisibleLeading(const MMod1d& leads, const std::unordered_map<uint64_t, int1d>& indices, MMod mon)
{
    auto key = mon.v_raw();
    auto p = indices.find(key);
    if (p != indices.end())
        for (int k : p->second)
            if (divisibleLF(leads[k], mon))
                return k;
    return -1;
}

Mod Groebner::Reduce(const CriMilnor& cp) const
{
    Mod result;

    Milnor tmp_a;
    Mod tmp_x1, tmp_x2;
    tmp_a.data.reserve(32);
    tmp_x1.data.reserve(64);
    tmp_x2.data.reserve(64);

    if (cp.i1 >= 0) {
        result.iaddmulP(cp.m1, gb_[cp.i1], tmp_a, tmp_x1, tmp_x2).iaddmulP(cp.m2, gb_[cp.i2], tmp_a, tmp_x1, tmp_x2);
    }
    else {
        result.iaddmulP(cp.m2, gb_[cp.i2], tmp_a, tmp_x1, tmp_x2);
    }

    size_t index;
    index = 0;
    while (index < result.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_, indices_, result.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.data[index], gb_[gb_index].data[0]);
            result.iaddmulP(m, gb_[gb_index], tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }
    return result;
}

Mod Groebner::Reduce(Mod x) const
{
    size_t index;
    index = 0;
    Milnor tmp_a;
    Mod tmp_x1, tmp_x2;
    while (index < x.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_, indices_, x.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(x.data[index], gb_[gb_index].data[0]);
            x.iaddmulP(m, gb_[gb_index], tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }

    return x;
}

void Groebner::AddRels(const Mod1d& rels, int deg_max, int1d& min_rels)
{
    using PMod1d = std::vector<const Mod*>;
    Milnor tmp_a;
    Mod tmp1, tmp2;

    if (deg_max > d_trunc_)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, PMod1d> rels_graded;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            const auto& lead = rels[i].GetLead();
            int d = lead.deg_m() + v_degs_[lead.v()];
            if (d <= deg_max)
                rels_graded[d].push_back(&rels[i]);
        }
    }

    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg_max && (d <= deg_max_rels || !criticals_.empty()); ++d) {
        auto pairs_d = Criticals(d);
        auto& rels_d = rels_graded[d];
        Mod1d rels_tmp(pairs_d.size() + rels_d.size());

        for (size_t i = 0; i < pairs_d.size(); ++i)
            pairs_d[i].Sij(*this, rels_tmp[i], tmp_a, tmp1, tmp2);
        size_t pairs_d_size = pairs_d.size();
        ut::for_each_seq(pairs_d.size(), [this, &pairs_d, &rels_tmp](size_t i) { rels_tmp[i] = Reduce(std::move(rels_tmp[i])); });
        ut::for_each_seq(rels_d.size(), [this, &rels_d, &rels_tmp, pairs_d_size](size_t i) { rels_tmp[pairs_d_size + i] = Reduce(*rels_d[i]); });

        /* Triangulate these relations */
        Mod1d gb_rels_d;
        for (size_t i = 0; i < rels_tmp.size(); ++i) {
            auto& rel = rels_tmp[i];
            for (const Mod& rel1 : gb_rels_d)
                if (std::binary_search(rel.data.begin(), rel.data.end(), rel1.GetLead()))
                    rel.iaddP(rel1, tmp1);
            if (rel) {
                if (i >= pairs_d_size)
                    min_rels.push_back(int(gb_.size() + gb_rels_d.size()));
                gb_rels_d.push_back(std::move(rel));
            }
        }

        /* Add these relations */
        for (auto& rel : gb_rels_d)
            push_back(std::move(rel));
    }
}

void Groebner::AddRels(const Mod1d& rels, int deg_max)
{
    int1d min_rels;
    AddRels(rels, deg_max, min_rels);
}

Mod subV(const Mod& x, const int1d& v_map)
{
    Mod result;
    for (auto& m : x.data)
        result.data.push_back(MMod(m.m(), (uint64_t)v_map[m.v()]));
    return result;
}

void Groebner::MinimizeOrderedGensRels(Mod1d& cells, int1d& min_rels)
{
    std::vector<size_t> redundant_vs;
    Mod1d data;
    std::map<int, int> map_rel_ind;
    cells.clear();
    for (size_t i = 0; i < gb_.size(); ++i) {
        if (!gb_[i].GetLead().m()) {
            redundant_vs.push_back(gb_[i].GetLead().v());
            ut::get(cells, gb_[i].GetLead().v()) = gb_[i] + gb_[i].GetLead();
        }
        else {
            map_rel_ind[(int)i] = (int)data.size();
            data.push_back(gb_[i]);
        }
    }

    int1d min_rels_new;
    for (int& i : min_rels)
        if (ut::has(map_rel_ind, i))
            min_rels_new.push_back(map_rel_ind.at(i));
    min_rels = std::move(min_rels_new);

    std::sort(redundant_vs.begin(), redundant_vs.end());
    auto range = ut::size_t_range(v_degs_.size());
    std::vector<size_t> remaining_v;
    std::set_symmetric_difference(redundant_vs.begin(), redundant_vs.end(), range.begin(), range.end(), std::back_inserter(remaining_v));

    int1d v_map(v_degs_.size(), -1);
    for (size_t i = 0; i < remaining_v.size(); ++i) {
        v_map[remaining_v[i]] = (int)i;
        ut::get(cells, remaining_v[i]) = MMod(MMilnor(), remaining_v[i]);
    }

    for (auto& rel : data)
        rel = subV(rel, v_map);
    for (auto& rel : cells)
        rel = subV(rel, v_map);

    int1d new_v_degs;
    for (auto& i : remaining_v)
        new_v_degs.push_back(v_degs_[i]);

    ResetRels();
    v_degs_ = std::move(new_v_degs);
    for (auto& rel : data)
        push_back(rel);
}

}  // namespace steenrod