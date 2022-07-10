#include "groebner.h"
#include "myio.h"
#include "utility.h"

namespace alg {

bool detail::HasGCD(const Mon& mon1, const Mon& mon2)
{
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            k++;
        else if (k->g_raw() < l->g_raw())
            l++;
        else
            return true;
    }
    return false;
}

int detail::DegLCM(const Mon& mon1, const Mon& mon2, const int1d& gen_degs)
{
    int result = 0;
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw()) {
            result += gen_degs[k->g()] * k->e();
            ++k;
        }
        else if (k->g_raw() < l->g_raw()) {
            result += gen_degs[l->g()] * l->e();
            ++l;
        }
        else {
            if (k->e() < l->e())
                result += gen_degs[l->g()] * l->e();
            else
                result += gen_degs[k->g()] * k->e();
            ++k;
            ++l;
        }
    }
    for (; k != mon1.end(); ++k)
        result += gen_degs[k->g()] * k->e();
    for (; l != mon2.end(); ++l)
        result += gen_degs[l->g()] * l->e();
    return result;
}

void CriPair::SetFromLM(CriPair& result, const Mon& lead1, const Mon& lead2, int i, int j)
{
    auto k = lead1.begin(), l = lead2.begin();
    while (k != lead1.end() && l != lead2.end()) {
        if (k->g_raw() > l->g_raw())
            result.m2.push_back(*k++);
        else if (k->g_raw() < l->g_raw())
            result.m1.push_back(*l++);
        else {
            if (k->e() < l->e())
                result.m1.push_back(GE(l->data - k->e()));
            else if (k->e() > l->e())
                result.m2.push_back(GE(k->data - l->e()));
            k++;
            l++;
        }
    }
    if (k != lead1.end())
        result.m2.insert(k, lead1.end());
    else
        result.m1.insert(l, lead2.end());
    result.i1 = i;
    result.i2 = j;
    result.trace_m2 = result.m2.Trace();
}

void CriPair::SijP(const Groebner& gb, Poly& result, Poly& tmp1, Poly& tmp2) const
{
    result.data.clear();
    mulP(gb[i1], m1, tmp1);
    result.iaddP(tmp1, tmp2);
    mulP(gb[i2], m2, tmp1);
    result.iaddP(tmp1, tmp2);
    return;
}

void GbCriPairs::Minimize(const Mon1d& leads, int d)
{
    if (buffer_redundent_pairs_.empty() || buffer_redundent_pairs_.begin()->first != d)
        return;
    if (buffer_min_pairs_.empty() || buffer_min_pairs_.begin()->first != d)
#ifdef NDEBUG
        return;
#else
        buffer_min_pairs_[d];
#endif
    /* Add to the Groebner basis of critical pairs */
    constexpr uint64_t NULL_J = ~uint64_t(0);
    auto& b_min_pairs_d = buffer_min_pairs_.begin()->second;
    for (uint64_t ij : buffer_redundent_pairs_.begin()->second) {
        uint64_t i, j;
        ut::UnBind(ij, i, j);
        while (j != NULL_J) {
            Mon gcd = GCD(leads[i], leads[j]);
            Mon m2 = leads[i] / gcd;
            if (j < (int)b_min_pairs_d.size()) {
                auto p = std::find_if(b_min_pairs_d[j].begin(), b_min_pairs_d[j].end(), [&m2](const CriPair& c) { return c.m2 == m2; });
                /* Remove it from `buffer_min_pairs_` */
                if (p != b_min_pairs_d[j].end()) {
                    p->i2 = -1;
                    break;
                }
            }

            /* Reduce (i, j) */
            MonTrace t_m2 = m2.Trace();
            auto c = gb_[j].begin();
            auto end = gb_[j].end();
            for (; c < end; ++c) {
                if (divisible(c->m2, m2, c->trace_m2, t_m2)) {
                    Mon m1 = leads[j] / gcd;
                    if (detail::HasGCD(c->m1, m1)) {
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

    /* Delete `bufbuffer_redundent_pairs_[d]` */
    buffer_redundent_pairs_.erase(d);
}

void GbCriPairs::AddToBuffers(const Mon1d& leads, const MonTrace1d& traces, const Mon& mon, const int1d& gen_degs)
{
    size_t lead_size = leads.size();
    if (gb_.size() < lead_size + 1)
        gb_.resize(lead_size + 1);
    std::vector<std::pair<int, CriPair>> new_pairs(leads.size());
    MonTrace t = mon.Trace();

    /* Populate `new_pairs` */
    for (size_t i = 0; i < lead_size; ++i) {
        new_pairs[i].first = -1;
        if (detail::HasGCD(leads[i], mon, traces[i], t) /* && pred */) {
            int d_pair = detail::DegLCM(leads[i], mon, gen_degs);
            if (d_pair <= deg_trunc_) {
                new_pairs[i].first = d_pair;
                CriPair::SetFromLM(new_pairs[i].second, leads[i], mon, (int)i, (int)lead_size);
            }
        }
    }

    /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
    for (size_t j = 1; j < new_pairs.size(); ++j) {
        if (new_pairs[j].first != -1) {
            for (size_t i = 0; i < j; ++i) {
                if (new_pairs[i].first != -1) {
                    if (divisible(new_pairs[i].second.m2, new_pairs[j].second.m2, new_pairs[i].second.trace_m2, new_pairs[j].second.trace_m2)) {
                        new_pairs[j].first = -1;
                        break;
                    }
                    else if (divisible(new_pairs[j].second.m2, new_pairs[i].second.m2, new_pairs[j].second.trace_m2, new_pairs[i].second.trace_m2))
                        new_pairs[i].first = -1;
                    else if (!detail::HasGCD(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                        int dij = detail::DegLCM(leads[i], leads[j], gen_degs);
                        if (dij <= deg_trunc_)
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
            buffer_min_pairs_[new_pairs[i].first][lead_size].push_back(std::move(new_pairs[i].second));
        }
    }
}

void GbCriPairs::init(const Mon1d& leads, const MonTrace1d& traces, const int1d& gen_degs, int t_min_buffer)
{
    Mon1d tmp_leads;
    for (auto& mon : leads) {
        AddToBuffers(tmp_leads, traces, mon, gen_degs);
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
}

Groebner::Groebner(int deg_trunc, int1d gen_degs, Poly1d polys, bool bDynamic) : criticals_(deg_trunc), data_(std::move(polys)), gen_degs_(std::move(gen_degs))
{
    for (int i = 0; i < (int)data_.size(); ++i) {
        leads_.push_back(data_[i].GetLead());
        traces_.push_back(data_[i].GetLead().Trace());
        indices_[Key(data_[i].GetLead())].push_back(i);
    }
    if (bDynamic)
        criticals_.init(leads_, traces_, gen_degs_, deg_trunc + 1);
}

int Groebner::IndexOfDivisibleLeading(const Mon& mon) const
{
    auto t = mon.Trace();
    for (int i = 0; i < (int)mon.size(); ++i) {
        for (int j = -1; j < i; ++j) {
            auto key = TypeIndexKey{mon[i].g() + (j == -1 ? 0 : ((mon[j].g() + 1) << 16))};
            auto p = indices_.find(key);
            if (p != indices_.end()) {
                for (int k : p->second) {
                    if (divisible(leads_[k], mon, traces_[k], t))
                        return k;
                }
            }
        }
    }
    return -1;
}

Poly Groebner::Reduce(Poly poly) const
{
    Poly tmp_prod, tmp;
    size_t index = 0;
    while (index < poly.data.size()) {
        int gb_index = IndexOfDivisibleLeading(poly.data[index]);
        if (gb_index != -1) {
            Mon q = poly.data[index] / data_[gb_index].GetLead();
            mulP(data_[gb_index], q, tmp_prod);
            poly.iaddP(tmp_prod, tmp);
        }
        else
            ++index;
    }
    return poly;
}

void Groebner::AddRels(const Poly1d& rels, int deg_max)
{
    using Poly = Poly;
    using Poly1d = std::vector<Poly>;
    using PPoly1d = std::vector<const Poly*>;
    Poly tmp1, tmp2;

    int d_trunc = deg_trunc();
    if (deg_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, PPoly1d> rels_graded;
    for (const auto& rel : rels) {
        if (rel) {
            int d = GetDeg(rel.GetLead(), gen_degs_);
            if (d <= deg_max)
                rels_graded[d].push_back(&rel);
        }
    }

    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg_max && (d <= deg_max_rels || !criticals_.empty()); ++d) {
        CriPair1d pairs_d = Criticals(d);
        auto& rels_d = rels_graded[d];
        Poly1d rels_tmp(pairs_d.size() + rels_d.size());

        for (size_t i = 0; i < pairs_d.size(); ++i)
            pairs_d[i].SijP(*this, rels_tmp[i], tmp1, tmp2);
        size_t pairs_d_size = pairs_d.size();
        ut::for_each_seq(pairs_d.size(), [this, &pairs_d, &rels_tmp](size_t i) { rels_tmp[i] = Reduce(std::move(rels_tmp[i])); });
        ut::for_each_seq(rels_d.size(), [this, &rels_d, &rels_tmp, pairs_d_size](size_t i) { rels_tmp[pairs_d_size + i] = Reduce(*rels_d[i]); });

        /* Triangulate these relations */
        Poly1d gb_rels_d;
        for (auto& rel : rels_tmp) {
            for (const Poly& rel1 : gb_rels_d)
                if (std::binary_search(rel.data.begin(), rel.data.end(), rel1.GetLead()))
                    rel.iaddP(rel1, tmp1);
            if (rel)
                gb_rels_d.push_back(std::move(rel));
        }

        /* Add these relations */
        for (auto& rel : gb_rels_d)
            push_back(std::move(rel));
    }
}

void powP(const Poly& poly, int n, const Groebner& gb, Poly& result, Poly& tmp)
{
    result.data.clear();
    result.data.push_back(Mon());
    if (n == 0)
        return;
    Poly power = poly;
    while (n) {
        if (n & 1) {
            result.imulP(power, tmp);
            result = gb.Reduce(std::move(result));
        }
        n >>= 1;
        if (n) {
            power.frobP(tmp);
            power = gb.Reduce(tmp);
        }
    }
}

} /* namespace alg */