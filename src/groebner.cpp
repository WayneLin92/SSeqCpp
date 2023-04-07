#include "groebner.h"
#include "myio.h"
#include "utility.h"

namespace alg2 {

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

template <typename FnGenDeg>
auto TplDegLCM(const Mon& mon1, const Mon& mon2, FnGenDeg _gen_deg)
{
    using ReturnType = decltype(_gen_deg(0));
    auto result = ReturnType();
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw()) {
            result += _gen_deg(k->g()) * k->e_masked();
            ++k;
        }
        else if (k->g_raw() < l->g_raw()) {
            result += _gen_deg(l->g()) * l->e_masked();
            ++l;
        }
        else {
            if (k->e() < l->e())
                result += _gen_deg(l->g()) * l->e_masked();
            else
                result += _gen_deg(k->g()) * k->e_masked();
            ++k;
            ++l;
        }
    }
    for (; k != mon1.end(); ++k)
        result += _gen_deg(k->g()) * k->e_masked();
    for (; l != mon2.end(); ++l)
        result += _gen_deg(l->g()) * l->e_masked();
    return result;
}

int detail::DegLCM(const Mon& mon1, const Mon& mon2, const int1d& gen_degs)
{
    return TplDegLCM(mon1, mon2, [&gen_degs](size_t i) { return gen_degs[i]; });
}

int detail::DegTLCM(const Mon& mon1, const Mon& mon2, const AdamsDeg1d& gen_degs)
{
    return TplDegLCM(mon1, mon2, [&gen_degs](size_t i) { return gen_degs[i].t; });
}

void detail::MutualQuotient(Mon& m1, Mon& m2, const Mon& lead1, const Mon& lead2)
{
    auto k = lead1.begin(), l = lead2.begin();
    while (k != lead1.end() && l != lead2.end()) {
        if (k->g_raw() > l->g_raw())
            m2.push_back(*k++);
        else if (k->g_raw() < l->g_raw())
            m1.push_back(*l++);
        else {
            if (k->e() < l->e())
                m1.push_back(GE(l->data - k->e_masked()));
            else if (k->e() > l->e())
                m2.push_back(GE(k->data - l->e_masked()));
            k++;
            l++;
        }
    }
    if (k != lead1.end())
        m2.insert(k, lead1.end());
    else
        m1.insert(l, lead2.end());
}

void CriPair::SetFromLM(CriPair& result, const Mon& lead1, const Mon& lead2, int i, int j)
{
    detail::MutualQuotient(result.m1, result.m2, lead1, lead2);
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
}

void CriPair::SijMod(const Groebner& gb, const GroebnerMod& gbm, Mod& result, Mod& tmp1, Mod& tmp2) const
{
    result.data.clear();
    if (i1 & FLAG_INDEX_X)
        tmp1 = Mod(gb[i1 ^ FLAG_INDEX_X] * m1, gbm[i2].GetLead().v);
    else
        mulP(m1, gbm[i1], tmp1);
    result.iaddP(tmp1, tmp2);
    mulP(m2, gbm[i2], tmp1);
    result.iaddP(tmp1, tmp2);
    return;
}

void GbCriPairs::Minimize(const Mon1d& leads, int d)
{
    if (buffer_redundent_pairs_.empty() || buffer_redundent_pairs_.begin()->first != d)
        return;
    if (buffer_min_pairs_.empty() || buffer_min_pairs_.begin()->first != d) {
#ifdef NDEBUG
        buffer_redundent_pairs_.erase(d);
        return;
#else
        buffer_min_pairs_[d];
#endif
    }

    /* Add to the Groebner basis of critical pairs */
    constexpr uint64_t NULL_J = ~uint64_t(0);
    auto& b_min_pairs_d = buffer_min_pairs_.begin()->second;
    for (uint64_t ij : buffer_redundent_pairs_.begin()->second) {
        uint64_t i, j;
        ut::UnBind(ij, i, j);
        while (j != NULL_J) {
            Mon gcd = GCD(leads[i], leads[j]);
            if (!gcd)
                break;
            Mon m2 = leads[i] / gcd;
            if (j < (int)b_min_pairs_d.size()) {
                auto p = std::find_if(b_min_pairs_d[j].begin(), b_min_pairs_d[j].end(), [&m2](const CriPair& c) { return c.m2 == m2; });
                /* Remove it from `buffer_min_pairs_` */
                if (p != b_min_pairs_d[j].end()) {
                    p->i2 = NULL_INDEX32;
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

void GbCriPairs::Minimize(const Mon1d& leadsx, const MMod1d& leads, int d)
{
    if (buffer_redundent_pairs_.empty() || buffer_redundent_pairs_.begin()->first != d)
        return;
    if (buffer_min_pairs_.empty() || buffer_min_pairs_.begin()->first != d) {
#ifdef NDEBUG
        buffer_redundent_pairs_.erase(d);
        return;
#else
        buffer_min_pairs_[d];
#endif
    }

    /* Add to the Groebner basis of critical pairs */
    constexpr uint64_t NULL_J = ~uint64_t(0);
    auto& b_min_pairs_d = buffer_min_pairs_.begin()->second;
    for (uint64_t ij : buffer_redundent_pairs_.begin()->second) {
        uint64_t i, j;
        ut::UnBind(ij, i, j);
        while (j != NULL_J) {
            Mon leadsi = (i & FLAG_INDEX_X) ? leadsx[i ^ FLAG_INDEX_X] : leads[i].m;
            Mon gcd = GCD(leadsi, leads[j].m);
            Mon m2 = leadsi / gcd;
            if (j < (int)b_min_pairs_d.size()) {
                auto p = std::find_if(b_min_pairs_d[j].begin(), b_min_pairs_d[j].end(), [&m2](const CriPair& c) { return c.m2 == m2; });
                /* Remove it from `buffer_min_pairs_` */
                if (p != b_min_pairs_d[j].end()) {
                    p->i2 = NULL_INDEX32;
                    break;
                }
            }

            /* Reduce (i, j) */
            MonTrace t_m2 = m2.Trace();
            auto c = gb_[j].begin();
            auto end = gb_[j].end();
            for (; c < end; ++c) {
                if (divisible(c->m2, m2, c->trace_m2, t_m2)) {
                    Mon m1 = leads[j].m / gcd;
                    if (detail::HasGCD(c->m1, m1)) {
                        j = NULL_J;
                        break;
                    }
                    else {
                        j = c->i1;
                        if (i > j || (j & FLAG_INDEX_X))
                            std::swap(i, j);
                        if (j & FLAG_INDEX_X)
                            j = NULL_J;
                        break;
                    }
                }
            }
#ifndef NDEBUG
            if (c == end)
                throw MyException(0xbcf149e1U, "Should not happen because gb_ is groebner");
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
        if (detail::HasGCD(leads[i], mon, traces[i], t)) {
            int d_pair = detail::DegLCM(leads[i], mon, gen_degs);
            if (d_pair <= deg_trunc_) {
                new_pairs[i].first = d_pair;
                CriPair::SetFromLM(new_pairs[i].second, leads[i], mon, (uint32_t)i, (uint32_t)lead_size);
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
                    else if (!detail::HasGCD(new_pairs[i].second.m1, new_pairs[j].second.m1) && detail::HasGCD(leads[i], leads[j], traces[i], traces[j])) {
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

void GbCriPairs::AddToBuffers(const Mon1d& leadsx, const MonTrace1d& tracesx, const MMod1d& leads, const MonTrace1d& traces, const MMod& mon, const int1d& gen_degs, const int1d& v_degs)
{
    size_t leadsx_size = leadsx.size();
    size_t leads_size = leads.size();
    if (gb_.size() < leads_size + 1)
        gb_.resize(leads_size + 1);
    std::vector<std::pair<int, CriPair>> new_pairs(leadsx_size + leads_size);
    MonTrace t = mon.m.Trace();

    /* Populate `new_pairs` */
    for (size_t i = 0; i < leadsx_size; ++i) {
        new_pairs[i].first = -1;
        int d_pair = detail::DegLCM(leadsx[i], mon.m, gen_degs) + v_degs[mon.v];
        if (d_pair <= deg_trunc_) {
            new_pairs[i].first = d_pair;
            CriPair::SetFromLM(new_pairs[i].second, leadsx[i], mon.m, uint32_t(i | FLAG_INDEX_X), (uint32_t)leads_size);
        }
    }
    for (size_t i = 0; i < leads_size; ++i) {
        auto& new_pairs_i = new_pairs[leadsx_size + i];
        new_pairs_i.first = -1;
        if (leads[i].v == mon.v) {
            int d_pair = detail::DegLCM(leads[i].m, mon.m, gen_degs) + v_degs[mon.v];
            if (d_pair <= deg_trunc_) {
                new_pairs_i.first = d_pair;
                CriPair::SetFromLM(new_pairs_i.second, leads[i].m, mon.m, (uint32_t)i, (uint32_t)leads_size);
            }
        }
    }

    /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
    for (size_t j = 0; j < new_pairs.size(); ++j) {
        if (new_pairs[j].first != -1) {
            for (size_t i = 0; i < j; ++i) {
                if (new_pairs[i].first != -1) {
                    if (divisible(new_pairs[i].second.m2, new_pairs[j].second.m2, new_pairs[i].second.trace_m2, new_pairs[j].second.trace_m2)) {
                        new_pairs[j].first = -1;
                        break;
                    }
                    else if (divisible(new_pairs[j].second.m2, new_pairs[i].second.m2, new_pairs[j].second.trace_m2, new_pairs[i].second.trace_m2))
                        new_pairs[i].first = -1;
                    else if (!detail::HasGCD(new_pairs[i].second.m1, new_pairs[j].second.m1) && j >= leadsx_size) {
                        Mon leadsi = i < leadsx_size ? leadsx[i] : leads[i - leadsx_size].m;
                        const size_t j1 = j - leadsx_size;
                        int dij = detail::DegLCM(leadsi, leads[j1].m, gen_degs) + v_degs[leads[j1].v];
                        if (dij <= deg_trunc_)
                            buffer_redundent_pairs_[dij].insert(ut::Bind(i < leadsx_size ? (i | FLAG_INDEX_X) : i - leadsx_size, j1));
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < new_pairs.size(); ++i) {
        if (new_pairs[i].first != -1) {
            gb_[leads_size].push_back(new_pairs[i].second);
            buffer_min_pairs_[new_pairs[i].first].resize(leads_size + 1);
            buffer_min_pairs_[new_pairs[i].first][leads_size].push_back(std::move(new_pairs[i].second));
        }
    }
}

template <typename T>
void RemoveSmallKey(T& map, int t_min_buffer)
{
    for (auto it = map.begin(); it != map.end();) {
        if (it->first < t_min_buffer)
            it = map.erase(it);
        else
            ++it;
    }
}

void GbCriPairs::init(const Mon1d& leads, const MonTrace1d& traces, const int1d& gen_degs, int t_min_buffer)
{
    Mon1d tmp_leads;
    for (auto& mon : leads) {
        AddToBuffers(tmp_leads, traces, mon, gen_degs);
        tmp_leads.push_back(mon);
    }
    RemoveSmallKey(buffer_min_pairs_, t_min_buffer);
    RemoveSmallKey(buffer_redundent_pairs_, t_min_buffer);
}

void GbCriPairs::init(const Mon1d& leadsx, const MonTrace1d& tracesx, const MMod1d& leads, const MonTrace1d& traces, const int1d& gen_degs, const int1d& v_degs, int t_min_buffer)
{
    MMod1d tmp_leads;
    for (auto& mon : leads) {
        AddToBuffers(leadsx, tracesx, tmp_leads, traces, mon, gen_degs, v_degs);
        tmp_leads.push_back(mon);
    }
    RemoveSmallKey(buffer_min_pairs_, t_min_buffer);
    RemoveSmallKey(buffer_redundent_pairs_, t_min_buffer);
}

Groebner::Groebner(int deg_trunc, int1d gen_degs, Poly1d data, bool bDynamic) : criticals_(deg_trunc), data_(std::move(data)), gen_degs_(std::move(gen_degs))
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

void Groebner::AddRels(const Poly1d& rels, int deg_max, int1d& min_rels)
{
    using PPoly1d = std::vector<const Poly*>;
    Poly tmp1, tmp2;

    int d_trunc = deg_trunc();
    if (deg_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, PPoly1d> rels_graded;
    std::map<int, int1d> indices_rels;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            int d = GetDeg(rels[i].GetLead(), gen_degs_);
            if (d <= deg_max) {
                rels_graded[d].push_back(&rels[i]);
                indices_rels[d].push_back((int)i);
            }
        }
    }

    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg_max && (d <= deg_max_rels || !criticals_.empty()); ++d) {
        CriPair1d pairs_d = Criticals(d);
        size_t pairs_d_size = pairs_d.size();
        auto& rels_d = rels_graded[d];
        Poly1d rels_tmp(pairs_d.size() + rels_d.size());

        for (size_t i = 0; i < pairs_d.size(); ++i)
            pairs_d[i].SijP(*this, rels_tmp[i], tmp1, tmp2);
        ut::for_each_seq(pairs_d.size(), [this, &pairs_d, &rels_tmp](size_t i) { rels_tmp[i] = Reduce(std::move(rels_tmp[i])); });
        ut::for_each_seq(rels_d.size(), [this, &rels_d, &rels_tmp, pairs_d_size](size_t i) { rels_tmp[pairs_d_size + i] = Reduce(*rels_d[i]); });

        /* Triangulate these relations */
        Poly1d gb_rels_d;
        for (size_t i = 0; i < rels_tmp.size(); ++i) {
            auto& rel = rels_tmp[i];
            for (const Poly& rel1 : gb_rels_d)
                if (std::binary_search(rel.data.begin(), rel.data.end(), rel1.GetLead()))
                    rel.iaddP(rel1, tmp1);
            if (rel) {
                gb_rels_d.push_back(std::move(rel));
                if (i >= pairs_d_size)
                    min_rels.push_back(indices_rels.at(d)[i - pairs_d_size]);
            }
        }

        /* Add these relations */
        for (auto& rel : gb_rels_d)
            push_back(std::move(rel));
    }
}

void Groebner::AddRels(const Poly1d& rels, int deg) {
    int1d min_rels;
    AddRels(rels, deg, min_rels);
}

void Groebner::ReducedGb()
{
    for (auto& rel : data_) {
        auto rel1 = rel + rel.GetLead();
        rel1 = Reduce(std::move(rel1));
        rel = rel1 + rel.GetLead();
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

/********************************* Modules ****************************************/
GroebnerMod::GroebnerMod(Groebner* pGb, int deg_trunc, int1d v_degs, Mod1d polys, bool bDynamic) : pGb_(pGb), criticals_(deg_trunc), data_(std::move(polys)), v_degs_(std::move(v_degs))
{
    for (int i = 0; i < (int)data_.size(); ++i) {
        leads_.push_back(data_[i].GetLead());
        traces_.push_back(data_[i].GetLead().m.Trace());
        indices_[Key(data_[i].GetLead())].push_back(i);
    }
    if (bDynamic)
        criticals_.init(pGb_->leads_, pGb_->traces_, leads_, traces_, pGb_->gen_degs(), v_degs_, deg_trunc + 1);
}

int GroebnerMod::IndexOfDivisibleLeading(const MMod& mon) const
{
    auto t = mon.m.Trace();
    for (int i = -1; i < (int)mon.m.size(); ++i) {
        auto key = TypeIndexKey{mon.v + (i == -1 ? 0 : ((mon.m[i].g() + 1) << 16))};
        auto p = indices_.find(key);
        if (p != indices_.end()) {
            for (int k : p->second) {
                if (divisible(leads_[k].m, mon.m, traces_[k], t))
                    return k;
            }
        }
    }
    return -1;
}

Mod GroebnerMod::Reduce(Mod poly) const
{
    Mod tmp_prod, tmp;
    Poly tmp_prod1;
    size_t index = 0;
    while (index < poly.data.size()) {
        int gbmod_index = IndexOfDivisibleLeading(poly.data[index]);
        if (gbmod_index != -1) {
            Mon q = poly.data[index].m / data_[gbmod_index].GetLead().m;
            mulP(q, data_[gbmod_index], tmp_prod);
            poly.iaddP(tmp_prod, tmp);
        }
        else {
            int gb_index = pGb_->IndexOfDivisibleLeading(poly.data[index].m);
            if (gb_index != -1) {
                Mon q = poly.data[index].m / pGb_->data()[gb_index].GetLead();
                mulP(pGb_->data()[gb_index], q, tmp_prod1);
                tmp_prod = Mod(tmp_prod1, poly.data[index].v);
                poly.iaddP(tmp_prod, tmp);
            }
            else
                ++index;
        }
    }
    return poly;
}

void GroebnerMod::AddRels(const Mod1d& rels, int deg_max)
{
    using PMod1d = std::vector<const Mod*>;
    Mod tmp1, tmp2;

    int d_trunc = deg_trunc();
    if (deg_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, PMod1d> rels_graded;
    for (const auto& rel : rels) {
        if (rel) {
            int d = GetDeg(rel.GetLead().m, pGb_->gen_degs()) + v_degs_[rel.GetLead().v];
            if (d <= deg_max)
                rels_graded[d].push_back(&rel);
        }
    }

    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg_max && (d <= deg_max_rels || !criticals_.empty()); ++d) {
        CriPair1d pairs_d = Criticals(d);
        auto& rels_d = rels_graded[d];
        Mod1d rels_tmp(pairs_d.size() + rels_d.size());

        for (size_t i = 0; i < pairs_d.size(); ++i)
            pairs_d[i].SijMod(*pGb_, *this, rels_tmp[i], tmp1, tmp2);
        size_t pairs_d_size = pairs_d.size();
        ut::for_each_seq(pairs_d.size(), [this, &pairs_d, &rels_tmp](size_t i) { rels_tmp[i] = Reduce(std::move(rels_tmp[i])); });
        ut::for_each_seq(rels_d.size(), [this, &rels_d, &rels_tmp, pairs_d_size](size_t i) { rels_tmp[pairs_d_size + i] = Reduce(*rels_d[i]); });

        /* Triangulate these relations */
        Mod1d gb_rels_d;
        for (auto& rel : rels_tmp) {
            for (const Mod& rel1 : gb_rels_d)
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

void GroebnerMod::ToSubMod(const Mod1d& rels, int deg_max, int1d& index_ind)
{
    int d_trunc = deg_trunc();
    if (deg_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");

    Mod tmp1, tmp2;

    const size_t rels_size = rels.size();
    for (Mod& x : data_)
        for (MMod m : x.data)
            m.v += (uint32_t)rels_size;
    v_degs_.resize(v_degs_.size() + rels_size);
    for (size_t i = v_degs_.size(); i-- > rels_size;)
        v_degs_[i] = v_degs_[i - rels_size];

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, Mod1d> rels_graded;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            int d = GetDeg(rels[i].GetLead().m, pGb_->gen_degs()) + v_degs_[(size_t)rels[i].GetLead().v + rels_size];
            if (!rels_graded.empty() && d < rels_graded.rbegin()->first)
                throw MyException(0x39b4cfd4U, "rels is not ordered.");
            if (d <= deg_max) {
                Mod rel = rels[i];
                v_degs_[i] = d;
                for (MMod& m : rel.data)
                    m.v += (uint32_t)rels_size;
                rel.iaddP(MMod(Mon(), (uint32_t)i), tmp1);
                rels_graded[d].push_back(rel);
            }
        }
    }

    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg_max && (d <= deg_max_rels || !criticals_.empty()); ++d) {
        // std::cout << "d=" << d << '\n';
        CriPair1d pairs_d = Criticals(d);
        auto& rels_d = rels_graded[d];
        Mod1d rels_tmp(pairs_d.size() + rels_d.size());

        for (size_t i = 0; i < pairs_d.size(); ++i)
            pairs_d[i].SijMod(*pGb_, *this, rels_tmp[i], tmp1, tmp2);
        size_t offset = pairs_d.size();
        ut::for_each_seq(offset, [this, &pairs_d, &rels_tmp](size_t i) { rels_tmp[i] = Reduce(std::move(rels_tmp[i])); });
        ut::for_each_seq(rels_d.size(), [this, &rels_d, &rels_tmp, offset](size_t i) { rels_tmp[offset + i] = Reduce(rels_d[i]); });

        /* Triangulate these relations */
        Mod1d gb_rels_d;
        for (auto& rel : rels_tmp) {
            for (const Mod& rel1 : gb_rels_d)
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

}  // namespace alg2