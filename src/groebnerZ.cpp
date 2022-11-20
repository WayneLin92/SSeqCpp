#include "groebnerZ.h"
#include "myio.h"
#include "utility.h"
#include <set>

namespace algZ {

void CriPair::SetFromLM(CriPair& result, const Mon& lead1, const Mon& lead2, int O1, int O2, int i, int j, const AdamsDeg1d& gen_degs)
{
    alg2::Mon m1m0, m1m1, m2m0, m2m1;
    alg2::detail::MutualQuotient(m1m0, m2m0, lead1.m0(), lead2.m0());
    alg2::detail::MutualQuotient(m1m1, m2m1, lead1.m1(), lead2.m1());
    int c_min = std::min(lead1.c(), lead2.c());
    result.m1 = Mon(lead2.c() - c_min, m1m0, m1m1, -1024);
    result.m1.SetFil(gen_degs);
    result.m2 = Mon(lead1.c() - c_min, m2m0, m2m1, -1024);
    result.m2.SetFil(gen_degs);
    result.i1 = i;
    result.i2 = j;
    result.O = std::min({result.m1.fil() + O1, result.m2.fil() + O2, FIL_MAX + 1});
    result.trace_m2 = result.m2.Trace();
}

void CriPair::SijP(const Groebner& gb, Poly& result, Poly& tmp) const
{
    result.data.clear();
    auto& gen_2tor_degs = gb.gen_2tor_degs();
    result.iaddmulP(m1, gb[i1], tmp, gen_2tor_degs);  //// TODO: improve this
    result.iaddmulP(m2, gb[i2], tmp, gen_2tor_degs);
}

void CriPair::SijMod(const Groebner& gb, const GroebnerMod& gbm, Mod& result, Mod& tmp) const
{
    result.data.clear();
    auto& gen_2tor_degs = gb.gen_2tor_degs();
    if (i2 >= gbm.size())
        throw MyException(0, "BUG");
    if (i1 & FLAG_INDEX_X) {
        int v = gbm[i2].GetLead().v;
        result.iaddmulP(m1, Mod(gb[i1 ^ FLAG_INDEX_X], v, gbm.v_degs()[v].s), tmp, gen_2tor_degs);
    }
    else
        result.iaddmulP(m1, gbm[i1], tmp, gen_2tor_degs);  //// TODO: improve this
    result.isubmulP(m2, gbm[i2], tmp, gen_2tor_degs);
}

void GbCriPairs::Minimize(const Mon1d& leads, const int1d& leads_O, int d, const AdamsDeg1d& gen_degs_)
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
            gcd.SetFil(gen_degs_);
            Mon m2 = div_unsigned(leads[i], gcd);
            int O = std::min({leads[j].fil() - gcd.fil() + leads_O[i], m2.fil() + leads_O[j], FIL_MAX + 1});
            if (b_min_pairs_d.find((int)j) != b_min_pairs_d.end()) {
                auto& pairs = b_min_pairs_d.at((int)j);
                auto p = std::find_if(pairs.begin(), pairs.end(), [&m2, O](const CriPair& c) { return c.m2 == m2 && (O > FIL_MAX || c.O <= O); });
                /* Remove it from `buffer_min_pairs_` */
                if (p != pairs.end()) {
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
                    if (c->O > FIL_MAX || c->O - c->m2.fil() >= O - m2.fil()) {
                        Mon m1 = div_unsigned(leads[j], gcd);
                        if (detail::HasGCD(c->m1, m1)) {
                            j = NULL_J;
                            break;
                        }
                        else {
                            j = c->i1;
                            if (i > j)
                                std::swap(i, j);
                            else if (i == j)
                                j = NULL_J;
                            break;
                        }
                    }
                }
            }
#ifndef NDEBUG
            if (c == end)
                throw MyException(0x3bc25c87U, "Should not happen because gb_ is groebner");
#endif
        }
    }

    /* Delete `bufbuffer_redundent_pairs_[d]` */
    buffer_redundent_pairs_.erase(d);
}

void GbCriPairs::Minimize(const Mon1d& leadsx, const int1d& leadsx_O, const MMod1d& leads, const int1d& leads_O, int d, const AdamsDeg1d& gen_degs_)
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
            Mon leadsi = (i & FLAG_INDEX_X) ? leadsx[i ^ FLAG_INDEX_X] : leads[i].m; /* Warning leads[i].m.fil contains v fil */
            Mon gcd = GCD(leadsi, leads[j].m);
            gcd.SetFil(gen_degs_);
            Mon m2 = div_unsigned(leadsi, gcd);
            m2.SetFil(gen_degs_);
            int Oi;
            if (i & FLAG_INDEX_X)
                Oi = leads[j].fil() - gcd.fil() + leadsx_O[i ^ FLAG_INDEX_X];
            else {
                int vj_fil = leads[j].fil() - GetDeg(leads[j].m, gen_degs_).s;
                Oi = leads[j].fil() - gcd.fil() + leads_O[i] - vj_fil;
            }
            int O = std::min({Oi, m2.fil() + leads_O[j], FIL_MAX + 1});
            if (b_min_pairs_d.find((int)j) != b_min_pairs_d.end()) {
                auto& pairs = b_min_pairs_d.at((int)j);
                auto p = std::find_if(pairs.begin(), pairs.end(), [&m2, O](const CriPair& c) { return c.m2 == m2 && (O > FIL_MAX || c.O <= O); });
                /* Remove it from `buffer_min_pairs_` */
                if (p != pairs.end()) {
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
                    if (c->O > FIL_MAX || c->O - c->m2.fil() >= O - m2.fil()) {
                        Mon m1 = div_unsigned(leads[j].m, gcd);
                        if (detail::HasGCD(c->m1, m1)) {
                            j = NULL_J;
                            break;
                        }
                        else {
                            j = c->i1;
                            if (i > j || (j & FLAG_INDEX_X))
                                std::swap(i, j);
                            if (i == j || j & FLAG_INDEX_X)
                                j = NULL_J;
                            break;
                        }
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

bool GEPair(const CriPair& p1, const CriPair& p2)
{
    return p1.O > FIL_MAX || p1.O - p1.m2.fil() >= p2.O - p2.m2.fil();
}

bool NEPair(const CriPair& p1, const CriPair& p2)
{
    return p1.O > FIL_MAX || p2.O > FIL_MAX || p1.O - p1.m2.fil() != p2.O - p2.m2.fil();
}

void GbCriPairs::AddToBuffers(const Mon1d& leads, const MonTrace1d& traces, const int1d& leads_O, const Mon& mon, int O, const AdamsDeg1d& gen_degs)
{
    int1d gen_degs_s(gen_degs.size()), gen_degs_t(gen_degs.size());
    for (size_t i = 0; i < gen_degs.size(); ++i) {
        gen_degs_s[i] = gen_degs[i].s;
        gen_degs_t[i] = gen_degs[i].t;
    }
    size_t lead_size = leads.size();
    if (gb_.size() < lead_size + 1)
        gb_.resize(lead_size + 1);
    new_pairs__.resize(leads.size());
    MonTrace t = mon.Trace();

    /* Populate `new_pairs__` */
    for (size_t i = 0; i < lead_size; ++i) {
        new_pairs__[i].first = -1;
        if (detail::HasGCD(leads[i], mon, traces[i], t)) {
            int d_pair = detail::DegLCM(leads[i], mon, gen_degs_t);
            if (d_pair <= deg_trunc_) {
                new_pairs__[i].first = d_pair;
                CriPair::SetFromLM(new_pairs__[i].second, leads[i], mon, leads_O[i], O, (uint32_t)i, (uint32_t)lead_size, gen_degs);
            }
        }
    }

    /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
    for (size_t j = 1; j < new_pairs__.size(); ++j) {
        auto& npj = new_pairs__[j];
        if (npj.first != -1) {
            for (size_t i = 0; i < j; ++i) {
                auto& npi = new_pairs__[i];
                if (npi.first != -1) {
                    if (divisible(npi.second.m2, npj.second.m2, npi.second.trace_m2, npj.second.trace_m2) && GEPair(npi.second, npj.second)) {
                        npj.first = -1;
                        break;
                    }
                    else if (divisible(npj.second.m2, npi.second.m2, npj.second.trace_m2, npi.second.trace_m2) && GEPair(npj.second, npi.second)) {
                        npi.first = -1;
                    }
                    else if (!detail::HasGCD(npi.second.m1, npj.second.m1) && detail::HasGCD(leads[i], leads[j], traces[i], traces[j]) && NEPair(npi.second, npj.second)) {
                        int dij = detail::DegLCM(leads[i], leads[j], gen_degs_t);
                        if (dij <= deg_trunc_) {
                            buffer_redundent_pairs_[dij].insert(ut::Bind(i, j));
                        }
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < new_pairs__.size(); ++i) {
        if (new_pairs__[i].first != -1) {
            gb_[lead_size].push_back(new_pairs__[i].second);
            buffer_min_pairs_[new_pairs__[i].first][(int)lead_size].push_back(std::move(new_pairs__[i].second));
        }
    }
}

void GbCriPairs::AddToBuffers(const Mon1d& leadsx, const MonTrace1d& tracesx, const int1d& leadsx_O, const MMod1d& leads, const MonTrace1d& traces, const int1d& leads_O, const MMod& mon, int O, const AdamsDeg1d& gen_degs, const AdamsDeg1d& v_degs)
{
    int1d gen_degs_s(gen_degs.size()), gen_degs_t(gen_degs.size());
    for (size_t i = 0; i < gen_degs.size(); ++i) {
        gen_degs_s[i] = gen_degs[i].s;
        gen_degs_t[i] = gen_degs[i].t;
    }

    size_t leadsx_size = leadsx.size();
    size_t leads_size = leads.size();
    if (gb_.size() < leads_size + 1)
        gb_.resize(leads_size + 1);
    new_pairs__.resize(leadsx_size + leads_size);
    MonTrace t = mon.m.Trace();

    /* Populate `new_pairs__` */
    for (size_t i = 0; i < leadsx_size; ++i) {
        new_pairs__[i].first = -1;
        int d_pair = detail::DegLCM(leadsx[i], mon.m, gen_degs_t) + v_degs[mon.v].t;
        if (d_pair <= deg_trunc_) {
            new_pairs__[i].first = d_pair;
            CriPair::SetFromLM(new_pairs__[i].second, leadsx[i], mon.m, leadsx_O[i] + v_degs[mon.v].s, O, uint32_t(i | FLAG_INDEX_X), (uint32_t)leads_size, gen_degs);
        }
    }
    for (size_t i = 0; i < leads_size; ++i) {
        auto& new_pairs_i = new_pairs__[leadsx_size + i];
        new_pairs_i.first = -1;
        if (leads[i].v == mon.v) {
            int d_pair = detail::DegLCM(leads[i].m, mon.m, gen_degs_t) + v_degs[mon.v].t;
            if (d_pair <= deg_trunc_) {
                new_pairs_i.first = d_pair;
                CriPair::SetFromLM(new_pairs_i.second, leads[i].m, mon.m, leads_O[i], O, (uint32_t)i, (uint32_t)leads_size, gen_degs);
            }
        }
    }

    /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
    for (size_t j = 0; j < new_pairs__.size(); ++j) {
        auto& npj = new_pairs__[j];
        if (npj.first != -1) {
            for (size_t i = 0; i < j; ++i) {
                auto& npi = new_pairs__[i];
                if (npi.first != -1) {
                    if (divisible(npi.second.m2, npj.second.m2, npi.second.trace_m2, npj.second.trace_m2) && GEPair(npi.second, npj.second)) {
                        npj.first = -1;
                        break;
                    }
                    else if (divisible(npj.second.m2, npi.second.m2, npj.second.trace_m2, npi.second.trace_m2) && GEPair(npj.second, npi.second)) {
                        npi.first = -1;
                    }
                    else if (!detail::HasGCD(npi.second.m1, npj.second.m1) && j >= leadsx_size && NEPair(npi.second, npj.second)) {
                        auto& leadsi = i < leadsx_size ? leadsx[i] : leads[i - leadsx_size].m;
                        const size_t j1 = j - leadsx_size;
                        int dij = detail::DegLCM(leadsi, leads[j1].m, gen_degs_t) + v_degs[leads[j1].v].t;
                        if (dij <= deg_trunc_) {
                            buffer_redundent_pairs_[dij].insert(ut::Bind(i < leadsx_size ? (i | FLAG_INDEX_X) : i - leadsx_size, j1));
                        }
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < new_pairs__.size(); ++i) {
        if (new_pairs__[i].first != -1) {
            gb_[leads_size].push_back(new_pairs__[i].second);
            buffer_min_pairs_[new_pairs__[i].first][(int)leads_size].push_back(std::move(new_pairs__[i].second));
        }
    }
}

void GbCriPairs::AddToBuffersX(const Mon1d& leadsx, const MonTrace1d& tracesx, const int1d& leadsx_O, const MMod1d& leads, const MonTrace1d& traces, const int1d& leads_O, const AdamsDeg1d& gen_degs, const AdamsDeg1d& v_degs,
                               size_t i_start)  // TODO: test
{
    for (size_t l = 0; l < leads.size(); ++l) {
        auto& mon = leads[l];
        auto& t = traces[l];

        int1d gen_degs_s(gen_degs.size()), gen_degs_t(gen_degs.size());
        for (size_t i = 0; i < gen_degs.size(); ++i) {
            gen_degs_s[i] = gen_degs[i].s;
            gen_degs_t[i] = gen_degs[i].t;
        }

        size_t leadsx_size = leadsx.size();
        new_pairs__.resize(leadsx_size - i_start);

        /* Populate `new_pairs__` */
        for (size_t i = i_start; i < leadsx_size; ++i) {
            auto& new_pairs_i = new_pairs__[i - i_start];
            new_pairs_i.first = -1;
            int d_pair = detail::DegLCM(leadsx[i], mon.m, gen_degs_t) + v_degs[mon.v].t;
            if (d_pair <= deg_trunc_) {
                new_pairs_i.first = d_pair;
                CriPair::SetFromLM(new_pairs_i.second, leadsx[i], mon.m, leadsx_O[i] + v_degs[mon.v].s, leads_O[l], uint32_t(i | FLAG_INDEX_X), (uint32_t)l, gen_degs);
            }
        }

        /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
        for (size_t j = size_t(i_start + 1); j < leadsx_size; ++j) {
            auto& new_pairs_j = new_pairs__[j - i_start];
            if (new_pairs_j.first != -1) {
                for (size_t i = i_start; i < j; ++i) {
                    auto& new_pairs_i = new_pairs__[i - i_start];
                    if (new_pairs_i.first != -1) {
                        if (divisible(new_pairs_i.second.m2, new_pairs_j.second.m2, new_pairs_i.second.trace_m2, new_pairs_j.second.trace_m2) && GEPair(new_pairs_i.second, new_pairs_j.second)) {
                            new_pairs_j.first = -1;
                            break;
                        }
                        else if (divisible(new_pairs_j.second.m2, new_pairs_i.second.m2, new_pairs_j.second.trace_m2, new_pairs_i.second.trace_m2) && GEPair(new_pairs_j.second, new_pairs_i.second))
                            new_pairs_i.first = -1;
                    }
                }
            }
        }
        for (size_t i = i_start; i < leadsx_size; ++i) {
            auto& new_pairs_i = new_pairs__[i - i_start];
            if (new_pairs_i.first != -1) {
                gb_[l].push_back(new_pairs_i.second);
                buffer_min_pairs_[new_pairs_i.first][(int)l].push_back(std::move(new_pairs_i.second));
            }
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

void GbCriPairs::init(const Mon1d& leads, const MonTrace1d& traces, const int1d& leads_O, const AdamsDeg1d& gen_degs, int t_min_buffer)
{
    Mon1d tmp_leads;
    for (size_t i = 0; i < leads.size(); ++i) {
        AddToBuffers(tmp_leads, traces, leads_O, leads[i], leads_O[i], gen_degs);
        tmp_leads.push_back(leads[i]);
    }
    RemoveSmallKey(buffer_min_pairs_, t_min_buffer);
    RemoveSmallKey(buffer_redundent_pairs_, t_min_buffer);
}

void GbCriPairs::init(const Mon1d& leadsx, const MonTrace1d& tracesx, const int1d& leadsx_O, const MMod1d& leads, const MonTrace1d& traces, const int1d& leads_O, const AdamsDeg1d& gen_degs, const AdamsDeg1d& v_degs, int t_min_buffer)
{
    MMod1d tmp_leads;
    for (size_t i = 0; i < leads.size(); ++i) {
        AddToBuffers(leadsx, tracesx, leadsx_O, tmp_leads, traces, leads_O, leads[i], leads_O[i], gen_degs, v_degs);
        tmp_leads.push_back(leads[i]);
    }
    RemoveSmallKey(buffer_min_pairs_, t_min_buffer);
    RemoveSmallKey(buffer_redundent_pairs_, t_min_buffer);
}

Groebner::Groebner(int deg_trunc, AdamsDeg1d gen_degs, Poly1d polys, bool bDynamic) : criticals_(deg_trunc), gen_degs_(std::move(gen_degs)), gen_2tor_degs_(gen_degs_.size(), FIL_MAX + 1)
{
    for (auto& p : polys) {
        AdamsDeg deg = GetDeg(p.GetLead(), gen_degs_);
        push_back_data_init(std::move(p), deg);
    }
    if (bDynamic)
        criticals_.init(leads_, traces_, leads_O_, gen_degs_, deg_trunc + 1);
}

void pop_indices_by_ub(int1d& indices, int ub)
{
    auto p = std::lower_bound(indices.begin(), indices.end(), ub);
    indices.erase(p, indices.end());
}

void Groebner::Pop(size_t gen_size, size_t rel_size)
{
    gen_degs_.resize(gen_size);
    gen_2tor_degs_.resize(gen_size);

    criticals_.Pop(rel_size);
    data_.resize(rel_size);
    leads_.resize(rel_size);
    traces_.resize(rel_size);
    leads_O_.resize(rel_size);
    for (auto& [_, indices] : leads_group_by_key_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& [_, indices] : leads_group_by_deg_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& indices : leads_group_by_last_gen_)
        pop_indices_by_ub(indices, (int)rel_size);
}

void Groebner::debug_print() const
{
    std::cout << "gen_degs_.size()=" << gen_degs_.size() << '\n';
    std::cout << "data_.size()=" << data_.size() << '\n';
    int sum = 0;
    for (auto& [_, indices] : leads_group_by_key_)
        sum += indices.size();
    std::cout << "leads_group_by_key_=" << sum << '\n';

    sum = 0;
    for (auto& [_, indices] : leads_group_by_deg_)
        sum += indices.size();
    std::cout << "leads_group_by_deg_=" << sum << '\n';

    sum = 0;
    for (auto& indices : leads_group_by_last_gen_)
        sum += indices.size();
    std::cout << "leads_group_by_last_gen_=" << sum << '\n';

    sum = 0;
    for (auto& cps : criticals_.gb_)
        sum += cps.size();
    std::cout << "criticals_=" << sum << '\n';
}

int Groebner::IndexOfDivisibleLeading(const Mon& mon, int eff_min) const
{
    auto t = mon.Trace();
    for (int i = -1; i < (int)mon.m0().size(); ++i) {
        for (int j = -1; j < (int)mon.m1().size(); ++j) {
            auto key = TypeIndexKey{(i == -1 ? 0 : mon.m0()[i].g() + 1) + (j == -1 ? 0 : ((mon.m1()[j].g() + 1) << 16))};
            auto p = leads_group_by_key_.find(key);
            if (p != leads_group_by_key_.end()) {
                for (int k : p->second) {
                    if (divisible(leads_[k], mon, traces_[k], t) && data_[k].EffNum() >= eff_min)
                        return k;
                }
            }
        }
    }
    return -1;
}

int Groebner::IndexOfDivisibleLeadingV2(const Mon& mon) const
{
    auto t = mon.Trace();
    int result = -1, eff = 0;
    for (int i = -1; i < (int)mon.m0().size(); ++i) {
        for (int j = -1; j < (int)mon.m1().size(); ++j) {
            auto key = TypeIndexKey{(i == -1 ? 0 : mon.m0()[i].g() + 1) + (j == -1 ? 0 : ((mon.m1()[j].g() + 1) << 16))};
            auto p = leads_group_by_key_.find(key);
            if (p != leads_group_by_key_.end()) {
                for (int k : p->second) {
                    if (divisible(leads_[k], mon, traces_[k], t)) {
                        if (data_[k].EffNum() > eff) {
                            result = k;
                            eff = data_[k].EffNum();
                        }
                    }
                }
            }
        }
    }
    return result;
}

Poly Groebner::Reduce(Poly poly) const
{
    Poly tmp_prod, tmp;
    size_t index = 0;
    while (index < poly.data.size() && !poly.data[index].IsUnKnown()) {
        if (poly.data[index].IsTrivial(gen_2tor_degs_)) {
            poly.data.erase(poly.data.begin() + index);
            continue;
        }
        int eff_min = poly.UnknownFil() > FIL_MAX ? FIL_MAX + 1 : poly.UnknownFil() - poly.data[index].fil();
        int gb_index = IndexOfDivisibleLeading(poly.data[index], eff_min);
        if (gb_index != -1) {
            Mon q = div_unsigned(poly.data[index], data_[gb_index].GetLead());
            int sign = get_sign(q, data_[gb_index].GetLead(), poly.data[index]);
            if (sign == 0)
                poly.data.erase(poly.data.begin() + index); /* Remove the monomial 2x^2 which is zero*/
            else if (sign == 1)
                poly.isubmulP(q, data_[gb_index], tmp, gen_2tor_degs_);
            else
                poly.iaddmulP(q, data_[gb_index], tmp, gen_2tor_degs_);
        }
        else {
            if (index > 0 && MultipleOf(poly.data[index], poly.data[0])) {
                int c = poly.data[index].c() - poly.data[0].c();
                Poly poly1 = poly;
                poly.isubmulP(Mon::twoTo(c), poly1, tmp, gen_2tor_degs_);
            }
            else
                ++index;
        }
    }
    return poly;
}

Poly Groebner::ReduceV2(Poly poly) const
{
    Poly tmp_prod, tmp;
    size_t index = 0;
    while (index < poly.data.size() && !poly.data[index].IsUnKnown()) {
        if (poly.data[index].IsTrivial(gen_2tor_degs_)) {
            poly.data.erase(poly.data.begin() + index);
            continue;
        }
        int gb_index = IndexOfDivisibleLeadingV2(poly.data[index]);
        if (gb_index != -1) {
            Mon q = div_unsigned(poly.data[index], data_[gb_index].GetLead());
            int sign = get_sign(q, data_[gb_index].GetLead(), poly.data[index]);
            if (sign == 0)
                poly.data.erase(poly.data.begin()); /* Remove the leading monomial which is zero*/
            else if (sign == 1)
                poly.isubmulP(q, data_[gb_index], tmp, gen_2tor_degs_);
            else
                poly.iaddmulP(q, data_[gb_index], tmp, gen_2tor_degs_);
        }
        else
            ++index;
    }
    return poly;
}

void Groebner::AddRels(Poly1d rels, int deg_max)
{
    int d_trunc = deg_trunc();
    if (deg_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");
    Poly tmp;

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, Poly1d> rels_graded;
    for (auto& rel : rels) {
        if (rel) {
            int d = GetDegT(rel.GetLead(), gen_degs_);
            if (d <= deg_max)
                rels_graded[d].push_back(std::move(rel));
        }
    }

    for (int d = 1; d <= deg_max && ((!rels_graded.empty() && d <= rels_graded.rbegin()->first) || !criticals_.empty()); ++d) {
        int nextd = criticals_.NextD();
        if (nextd != -1 && nextd < d) {
            if (nextd < d - 1)
                throw MyException(0xc3b447bcU, "buffer_min_pairs_ contains degree < d - 1");
            --d;
        }

        CriPair1d pairs_d = Criticals(d);
        size_t pairs_d_size = pairs_d.size();
        auto& rels_d = rels_graded[d];
        Poly1d rels_tmp(pairs_d_size + rels_d.size());
        for (size_t i = 0; i < pairs_d_size; ++i)
            pairs_d[i].SijP(*this, rels_tmp[i], tmp);
        for (size_t i = 0; i < rels_d.size(); ++i)
            rels_tmp[pairs_d_size + i] = std::move(rels_d[i]);
        std::sort(rels_tmp.begin(), rels_tmp.end(), [](const Poly& p1, const Poly& p2) { return p1.EffNum() > p2.EffNum(); });

        for (auto& p : rels_tmp) {
            p = Reduce(std::move(p));
            if (p && !p.GetLead().IsUnKnown()) {
                AdamsDeg deg = GetDeg(p.GetLead(), gen_degs_);
                if (deg.t <= d_trunc) {
                    if (deg_max == d_trunc && deg.t > d)
                        rels_graded[deg.t].push_back(std::move(p));
                    else
                        push_back_data(std::move(p), deg);
                }
            }
        }
    }
    criticals_.ClearBuffer();
}

void Groebner::SimplifyRels()
{
    Poly1d data1 = std::move(data_);
    ResetRels();
    AddRels(std::move(data1), criticals_.deg_trunc());
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
        if (n)
            power = gb.Reduce(power * power);
    }
}

Mon1d Groebner::GenBasis(AdamsDeg deg, const std::map<AdamsDeg, Mon1d>& basis) const
{
    Mon1d result;
    for (size_t gen_id = 0; gen_id < gen_degs_.size(); ++gen_id) {
        AdamsDeg d1 = deg - gen_degs_[gen_id];
        auto p = basis.find(d1);
        if (p != basis.end()) {
            for (auto& m : p->second) {
                if (!m || (int)gen_id >= m.backg()) {
                    Mon mon = mul_unsigned(m, Gen((uint32_t)gen_id));
                    if (gen_id >= leads_group_by_last_gen_.size() || std::none_of(leads_group_by_last_gen_[gen_id].begin(), leads_group_by_last_gen_[gen_id].end(), [this, &mon](int i) { return divisible(this->leads()[i], mon); }))
                        result.push_back(std::move(mon));
                }
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

Poly1d Groebner::RelsLF(AdamsDeg deg) const
{
    Poly1d result;
    auto p = leads_group_by_deg_.find(deg);
    if (p != leads_group_by_deg_.end())
        for (int i : p->second)
            result.push_back(data_[i].LF());
    return result;
}

/********************************* Modules ****************************************/

GroebnerMod::GroebnerMod(Groebner* pGb, int deg_trunc, AdamsDeg1d v_degs, Mod1d polys, bool bDynamic) : pGb_(pGb), criticals_(deg_trunc), v_degs_(std::move(v_degs)), old_pGb_size_(pGb->leads_.size())
{
    for (auto& p : polys) {
        AdamsDeg deg = GetDeg(p.GetLead(), pGb_->gen_degs(), v_degs_);
        push_back_data_init(std::move(p), deg);
    }
    if (bDynamic)
        criticals_.init(pGb_->leads_, pGb_->traces_, pGb_->leads_O_, leads_, traces_, leads_O_, pGb_->gen_degs(), v_degs_, deg_trunc + 1);
}

MMod1d GroebnerMod::GenBasis(AdamsDeg deg, const std::map<AdamsDeg, Mon1d>& basis) const
{
    MMod1d result;
    for (size_t v = 0; v < v_degs_.size(); ++v) {
        AdamsDeg d1 = deg - v_degs_[v];
        auto p = basis.find(d1);
        if (p != basis.end()) {
            for (auto& m : p->second) {
                if (v >= leads_group_by_v_.size() || std::none_of(leads_group_by_v_[v].begin(), leads_group_by_v_[v].end(), [this, &m](int i) { return divisible(this->leads()[i].m, m); }))
                    result.emplace_back(m, (int)v, v_degs_[v].s);
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

Mod1d GroebnerMod::RelsLF(AdamsDeg deg) const
{
    Mod1d result;
    auto p = leads_group_by_deg_.find(deg);
    if (p != leads_group_by_deg_.end())
        for (int i : p->second)
            result.push_back(data_[i].LF());
    return result;
}

void GroebnerMod::Pop(size_t gen_size, size_t rel_size)
{
    old_pGb_size_ = pGb_->leads_.size();
    v_degs_.resize(gen_size);

    criticals_.Pop(old_pGb_size_, rel_size);
    data_.resize(rel_size);
    leads_.resize(rel_size);
    traces_.resize(rel_size);
    leads_O_.resize(rel_size);
    for (auto& [_, indices] : leads_group_by_key_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& [_, indices] : leads_group_by_deg_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& indices : leads_group_by_v_)
        pop_indices_by_ub(indices, (int)rel_size);
}

void GroebnerMod::debug_print() const
{
    std::cout << "v_degs_.size()=" << v_degs_.size() << '\n';
    std::cout << "data_.size()=" << data_.size() << '\n';
    int sum = 0;
    for (auto& [_, indices] : leads_group_by_key_)
        sum += indices.size();
    std::cout << "leads_group_by_key_=" << sum << '\n';

    sum = 0;
    for (auto& [_, indices] : leads_group_by_deg_)
        sum += indices.size();
    std::cout << "leads_group_by_deg_=" << sum << '\n';

    sum = 0;
    for (auto& indices : leads_group_by_v_)
        sum += indices.size();
    std::cout << "leads_group_by_v_=" << sum << '\n';

    sum = 0;
    for (auto& cps : criticals_.gb_)
        sum += cps.size();
    std::cout << "criticals_=" << sum << '\n';
}

int GroebnerMod::IndexOfDivisibleLeading(const MMod& mon, int eff_min) const
{
    auto t = mon.m.Trace();
    int i_max = int(mon.m.m0().size() + mon.m.m1().size());
    for (int i = -2; i < i_max; ++i) {
        uint32_t backg = 0;
        if (i == -1)
            backg = 1;
        else if (i >= 0)
            backg = i < (int)mon.m.m0().size() ? mon.m.m0()[i].g() + 1 : mon.m.m1()[(size_t)i - mon.m.m0().size()].g() + 1;
        auto key = TypeIndexKey{mon.v + (backg << 16)};
        auto p = leads_group_by_key_.find(key);
        if (p != leads_group_by_key_.end()) {
            for (int k : p->second) {
                if (divisible(leads_[k].m, mon.m, traces_[k], t) && data_[k].EffNum() >= eff_min)
                    return k;
            }
        }
    }
    return -1;
}

int GroebnerMod::IndexOfDivisibleLeadingV2(const MMod& mon) const
{
    int result = -1, eff = 0;
    auto t = mon.m.Trace();
    int i_max = int(mon.m.m0().size() + mon.m.m1().size());
    for (int i = -2; i < i_max; ++i) {
        uint32_t backg = 0;
        if (i == -1)
            backg = 1;
        else if (i >= 0)
            backg = i < (int)mon.m.m0().size() ? mon.m.m0()[i].g() + 1 : mon.m.m1()[(size_t)i - mon.m.m0().size()].g() + 1;
        auto key = TypeIndexKey{mon.v + (backg << 16)};
        auto p = leads_group_by_key_.find(key);
        if (p != leads_group_by_key_.end()) {
            for (int k : p->second) {
                if (divisible(leads_[k].m, mon.m, traces_[k], t)) {
                    if (data_[k].EffNum() > eff) {
                        result = k;
                        eff = data_[k].EffNum();
                    }
                }
            }
        }
    }
    return result;
}

Mod GroebnerMod::Reduce(Mod x) const
{
    Mod tmp_prodm, tmpm1, tmpm2;
    Poly tmp_prod, tmp1, tmp2;
    size_t index = 0;
    while (index < x.data.size() && !x.data[index].IsUnKnown()) {
        if (x.data[index].m.IsTrivial(pGb_->gen_2tor_degs())) {
            x.data.erase(x.data.begin() + index);
            continue;
        }
        int eff_min = x.UnknownFil() > FIL_MAX ? FIL_MAX + 1 : x.UnknownFil() - x.data[index].fil();
        int gbmod_index = IndexOfDivisibleLeading(x.data[index], eff_min);
        if (gbmod_index != -1) {
            Mon q = div_unsigned(x.data[index].m, data_[gbmod_index].GetLead().m);
            int sign = get_sign(q, data_[gbmod_index].GetLead().m, x.data[index].m);
            if (sign == 1)
                x.isubmulP(q, data_[gbmod_index], tmpm1, pGb_->gen_2tor_degs());
            else if (sign == -1)
                x.iaddmulP(q, data_[gbmod_index], tmpm1, pGb_->gen_2tor_degs());
            else
                throw MyException(0x8dd4cebcU, "Something is wrong.");
        }
        else {
            int gb_index = pGb_->IndexOfDivisibleLeading(x.data[index].m, eff_min);
            if (gb_index != -1) {
                Mon q = div_unsigned(x.data[index].m, pGb_->data()[gb_index].GetLead()); /* q has extra filtration from v */
                int sign = get_sign(q, pGb_->data()[gb_index].GetLead(), x.data[index].m);
                int v = x.data[index].v;
                if (sign == 1)
                    x.isubmulP(q, Mod(pGb_->data()[gb_index], v, 0), tmpm1, pGb_->gen_2tor_degs());
                else if (sign == -1)
                    x.iaddmulP(q, Mod(pGb_->data()[gb_index], v, 0), tmpm1, pGb_->gen_2tor_degs());
                else
                    throw MyException(0x3e2f96c0U, "Something is wrong.");
            }
            else {
                if (index > 0 && MultipleOf(x.data[index], x.data[0])) {
                    int c = x.data[index].c() - x.data[0].c();
                    Mod x1 = x;
                    x.isubmulP(Mon::twoTo(c), x1, tmpm1, pGb_->gen_2tor_degs());
                }
                else
                    ++index;
            }
        }
    }
    return x;
}

Mod GroebnerMod::ReduceV2(Mod x) const
{
    Mod tmp_prodm, tmpm1, tmpm2;
    Poly tmp_prod, tmp1, tmp2;
    size_t index = 0;
    while (index < x.data.size() && !x.data[index].IsUnKnown()) {
        if (x.data[index].m.IsTrivial(pGb_->gen_2tor_degs())) {
            x.data.erase(x.data.begin() + index);
            continue;
        }
        int gbmod_index = IndexOfDivisibleLeadingV2(x.data[index]);
        int gb_index = pGb_->IndexOfDivisibleLeadingV2(x.data[index].m);
        if (gbmod_index != -1 && (gb_index == -1 || data_[gbmod_index].EffNum() >= pGb_->data()[gb_index].EffNum())) {
            Mon q = div_unsigned(x.data[index].m, data_[gbmod_index].GetLead().m);
            int sign = get_sign(q, data_[gbmod_index].GetLead().m, x.data[index].m);
            if (sign == 1)
                x.isubmulP(q, data_[gbmod_index], tmpm1, pGb_->gen_2tor_degs());
            else if (sign == -1)
                x.iaddmulP(q, data_[gbmod_index], tmpm1, pGb_->gen_2tor_degs());
            else
                throw MyException(0x8dd4cebcU, "Something is wrong.");
        }
        else if (gb_index != -1) {
            Mon q = div_unsigned(x.data[index].m, pGb_->data()[gb_index].GetLead()); /* q has extra filtration from v */
            int sign = get_sign(q, pGb_->data()[gb_index].GetLead(), x.data[index].m);
            int v = x.data[index].v;
            if (sign == 1)
                x.isubmulP(q, Mod(pGb_->data()[gb_index], v, 0), tmpm1, pGb_->gen_2tor_degs());  //// TODO: improve
            else if (sign == -1)
                x.iaddmulP(q, Mod(pGb_->data()[gb_index], v, 0), tmpm1, pGb_->gen_2tor_degs());
            else
                throw MyException(0x3e2f96c0U, "Something is wrong.");
        }
        else
            ++index;
    }
    return x;
}

void GroebnerMod::AddRels(Mod1d rels, int deg_max)
{
    int d_trunc = deg_trunc();
    if (deg_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");
    if (old_pGb_size_ != pGb_->leads_.size()) {
        criticals_.AddToBuffersX(pGb_->leads_, pGb_->traces_, pGb_->leads_O_, leads_, traces_, leads_O_, pGb_->gen_degs(), v_degs_, old_pGb_size_);
        old_pGb_size_ = pGb_->leads_.size();
    }

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, Mod1d> rels_graded;
    for (auto& rel : rels) {
        if (rel) {
            int d = GetDeg(rel.GetLead(), pGb_->gen_degs(), v_degs_).t;
            if (d <= deg_max)
                rels_graded[d].push_back(std::move(rel));
        }
    }

    Mod tmp;
    for (int d = 0; d <= deg_max && ((!rels_graded.empty() && d <= rels_graded.rbegin()->first) || !criticals_.empty()); ++d) {
        int nextd = criticals_.NextD();
        if (nextd != -1 && nextd < d) {
            if (nextd < d - 1)
                throw MyException(0xc3b447bcU, "buffer_min_pairs_ contains degree < d - 1");
            --d;
        }

        CriPair1d pairs_d = Criticals(d);
        size_t pairs_d_size = pairs_d.size();
        auto& rels_d = rels_graded[d];
        Mod1d rels_tmp(pairs_d_size + rels_d.size());
        for (size_t i = 0; i < pairs_d.size(); ++i)
            pairs_d[i].SijMod(*pGb_, *this, rels_tmp[i], tmp);
        for (size_t i = 0; i < rels_d.size(); ++i)
            rels_tmp[pairs_d_size + i] = std::move(rels_d[i]);
        std::sort(rels_tmp.begin(), rels_tmp.end(), [](const Mod& p1, const Mod& p2) { return p1.EffNum() > p2.EffNum(); });

        for (auto& p : rels_tmp) {
            p = Reduce(std::move(p));
            if (p && !p.GetLead().IsUnKnown()) {
                AdamsDeg deg = GetDeg(p.GetLead(), pGb_->gen_degs(), v_degs_);
                if (deg.t <= d_trunc) {
                    if (deg_max == d_trunc && deg.t > d)
                        rels_graded[deg.t].push_back(std::move(p));
                    else
                        push_back_data(std::move(p), deg);
                }
            }
        }
    }
    criticals_.ClearBuffer();
}

void GroebnerMod::SimplifyRels()
{
    Mod1d data1 = std::move(data_);
    ResetRels();
    AddRels(std::move(data1), criticals_.deg_trunc());
}

// void GroebnerMod::ToSubMod(const Mod1d& rels, int deg_max, int1d& index_ind)
//{
//     int d_trunc = deg_trunc();
//     if (deg_max > d_trunc)
//         throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");
//     Mod tmp;
//
//     const size_t rels_size = rels.size();
//     for (Mod& x : data_)
//         for (MMod m : x.data)
//             m.v += (uint32_t)rels_size;
//     v_degs_.resize(v_degs_.size() + rels_size);
//     for (size_t i = v_degs_.size(); i-- > rels_size;)
//         v_degs_[i] = v_degs_[i - rels_size];
//
//     /* Calculate the degrees of `rels` and group them by degree */
//     std::map<int, Mod1d> rels_graded;
//     for (size_t i = 0; i < rels.size(); ++i) {
//         if (rels[i]) {
//             AdamsDeg d = GetDeg(rels[i].GetLead().m, pGb_->gen_degs()) + v_degs_[(size_t)rels[i].GetLead().v + rels_size];
//             if (!rels_graded.empty() && d.t < rels_graded.rbegin()->first)
//                 throw MyException(0x39b4cfd4U, "rels is not ordered.");
//             if (d.t <= deg_max) {
//                 Mod rel = rels[i];
//                 v_degs_[i] = d;
//                 for (MMod& m : rel.data)
//                     m.v += (uint32_t)rels_size;
//                 rel.iaddP(MMod(Mon(), (uint32_t)i, v_degs_[i].s), tmp);
//                 rels_graded[d.t].push_back(rel);
//             }
//         }
//     }
//
//     int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
//     for (int d = 1; d <= deg_max && (d <= deg_max_rels || !criticals_.empty()); ++d) {
//         // std::cout << "d=" << d << '\n';
//         CriPair1d pairs_d = Criticals(d);
//         auto& rels_d = rels_graded[d];
//         Mod1d rels_tmp(pairs_d.size() + rels_d.size());
//
//         for (size_t i = 0; i < pairs_d.size(); ++i)
//             pairs_d[i].SijMod(*pGb_, *this, rels_tmp[i], tmp);
//         size_t offset = pairs_d.size();
//         ut::for_each_seq(offset, [this, &pairs_d, &rels_tmp](size_t i) { rels_tmp[i] = Reduce(std::move(rels_tmp[i])); });
//         ut::for_each_seq(rels_d.size(), [this, &rels_d, &rels_tmp, offset](size_t i) { rels_tmp[offset + i] = Reduce(rels_d[i]); });
//
//         /* Triangulate these relations */
//         Mod1d gb_rels_d;
//         for (auto& rel : rels_tmp) {
//             for (const Mod& rel1 : gb_rels_d)
//                 if (std::binary_search(rel.data.begin(), rel.data.end(), rel1.GetLead()))
//                     rel.iaddP(rel1, tmp);
//             if (rel)
//                 gb_rels_d.push_back(std::move(rel));
//         }
//
//         /* Add these relations */
//         for (auto& rel : gb_rels_d)
//             push_back_data(std::move(rel));
//     }
// }

}  // namespace algZ
