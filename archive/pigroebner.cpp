#include "groebnerZV2.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
// #include "benchmark.h"////
#include <set>

namespace algZ {

int NextO(const ut::map_seq2d<int, 0>& possEinf, int t_max, int stem, int O_min)
{
    int s_max = (stem + 3) / 2;  ////
    for (int s = O_min; s <= s_max; ++s) {
        if (stem + s > t_max || possEinf(stem, s))
            return s;
    }
    return algZ::FIL_MAX + 1;
}

template <typename T>
bool ExtendO(const ut::map_seq2d<int, 0>& possEinf, int t_max, int stem, T& rel)
{
    auto& m = rel.data.back();
    if (m.IsUnKnown()) {
        int s = algZ::NextO(possEinf, t_max, stem, m.fil());
        if (s > m.fil()) {
            if (s <= algZ::FIL_MAX)
                rel.data.back() = T::O(s);
            else
                rel.data.pop_back();
            return true;
        }
    }
    return false;
}

PiGroebner::PiGroebner(int deg_trunc, AdamsDeg1d gen_degs, Poly1d polys) : t_trunc_(deg_trunc), gen_degs_(std::move(gen_degs)), nodes_gen_2tor_degs_({{}})
{
    for (AdamsDeg deg : gen_degs_) {
        if (deg == AdamsDeg(1, 1))
            nodes_gen_2tor_degs_.back().push_back(FIL_MAX + 1);
        else
            nodes_gen_2tor_degs_.back().push_back((deg.stem() + 5) / 2);
    }
    for (auto& p : polys) {
        AdamsDeg deg = GetDeg(p.GetLead(), gen_degs_);
        push_back_data(std::move(p), deg);
    }
}

/* remove i from indices if i<=ub */
void pop_indices_by_ub(int1d& indices, int ub)
{
    auto p = std::lower_bound(indices.begin(), indices.end(), ub);
    indices.erase(p, indices.end());
}

void PiGroebner::Pop(size_t gen_size, size_t rel_size)
{
    gen_degs_.resize(gen_size);
    gen_Einf.resize(gen_size);
    nodes_gen_2tor_degs_.back().resize(gen_size);

    data_.resize(rel_size);
    leads_.resize(rel_size);
    traces_.resize(rel_size);
    leads_O_.resize(rel_size);
    for (auto& [_, indices] : leads_group_by_key_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& [_, indices] : leads_group_by_deg_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& [_, indices] : leads_group_by_t_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& indices : leads_group_by_last_gen_)
        pop_indices_by_ub(indices, (int)rel_size);
}

void PiGroebner::debug_print() const
{
    std::cout << "gen_degs_.size()=" << gen_degs_.size() << '\n';
    std::cout << "data_.size()=" << data_.size() << '\n';
    size_t sum = 0;
    for (auto& [_, indices] : leads_group_by_key_)
        sum += indices.size();
    std::cout << "leads_group_by_key_=" << sum << '\n';

    sum = 0;
    for (auto& [_, indices] : leads_group_by_deg_)
        sum += indices.size();
    std::cout << "leads_group_by_deg_=" << sum << '\n';

    sum = 0;
    for (auto& [_, indices] : leads_group_by_t_)
        sum += indices.size();
    std::cout << "leads_group_by_t_=" << sum << '\n';

    sum = 0;
    for (auto& indices : leads_group_by_last_gen_)
        sum += indices.size();
    std::cout << "leads_group_by_last_gen_=" << sum << '\n';
}

int PiGroebner::IndexOfDivisibleLeading(const Mon& mon, int eff_min) const
{
    auto t = mon.Trace();
    for (int i = mon.c() > 0 ? -1 : 0; i < (int)mon.m().size(); ++i) {
        for (int j = mon.c() > 0 ? -2 : -1; j < i; ++j) {
            uint32_t g1 = i == -1 ? 0 : mon.m()[i].g();
            uint32_t g2 = 0;
            if (j == -1 && mon.c() > 0)
                g2 = 1;
            else if (j >= 0)
                g2 = mon.m()[j].g() + 1;
            auto key = TypeIndexKey{g1 + (g2 << 16)};
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

int PiGroebner::IndexOfDivisibleLeadingV2(const Mon& mon) const
{
    auto t = mon.Trace();
    int result = -1, eff = 0;
    for (int i = mon.c() > 0 ? -1 : 0; i < (int)mon.m().size(); ++i) {
        for (int j = mon.c() > 0 ? -2 : -1; j < i; ++j) {
            uint32_t g1 = i == -1 ? 0 : mon.m()[i].g();
            uint32_t g2 = 0;
            if (j == -1 && mon.c() > 0)
                g2 = 1;
            else if (j >= 0)
                g2 = mon.m()[j].g() + 1;
            auto key = TypeIndexKey{g1 + (g2 << 16)};
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

Poly PiGroebner::Reduce(Poly poly) const
{
    Poly tmp_prod, tmp;
    size_t index = 0;
    while (index < poly.data.size() && !poly.data[index].IsUnKnown()) {
        if (poly.data[index].IsTrivial(gen_2tor_degs())) {
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
                poly.isubmulP(q, data_[gb_index], tmp, gen_2tor_degs());
            else
                poly.iaddmulP(q, data_[gb_index], tmp, gen_2tor_degs());
        }
        else
            ++index;
    }
    return poly;
}

Poly PiGroebner::ReduceForGbRel(Poly poly) const
{
    Poly tmp_prod, tmp;
    size_t index = 0;
    while (index < poly.data.size() && !poly.data[index].IsUnKnown()) {
        if (poly.data[index].IsTrivial(gen_2tor_degs())) {
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
                poly.isubmulP(q, data_[gb_index], tmp, gen_2tor_degs());
            else
                poly.iaddmulP(q, data_[gb_index], tmp, gen_2tor_degs());
        }
        else {
            if (index > 0 && MultipleOf(poly.data[index], poly.data[0])) {
                int c = poly.data[index].c() - poly.data[0].c();
                Poly poly1 = poly;
                poly.isubmulP(Mon::twoTo(c), poly1, tmp, gen_2tor_degs());
            }
            else
                ++index;
        }
    }
    return poly;
}

Poly PiGroebner::ReduceV2(Poly poly) const
{
    Poly tmp_prod, tmp;
    size_t index = 0;
    while (index < poly.data.size() && !poly.data[index].IsUnKnown()) {
        if (poly.data[index].IsTrivial(gen_2tor_degs())) {
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
                poly.isubmulP(q, data_[gb_index], tmp, gen_2tor_degs());
            else
                poly.iaddmulP(q, data_[gb_index], tmp, gen_2tor_degs());
        }
        else
            ++index;
    }
    return poly;
}

void PiGroebner::AddRel(Poly rel, int t_max, const ut::map_seq2d<int, 0>& possEinf)
{
    auto deg = GetDeg(rel.GetLead(), gen_degs_);
    Poly prod, tmp;
    for (auto& [d, pi_basis_d] : nodes_basis_.front()) {
        if (deg.t + d.t > t_trunc_)
            break;
        auto deg_prod = deg + d;
        for (auto& m : GetRecentPiBasis(d)->pi_basis) {
            Mon prod_lead = mul_unsigned(m, rel.GetLead());
            if (!ut::has(GetRecentPiBasis(deg_prod)->pi_basis, prod_lead)) {
                prod.data.clear();
                prod.iaddmulP(m, rel, tmp, gen_2tor_degs());
                prod = ReduceForGbRel(std::move(prod));
                if (IsValidRel(prod)) {
                    deg_prod = GetDeg(prod.GetLead(), gen_degs_);
                    if (deg_prod.t <= t_trunc_) {
                        if (t_max == t_trunc_ && deg_prod.t > t) {
                            /*if (p.Str() == "x_8x_{51}+x_{13}x_{43}+O(29)") {
                                std::cout << (p.data[0] < p.data[1]) << '\n';
                                std::cout << "debug\n";
                            }*/
                            rels_graded[deg.t].push_back(std::move(rel));
                        }
                        else {
                            /*if (p.Str() == "x_8x_{51}+x_{13}x_{43}+O(29)") {
                                std::cout << (p.data[0] < p.data[1]) << '\n';
                                std::cout << "debug\n";
                            }*/
                            ExtendO(possEinf, d_trunc, deg.stem(), rel);
                            push_back_data(std::move(rel), deg);
                        }
                    }
                }
            }
        }
    }
}

void PiGroebner::AddRels(Poly1d rels, int t_max, const ut::map_seq2d<int, 0>& possEinf)
{
    int d_trunc = deg_trunc();
    if (t_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");
    Poly prod, tmp;

    Poly2d rels_graded(t_max + 1), data_graded(t_max + 1);
    for (auto& rel : rels) {
        if (IsValidRel(rel)) {
            int t = GetDegT(rel.GetLead(), gen_degs_);
            if (t <= t_max)
                rels_graded[t].push_back(std::move(rel));
        }
    }

    std::vector<AdamsDeg1d> degs_group_by_t(t_max + 1);
    for (auto& [deg, _] : nodes_basis_.front())
        if (deg.t <= t_max)
            degs_group_by_t[deg.t].push_back(deg);

    for (int t = 1; t <= t_max; ++t) {
        for (int t1 = 1; t1 <= t; ++t1) {
            int t2 = t - t1;
            for (AdamsDeg d : degs_group_by_t[t2]) {
                for (auto& m : GetRecentPiBasis(d)->pi_basis) {
                    for (auto& rel : rels_graded[t1]) {
                        Mon prod_lead = mul_unsigned(m, rel.GetLead());
                        int eff_min = rel.UnknownFil() > FIL_MAX ? FIL_MAX + 1 : rel.UnknownFil() - rel.data[0].fil();
                        int gb_index = IndexOfDivisibleLeading(prod_lead, eff_min);
                        if (gb_index != -1) {
                            prod.data.clear();
                            prod.iaddmulP(m, rel, tmp, gen_2tor_degs());
                            data_graded[t].push_back(std::move(prod));
                        }
                    }
                }
            }
        }

        auto& data_graded_t = data_graded[t];
        ut::for_each_par128(data_graded_t.size(), [this, &data_graded_t](size_t i) { data_graded_t[i] = Reduce(std::move(data_graded_t[i])); });
        std::sort(data_graded_t.begin(), data_graded_t.end(), [](const Poly& p1, const Poly& p2) { return p1.UnknownFil() > p2.UnknownFil(); });

        for (auto& rel : data_graded_t) {
            rel = ReduceForGbRel(std::move(rel));
            if (IsValidRel(rel)) {
                AdamsDeg deg = GetDeg(rel.GetLead(), gen_degs_);
                if (deg.t <= d_trunc) {
                    if (deg.t > t) {
                        /*if (p.Str() == "x_8x_{51}+x_{13}x_{43}+O(29)") {
                            std::cout << (p.data[0] < p.data[1]) << '\n';
                            std::cout << "debug\n";
                        }*/
                        data_graded[deg.t].push_back(std::move(rel));
                    }
                    else {
                        /*if (p.Str() == "x_8x_{51}+x_{13}x_{43}+O(29)") {
                            std::cout << (p.data[0] < p.data[1]) << '\n';
                            std::cout << "debug\n";
                        }*/
                        bool changed = ExtendO(possEinf, d_trunc, deg.stem(), rel);
                        if(changed)
                            rels_graded[t].push_back(rel);
                        push_back_data(std::move(rel), deg);
                    }
                }
            }
        }
    }
}

void PiGroebner::SimplifyRels(const ut::map_seq2d<int, 0>& possEinf)
{
    Poly1d data1 = std::move(data_);
    ResetRels();
    AddRels(std::move(data1), criticals_.deg_trunc(), possEinf);
}

/**
 * Sort the sequence and combine the coefficients of monomials
 */
template <typename T>
void Merge1(std::vector<T>& data, std::vector<T>& tmp)
{
    auto cmp = [](const T& m1, const T& m2) { return m2 < m1; };
    tmp.clear();
    std::make_heap(data.begin(), data.end(), cmp);
    while (!data.empty()) {
        if (data.front().IsUnKnown()) {
            tmp.push_back(data.front());
            break;
        }
        std::pop_heap(data.begin(), data.end(), cmp);
        if (data.size() > 1 && data.back() == data.front()) {
            data.pop_back();
            std::pop_heap(data.begin(), data.end(), cmp);
            if (data.back().Is2Torsion() || data.back().fil() >= FIL_MAX)
                data.pop_back();
            else {
                data.back().imul2();
                std::push_heap(data.begin(), data.end(), cmp);
            }
        }
        else {
            tmp.push_back(data.back());
            data.pop_back();
        }
    }
    ut::copy(tmp, data);
}

void PiGroebner::SimplifyRelsReorder(const ut::map_seq2d<int, 0>& possEinf)
{
    Poly1d data1 = std::move(data_);
    ResetRels();
    Mon1d tmp;
    for (auto& p : data1)
        Merge1(p.data, tmp);
    AddRels(std::move(data1), criticals_.deg_trunc(), possEinf);
}

void PiGroebner::Sync()
{
    int t_max = ssS0_.t_max;
    auto& gb = ssS0_.gb;
    auto& basis = ssS0_.basis;
    auto& nodes_ss = ssS0_.nodes_ss;
    auto& pi_gb = ssS0_.pi_gb;
    auto& pi_basis = ssS0_.pi_basis;
    auto& pi_gen_Einf = ssS0_.pi_gen_Einf;

    int count_ss_old = count_ss;
    auto deduce_out = myio::Logger::fout2();

    Poly tmp;
    for (int t = 1; t <= t_max; ++t) {
        for (int s = 0; s <= t; ++s) {
            AdamsDeg deg(s, t);
            if (deg.stem() < deg_min.stem() || deg.s < deg_min.s)
                continue;
            algZ::Mon1d pi_basis_d = GenBasis(pi_gb, deg, pi_basis);
            if (basis.find(deg) != basis.end()) {
                /* Add new boundaries to ss */
                int2d boundaries = GetS0GbEinf(deg);
                const int r = deg.s + 1;
                const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                for (int1d& boundary : boundaries) {
                    if (!boundaries.empty()) {
                        const int count1 = SetS0DiffLeibnizV2(deg_src, {}, boundary, r);
                        if (count1 > 0 && depth == 0) {
                            std::cout << "\033[38;2;128;128;255mHomotopy to SS:  S0  " << deg.StrAdams() << "  " << boundary << " is a boundary\n\033[0m";
                            deduce_out << "Nontrivial: Homotopy to SS:  S0  " << deg.StrAdams() << "  " << boundary << " is a boundary\n";
                        }
                        count_ss += count1;
                    }
                }

                /* Construct the projection map */
                auto& basis_d = basis.at(deg);
                int2d projs;
                {
                    auto& sc = GetRecentStaircase(nodes_ss, deg);
                    size_t first_PC = GetFirstIndexOnLevel(sc, kLevelMax / 2);
                    for (auto& m : pi_basis_d) {
                        int1d proj = Poly2Indices(gb.Reduce(Proj(m, pi_gen_Einf)), basis_d);
                        proj = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_PC, proj);
                        projs.push_back(proj);
                    }
                }
                int2d image, kernel, g;
                lina::SetLinearMap(projs, image, kernel, g);

                /* Compute pi_basis */
                int1d non_leads = lina::AddVectors(ut::int_range((int)pi_basis_d.size()), lina::GetLeads(kernel));
                algZ::Mon1d pi_basis_d_v2;
                int2d pi_basis_d_Einf;
                for (int i : non_leads) {
                    pi_basis_d_v2.push_back(pi_basis_d[i]);
                    pi_basis_d_Einf.push_back(projs[i]);
                }

                /* Add new relations in homotopy */
                algZ::Poly1d pi_rels;
                for (auto& k : kernel) {
                    pi_rels.push_back(algZ::Indices2Poly(k, pi_basis_d) + algZ::Mon::O(deg.s + 1));
                    if (depth == 0) {
                        std::cout << "\033[38;2;220;255;220mSS to Homotopy:  S0  " << pi_rels.back() << "=0\n\033[0m";
                        deduce_out << "Nontrivial: SS to Homotopy:  S0  " << pi_rels.back() << "=0\n";
                    }
                }
                pi_gb.AddRels(pi_rels, t_max, ssS0_.basis_ss_possEinf);
                count_homotopy += (int)pi_rels.size();

                /* Add new generators in homotopy */
                int2d Einf;
                {
                    auto& sc = GetRecentStaircase(nodes_ss, deg);
                    size_t first_PC = GetFirstIndexOnLevel(sc, kLevelMax / 2);
                    size_t last_PC = GetFirstIndexOnLevel(sc, kLevelPC + 1);
                    for (size_t i = first_PC; i < last_PC; ++i)
                        Einf.push_back(sc.basis_ind[i]);
                }
                int2d new_generators = QuotientSpace(Einf, image); /* image is supposed to be a subspace of Einf */
                size_t gen_size_old = pi_gb.gen_degs().size();
                for (size_t i = 0; i < new_generators.size(); ++i) {
                    pi_gb.AddGen(deg);
                    pi_gen_Einf.push_back(Indices2Poly(new_generators[i], basis_d));
                    pi_basis_d_v2.push_back(algZ::Mon::Gen(uint32_t(gen_size_old + i), 1, deg.s, deg.stem() % 2 == 0));
                    pi_basis_d_Einf.push_back(new_generators[i]);
                    if (depth == 0) {
                        std::cout << "\033[38;2;128;255;128mSS to Homotopy:  S0  x_{" << std::to_string(pi_gb.gen_degs().size() - 1) << "} detected by " << pi_gen_Einf.back() << "\n\033[0m";
                        deduce_out << "Nontrivial: SS to Homotopy:  S0  x_{" << std::to_string(pi_gb.gen_degs().size() - 1) << "} detected by " << pi_gen_Einf.back() << "\n";
                    }
                }
                count_homotopy += (int)new_generators.size();

                if (deg.stem() % 2 == 1) {
                    algZ::Poly1d pi_trivial_rels;
                    for (size_t i = 0; i < new_generators.size(); ++i)
                        pi_trivial_rels.push_back(algZ::Mon::two_x_square(uint32_t(gen_size_old + i), deg.s));
                    pi_gb.AddRels(pi_trivial_rels, t_max, ssS0_.basis_ss_possEinf);
                }

                auto indices = ut::size_t_range(pi_basis_d_v2.size());
                std::sort(indices.begin(), indices.end(), [&pi_basis_d_v2](size_t i, size_t j) { return pi_basis_d_v2[i] < pi_basis_d_v2[j]; });
                pi_basis.back()[deg].pi_basis.clear();
                pi_basis.back()[deg].Einf.clear();
                for (auto i : indices) {
                    pi_basis.back().at(deg).pi_basis.push_back(pi_basis_d_v2[i]);
                    pi_basis.back().at(deg).Einf.push_back(std::move(pi_basis_d_Einf[i]));
                }
            }
            else {
                /* Add new relations in homotopy */
                algZ::Poly1d pi_rels;
                for (auto& m : pi_basis_d)
                    pi_rels.push_back(m + algZ::Mon::O(deg.s + 1));
                pi_gb.AddRels(pi_rels, t_max, ssS0_.basis_ss_possEinf);
            }
        }
    }
    if (count_ss_old != count_ss)
        UpdatePossEinf(ssS0_.nodes_ss, ssS0_.basis_ss_possEinf);
}

bool PiGroebner::IsNewBaseByLastGen(const Mon& mon, uint32_t last_gen) const
{
    if (last_gen >= leads_group_by_last_gen_.size() || std::none_of(leads_group_by_last_gen_[last_gen].begin(), leads_group_by_last_gen_[last_gen].end(), [this, &mon](int i) { return divisible(this->leads()[i], mon); }))
        return true;
    return false;
}

Poly1d PiGroebner::RelsLF(AdamsDeg deg) const
{
    Poly1d result;
    auto p = leads_group_by_deg_.find(deg);
    if (p != leads_group_by_deg_.end())
        for (int i : p->second)
            result.push_back(data_[i].LF());
    return result;
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

/********************************* Modules ****************************************/

PiGroebnerMod::PiGroebnerMod(PiGroebner* pGb, int deg_trunc, AdamsDeg1d v_degs, Mod1d polys, bool bDynamic) : pGb_(pGb), criticals_(deg_trunc), v_degs_(std::move(v_degs)), old_pGb_size_(pGb->leads_.size())
{
    for (auto& p : polys) {
        AdamsDeg deg = GetDeg(p.GetLead(), pGb_->gen_degs(), v_degs_);
        push_back_data_init(std::move(p), deg);
    }
}

bool PiGroebnerMod::IsNewBaseByV(const Mon& mon, uint32_t v) const
{
    if (v >= leads_group_by_v_.size() || std::none_of(leads_group_by_v_[v].begin(), leads_group_by_v_[v].end(), [this, &mon](int i) { return divisible(this->leads()[i].m, mon); }))
        return true;
    return false;
}

Mod1d PiGroebnerMod::RelsLF(AdamsDeg deg) const
{
    Mod1d result;
    auto p = leads_group_by_deg_.find(deg);
    if (p != leads_group_by_deg_.end())
        for (int i : p->second)
            result.push_back(data_[i].LF());
    return result;
}

void PiGroebnerMod::Pop(size_t gen_size, size_t rel_size)
{
    old_pGb_size_ = pGb_->leads_.size();
    v_degs_.resize(gen_size);

    data_.resize(rel_size);
    leads_.resize(rel_size);
    traces_.resize(rel_size);
    leads_O_.resize(rel_size);
    for (auto& [_, indices] : leads_group_by_key_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& [_, indices] : leads_group_by_deg_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& [_, indices] : leads_group_by_t_)
        pop_indices_by_ub(indices, (int)rel_size);
    for (auto& indices : leads_group_by_v_)
        pop_indices_by_ub(indices, (int)rel_size);
}

void PiGroebnerMod::debug_print() const
{
    std::cout << "v_degs_.size()=" << v_degs_.size() << '\n';
    std::cout << "data_.size()=" << data_.size() << '\n';
    size_t sum = 0;
    for (auto& [_, indices] : leads_group_by_key_)
        sum += indices.size();
    std::cout << "leads_group_by_key_=" << sum << '\n';

    sum = 0;
    for (auto& [_, indices] : leads_group_by_deg_)
        sum += indices.size();
    std::cout << "leads_group_by_deg_=" << sum << '\n';

    sum = 0;
    for (auto& [_, indices] : leads_group_by_t_)
        sum += indices.size();
    std::cout << "leads_group_by_t_=" << sum << '\n';

    sum = 0;
    for (auto& indices : leads_group_by_v_)
        sum += indices.size();
    std::cout << "leads_group_by_v_=" << sum << '\n';
}

int PiGroebnerMod::IndexOfDivisibleLeading(const MMod& mon, int eff_min) const
{
    const auto t = mon.m.Trace();
    const int i_max = int(mon.m.m().size());
    for (int i = mon.m.c() > 0 ? -2 : -1; i < i_max; ++i) {
        uint32_t backg = 0;
        if (i == -1 && mon.m.c() > 0)
            backg = 1;
        else if (i >= 0)
            backg = mon.m.m()[i].g() + 1;
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

int PiGroebnerMod::IndexOfDivisibleLeadingV2(const MMod& mon) const
{
    int result = -1, eff = 0;
    auto t = mon.m.Trace();
    const int i_max = int(mon.m.m().size());
    for (int i = mon.m.c() > 0 ? -2 : -1; i < i_max; ++i) {
        uint32_t backg = 0;
        if (i == -1 && mon.m.c() > 0)
            backg = 1;
        else if (i >= 0)
            backg = mon.m.m()[i].g() + 1;
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

Mod PiGroebnerMod::Reduce(Mod x) const
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
            else
                ++index;
        }
    }
    return x;
}

Mod PiGroebnerMod::ReduceForGbRel(Mod x) const
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

Mod PiGroebnerMod::ReduceV2(Mod x) const
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

void PiGroebnerMod::AddRels(Mod1d rels, int t_max, const ut::map_seq2d<int, 0>& possEinf)
{
    int d_trunc = deg_trunc();
    if (t_max > d_trunc)
        throw MyException(0x42e4ce5dU, "deg is bigger than the truncation degree.");
    if (old_pGb_size_ != pGb_->leads_.size()) {
        criticals_.AddToBuffersX(pGb_->leads_, pGb_->traces_, pGb_->leads_O_, leads_, traces_, leads_O_, pGb_->gen_degs(), v_degs_, old_pGb_size_);
        old_pGb_size_ = pGb_->leads_.size();
    }

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, Mod1d> rels_graded;
    for (auto& rel : rels) {
        if (rel) {
            AdamsDeg d = GetDeg(rel.GetLead(), pGb_->gen_degs(), v_degs_);
            if (d.t <= t_max)
                rels_graded[d.stem()].push_back(std::move(rel));
        }
    }

    Mod tmp;
    for (int t = 0; t <= t_max && ((!rels_graded.empty() && t <= rels_graded.rbegin()->first) || !criticals_.empty()); ++t) {
        int next_stem = criticals_.NextD();
        if (next_stem != -1 && next_stem < t) {
            if (next_stem < t - 1)
                throw MyException(0x3dc040d6U, "buffer_min_pairs_ contains degree < t - 1");
            --t;
        }

        CriPair1d pairs_d = Criticals(t);
        size_t pairs_d_size = pairs_d.size();
        auto& rels_d = rels_graded[t];
        Mod1d rels_tmp(pairs_d_size + rels_d.size());
        for (size_t i = 0; i < pairs_d.size(); ++i)
            pairs_d[i].SijMod(*pGb_, *this, rels_tmp[i], tmp);
        for (size_t i = 0; i < rels_d.size(); ++i)
            rels_tmp[pairs_d_size + i] = std::move(rels_d[i]);
        ut::for_each_par128(rels_tmp.size(), [this, &rels_tmp](size_t i) { rels_tmp[i] = Reduce(std::move(rels_tmp[i])); });
        std::sort(rels_tmp.begin(), rels_tmp.end(), [](const Mod& p1, const Mod& p2) { return p1.UnknownFil() > p2.UnknownFil(); });

        for (auto& rel : rels_tmp) {
            rel = ReduceForGbRel(std::move(rel));
            if (rel && !rel.GetLead().IsUnKnown()) {
                AdamsDeg deg = GetDeg(rel.GetLead(), pGb_->gen_degs(), v_degs_);
                if (deg.t <= d_trunc) {
                    if (t == d_trunc && deg.t > t) {
                        /*if (p.Str() == "x_5v_6+O(16)" || p.Str() == "x_1v_{13}")
                            std::cout << "debug\n";*/
                        rels_graded[deg.t].push_back(std::move(rel));
                    }
                    else {
                        /*if (p.Str() == "x_8x_{51}+x_{13}x_{43}+O(29)") {
                            std::cout << (p.data[0] < p.data[1]) << '\n';
                            std::cout << "debug\n";
                        }*/
                        ExtendO(possEinf, d_trunc, deg.stem(), rel);
                        push_back_data(std::move(rel), deg);
                    }
                }
            }
        }
    }
}

void PiGroebnerMod::SimplifyRels(const ut::map_seq2d<int, 0>& possEinf)
{
    Mod1d data1 = std::move(data_);
    ResetRels();
    AddRels(std::move(data1), criticals_.deg_trunc(), possEinf);
}

void PiGroebnerMod::SimplifyRelsReorder(const ut::map_seq2d<int, 0>& possEinf)
{
    Mod1d data1 = std::move(data_);
    ResetRels();
    MMod1d tmp;
    for (auto& p : data1)
        Merge1(p.data, tmp);
    AddRels(std::move(data1), criticals_.deg_trunc(), possEinf);
}

}  // namespace algZ
