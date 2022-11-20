#include "algebras/linalg.h"
#include "main.h"

Poly Proj(const algZ::Mon& mon, const Poly1d& map)
{
    Poly tmp;
    Poly result = pow(map[0], mon.c());
    result.imulP(subs(mon.m0(), map), tmp);
    result.imulP(subs(mon.m1(), map), tmp);
    return result;
}

Mod Proj(const algZ::MMod& mon, const Poly1d& map, const Mod1d& map_v)
{
    return Proj(mon.m, map) * map_v[mon.v];
}

int NextSExt(const Staircases1d& basis_ss, int t_max, int stem, int s_min)
{
    int s_max = (stem + 3) / 2;
    for (int s = s_min; s <= s_max; ++s) {
        if (stem + s > t_max || Diagram::PossEinf(basis_ss, AdamsDeg(s, stem + s)))
            return s;
    }
    return algZ::FIL_MAX + 1;
}

/* This is for use in GetKernel */
int NextSExtV2(const Staircases1d& basis_ss, int t_max, int stem, int s_min, std::unordered_map<int, int>& num_leads)
{
    int s_max = (stem + 3) / 2;
    for (int s = s_min; s <= s_max; ++s) {
        if (stem + s > t_max || Diagram::PossEinf(basis_ss, AdamsDeg(s, stem + s)) > num_leads[s])
            return s;
    }
    return algZ::FIL_MAX + 1;
}

template <typename T>
int ExtendRel(const Staircases1d& basis_ss, int stem, int t_max, const T& rel, T& rel_extended)
{
    auto& m = rel.data.back();
    if (m.IsUnKnown()) {
        int s = NextSExt(basis_ss, t_max, stem, m.fil());
        if (s > m.fil()) {
            rel_extended = rel;
            if (s <= algZ::FIL_MAX)
                rel_extended.data.back() = T::O(s);
            else
                rel_extended.data.pop_back();
            return 1;
        }
    }
    return 0;
}

template <typename T>
int ExtendRelV2(const Staircases1d& basis_ss, int stem, int t_max, const T& rel, T& rel_extended, std::unordered_map<int, int>& num_leads)
{
    auto& m = rel.data.back();
    if (m.IsUnKnown()) {
        int s = NextSExtV2(basis_ss, t_max, stem, m.fil(), num_leads);
        if (s > m.fil()) {
            rel_extended = rel;
            if (s <= algZ::FIL_MAX)
                rel_extended.data.back() = T::O(s);
            else
                rel_extended.data.pop_back();
            return 1;
        }
    }
    return 0;
}

template <typename T1, typename T2, typename GB1, typename GB2>
void GetKernel(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2, std::vector<T1>& kernel)
{
    std::vector<T2> image;
    std::vector<T1> g;
    T1 tmp1{};
    T2 tmp2{};
    /* f(g[i]) = image[i] */
    for (size_t i = 0; i < fx.size(); ++i) {
        T1 src = x[i];
        T2 tgt = fx[i];
        {
            size_t index = 0;
            while (index < tgt.data.size()) {
                size_t d_index = 1;
                for (size_t j = 0; j < image.size(); j++) {
                    if (MultipleOf(tgt.data[index], image[j].GetLead())) {
                        int c = tgt.data[index].c() - image[j].GetLead().c();
                        tgt.isubmulP(algZ::Mon::twoTo(c), image[j], tmp2, gb2.gen_2tor_degs());
                        tgt = gb2.Reduce(std::move(tgt));
                        src.isubmulP(algZ::Mon::twoTo(c), g[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.Reduce(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
        }
        if (tgt) {
            image.push_back(std::move(tgt));
            g.push_back(std::move(src));
        }
        else {
            size_t index = 0;
            while (index < src.data.size()) {
                size_t d_index = 1;
                for (size_t j = 0; j < kernel.size(); j++) {
                    if (MultipleOf(src.data[index], kernel[j].GetLead())) {
                        int c = src.data[index].c() - kernel[j].GetLead().c();
                        src.isubmulP(algZ::Mon::twoTo(c), kernel[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.Reduce(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
            if (src)
                kernel.push_back(std::move(src));
        }
    }
}

template <typename T1, typename T2, typename GB1, typename GB2>
void GetImage(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2, std::vector<T2>& image, int& O1, int& O2, T1& gO)
{
    /* Sort by certainty */
    auto indices = ut::size_t_range(fx.size());
    std::stable_sort(indices.begin(), indices.end(), [&fx](size_t i, size_t j) { return fx[i].UnknownFil() > fx[j].UnknownFil(); });

    std::vector<T1> g;
    T1 tmp1{};
    T2 tmp2{};
    /* f(g[i]) = image[i] */
    for (size_t i : indices) {
        T1 src = x[i];
        T2 tgt = fx[i];
        {
            size_t index = 0;
            while (index < tgt.data.size() && !tgt.data[index].IsUnKnown()) {
                size_t d_index = 1;
                for (size_t j = 0; j < image.size(); j++) {
                    if (MultipleOf(tgt.data[index], image[j].GetLead())) {
                        int c = tgt.data[index].c() - image[j].GetLead().c();
                        tgt.isubmulP(algZ::Mon::twoTo(c), image[j], tmp2, gb2.gen_2tor_degs());
                        tgt = gb2.ReduceV2(std::move(tgt));
                        src.isubmulP(algZ::Mon::twoTo(c), g[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.ReduceV2(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
        }
        if (tgt) {
            if (!tgt.GetLead().IsUnKnown()) {
                image.push_back(std::move(tgt));
                g.push_back(std::move(src));
            }
            else {
                int O = tgt.UnknownFil();
                if (O < O1) {
                    O2 = O1;
                    O1 = O;
                    gO = std::move(src);
                }
                else if (O == O1) {
                    O2 = O1;
                    gO = T1::O(0);
                }
                else if (O < O2)
                    O2 = O;
            }
        }
    }
}

/* Group the leads by s and return the numbers of leads for each s */
template <typename T1, typename T2, typename GB1, typename GB2>
auto GetImageLeads(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2) -> std::unordered_map<int, int>
{
    std::unordered_map<int, int> result;

    /* Sort by certainty */
    auto indices = ut::size_t_range(fx.size());
    std::stable_sort(indices.begin(), indices.end(), [&fx](size_t i, size_t j) { return fx[i].UnknownFil() > fx[j].UnknownFil(); });

    /* f(g[i]) = image[i] */
    std::vector<T2> image;
    std::vector<T1> g;
    T1 tmp1{};
    T2 tmp2{};
    for (size_t i : indices) {
        T1 src = x[i];
        T2 tgt = fx[i];
        {
            size_t index = 0;
            while (index < tgt.data.size() && !tgt.data[index].IsUnKnown()) {
                size_t d_index = 1;
                for (size_t j = 0; j < image.size(); j++) {
                    if (tgt.data[index] == image[j].GetLead()) {
                        tgt.isubP(image[j], tmp2, gb2.gen_2tor_degs());
                        tgt = gb2.ReduceV2(std::move(tgt));
                        src.isubP(g[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.ReduceV2(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
        }
        if (tgt) {
            if (!tgt.GetLead().IsUnKnown()) {
                ++result[tgt.GetLead().fil()];
                image.push_back(std::move(tgt));
                g.push_back(std::move(src));
            }
        }
    }
    return result;
}

template <typename T, typename GB>
T Residue(const std::vector<T>& space, T x, const GB& gb)
{
    T tmp{};
    x = gb.ReduceV2(std::move(x));
    size_t index = 0;
    while (index < x.data.size() && !x.data[index].IsUnKnown()) {
        size_t d_index = 1;
        for (size_t j = 0; j < space.size(); j++) {
            if (MultipleOf(x.data[index], space[j].GetLead())) {
                int c = x.data[index].c() - space[j].GetLead().c();
                x.isubmulP(algZ::Mon::twoTo(c), space[j], tmp, gb.gen_2tor_degs());
                x = gb.ReduceV2(std::move(x));
                d_index = 0;
                break;
            }
        }
        index += d_index;
    }
    return x;
}

int Diagram::PossEinf(const Staircases1d& basis_ss, AdamsDeg deg)
{
    if (basis_ss.front().find(deg) != basis_ss.front().end()) {
        const Staircase& sc = GetRecentStaircase(basis_ss, deg);
        size_t i_start_perm = GetFirstIndexOnLevel(sc, kLevelMax / 2);
        size_t i_stable = GetFirstIndexOfFixedLevels(basis_ss, deg, kLevelPC + 1);
        return int(i_stable - i_start_perm);
    }
    return 0;
}

int Diagram::PossMoreEinf(const Staircases1d& basis_ss, AdamsDeg deg)
{
    if (basis_ss.front().find(deg) != basis_ss.front().end()) {
        const Staircase& sc = GetRecentStaircase(basis_ss, deg);
        size_t i_end_perm = GetFirstIndexOnLevel(sc, kLevelPC + 1);
        size_t i_stable = GetFirstIndexOfFixedLevels(basis_ss, deg, kLevelPC + 1);
        return int(i_stable - i_end_perm);
    }
    return 0;
}

int1d Diagram::PossMoreEinfFirstS_S0() const
{
    int1d result;
    for (int i = 0; i <= ssS0_.t_max; ++i)
        result.push_back(ssS0_.t_max - i);
    for (auto& [deg, _] : ssS0_.basis_ss.front()) {
        if (PossMoreEinf(ssS0_.basis_ss, deg)) {
            if (deg.s < result[deg.stem()])
                result[deg.stem()] = deg.s;
        }
    }
    return result;
}

int1d Diagram::PossMoreEinfFirstS_Cof(size_t iCof) const
{
    int1d result;
    for (int i = 0; i <= ssCofs_[iCof].t_max; ++i)
        result.push_back(ssCofs_[iCof].t_max - i);
    for (auto& [deg, _] : ssCofs_[iCof].basis_ss.front()) {
        if (PossMoreEinf(ssCofs_[iCof].basis_ss, deg)) {
            if (deg.s < result[deg.stem()])
                result[deg.stem()] = deg.s;
        }
    }
    return result;
}

int Diagram::ExtendRelS0(int stem, const algZ::Poly& rel, algZ::Poly& rel_extended) const
{
    return ExtendRel(ssS0_.basis_ss, stem, ssS0_.t_max, rel, rel_extended);
}

int Diagram::ExtendRelCof(size_t iCof, int stem, const algZ::Mod& rel, algZ::Mod& rel_extended) const
{
    return ExtendRel(ssCofs_[iCof].basis_ss, stem, ssCofs_[iCof].t_max, rel, rel_extended);
}

void Diagram::ExtendRelS0(int stem, algZ::Poly& rel) const
{
    if (rel) {
        algZ::Poly rel1;
        if (ExtendRelS0(stem, rel, rel1))
            rel = std::move(rel1);
    }
}

void Diagram::ExtendRelCof(size_t iCof, int stem, algZ::Mod& rel) const
{
    if (rel) {
        algZ::Mod rel1;
        if (ExtendRelCof(iCof, stem, rel, rel1))
            rel = std::move(rel1);
    }
}

void Diagram::ExtendRelS0V2(int stem, algZ::Poly& rel, std::unordered_map<int, int>& num_leads) const
{
    if (rel) {
        algZ::Poly rel1;
        if (ExtendRelV2(ssS0_.basis_ss, stem, ssS0_.t_max, rel, rel1, num_leads))
            rel = std::move(rel1);
    }
}

void Diagram::ExtendRelCofV2(size_t iCof, int stem, algZ::Mod& rel, std::unordered_map<int, int>& num_leads) const
{
    if (rel) {
        algZ::Mod rel1;
        if (ExtendRelV2(ssCofs_[iCof].basis_ss, stem, ssCofs_[iCof].t_max, rel, rel1, num_leads))
            rel = std::move(rel1);
    }
}

int2d Diagram::GetS0GbEinf(AdamsDeg deg) const
{
    auto& gb = ssS0_.gb;
    auto& pi_gb = ssS0_.pi_gb;
    auto& pi_gen_Einf = ssS0_.pi_gen_Einf;

    int2d result;
    Poly tmp;
    algZ::Poly1d gb_Einf = pi_gb.RelsLF(deg);
    for (auto& f_Einf : gb_Einf) {
        Poly LF;
        for (auto& m : f_Einf.data)
            LF.iaddP(Proj(m, pi_gen_Einf), tmp);
        LF = gb.Reduce(std::move(LF));
        if (LF)
            result.push_back(Poly2Indices(LF, ssS0_.basis.at(deg)));
        else
            result.push_back(int1d{});
    }

    return result;
}

int2d Diagram::GetCofGbEinf(int iCof, AdamsDeg deg) const
{
    auto& gb = ssCofs_[iCof].gb;
    auto& pi_gb = ssCofs_[iCof].pi_gb;
    auto& pi_gen_Einf = ssCofs_[iCof].pi_gen_Einf;

    int2d result;
    Mod tmp;
    algZ::Mod1d gb_Einf = pi_gb.RelsLF(deg);
    for (auto& f_Einf : gb_Einf) {
        Mod LF;
        for (auto& m : f_Einf.data)
            LF.iaddP(Proj(m, ssS0_.pi_gen_Einf, pi_gen_Einf), tmp);
        LF = gb.Reduce(std::move(LF));
        if (LF)
            result.push_back(Mod2Indices(LF, ssCofs_[iCof].basis.at(deg)));
        else
            result.push_back(int1d{});
    }

    return result;
}

std::map<AdamsDeg, int2d> Diagram::GetS0GbEinf() const
{
    std::map<AdamsDeg, int2d> result;
    for (auto& [deg, _] : ssS0_.pi_gb.leads_group_by_deg())
        result[deg] = GetS0GbEinf(deg);
    return result;
}

std::map<AdamsDeg, int2d> Diagram::GetCofGbEinf(int iCof) const
{
    std::map<AdamsDeg, int2d> result;
    for (auto& [deg, _] : ssCofs_[iCof].pi_gb.leads_group_by_deg())
        result[deg] = GetCofGbEinf(iCof, deg);
    return result;
}

void Diagram::AddPiRelsS0(algZ::Poly1d rels)
{
    ssS0_.pi_gb.AddRels(std::move(rels), ssS0_.t_max);
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
        ssCofs_[iCof].pi_gb.AddRels({}, ssCofs_[iCof].t_max);
}

void Diagram::AddPiRelsCof(size_t iCof, algZ::Mod1d rels)
{
    algZ::Poly1d relsS0;
    for (auto& rel : rels) {
        auto relS0 = algZ::subsMod(rel, ssCofs_[iCof].pi_f_top_cell.back(), ssCofs_[iCof].pi_gb.v_degs());
        if (algZ::IsValidRel(relS0))
            relsS0.push_back(std::move(relS0));
    }
    ssCofs_[iCof].pi_gb.AddRels(std::move(rels), ssCofs_[iCof].t_max);
    AddPiRelsS0(relsS0);
}

void Diagram::AddPiRelsCof2S0(size_t iCof)
{
    algZ::Poly1d relsS0;
    for (auto& rel : ssCofs_[iCof].pi_gb.data()) {
        auto relS0 = algZ::subsMod(rel, ssCofs_[iCof].pi_f_top_cell.back(), ssCofs_[iCof].pi_gb.v_degs());
        if (algZ::IsValidRel(relS0))
            relsS0.push_back(std::move(relS0));
    }
    AddPiRelsS0(relsS0);
}

void Diagram::SyncS0Homotopy(int& count_ss, int& count_homotopy, int depth)
{
    int t_max = ssS0_.t_max;
    auto& gb = ssS0_.gb;
    auto& basis = ssS0_.basis;
    auto& basis_ss = ssS0_.basis_ss;
    auto& pi_gb = ssS0_.pi_gb;
    auto& pi_basis = ssS0_.pi_basis;
    auto& pi_gen_Einf = ssS0_.pi_gen_Einf;

    Poly tmp;
    for (int t = 1; t <= t_max; ++t) {
        for (int s = 0; s <= t; ++s) {
            AdamsDeg deg(s, t);
            algZ::Mon1d pi_basis_d = pi_gb.GenBasis(deg, pi_basis);
            if (basis.find(deg) != basis.end()) {
                /* Add new boundaries to ss */
                int2d boundaries = GetS0GbEinf(deg);
                const int r = deg.s + 1;
                const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                for (int1d& boundary : boundaries) {
                    const int count1 = SetS0DiffLeibnizV2(deg_src, {}, boundary, r);
                    if (count1 > 0 && depth == 0)
                        std::cout << "Homotopy to SS:  iSS=0  " << deg.StrAdams() << "  " << boundary << " is a boundary           \n";
                    count_ss += count1;
                }

                /* Construct the projection map */
                auto& basis_d = basis.at(deg);
                int2d projs;
                {
                    auto& sc = GetRecentStaircase(basis_ss, deg);
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
                    if (depth == 0)
                        std::cout << "SS to Homotopy:  iSS=0  " << pi_rels.back() << "=0          \n";
                }
                pi_gb.AddRels(pi_rels, t_max);
                count_homotopy += (int)pi_rels.size();

                /* Add new generators in homotopy */
                int2d Einf;
                {
                    auto& sc = GetRecentStaircase(basis_ss, deg);
                    size_t first_PC = GetFirstIndexOnLevel(sc, kLevelMax / 2);
                    size_t last_PC = GetFirstIndexOnLevel(sc, kLevelPC + 1);
                    for (size_t i = first_PC; i < last_PC; ++i)
                        Einf.push_back(sc.basis_ind[i]);
                }
                int2d new_generators = lina::QuotientSpace(Einf, image); /* image is supposed to be a subspace of Einf */
                size_t gen_size_old = pi_gb.gen_degs().size();
                for (size_t i = 0; i < new_generators.size(); ++i) {
                    pi_gb.AddGen(deg);
                    pi_gen_Einf.push_back(Indices2Poly(new_generators[i], basis_d));
                    pi_basis_d_v2.push_back(algZ::Mon::Gen(uint32_t(gen_size_old + i), 1, deg.s, deg.stem() % 2 == 0));
                    pi_basis_d_Einf.push_back(new_generators[i]);
                    if (depth == 0)
                        std::cout << "SS to Homotopy:  iSS=0  x_{" << std::to_string(pi_gb.gen_degs().size() - 1) << "} detected by " << pi_gen_Einf.back() << "          \n";
                }
                count_homotopy += (int)new_generators.size();

                if (deg.stem() % 2 == 1) {
                    algZ::Poly1d pi_trivial_rels;
                    for (size_t i = 0; i < new_generators.size(); ++i)
                        pi_trivial_rels.push_back(algZ::Mon::two_x_square(uint32_t(gen_size_old + i), deg.s));
                    pi_gb.AddRels(pi_trivial_rels, t_max);
                }

                pi_basis[deg] = pi_basis_d_v2;
                ssS0_.pi_basis_Einf[deg] = pi_basis_d_Einf;
            }
            else {
                /* Add new relations in homotopy */
                algZ::Poly1d pi_rels;
                for (auto& m : pi_basis_d)
                    pi_rels.push_back(m + algZ::Mon::O(deg.s + 1));
                pi_gb.AddRels(pi_rels, t_max);
            }
        }
    }
}

void Diagram::SyncCofHomotopy(int iCof, int& count_ss, int& count_homotopy, int depth)
{
    int t_max = ssCofs_[iCof].t_max;
    auto& gb = ssCofs_[iCof].gb;
    auto& basis = ssCofs_[iCof].basis;
    auto& basis_ss = ssCofs_[iCof].basis_ss;
    auto& pi_gb = ssCofs_[iCof].pi_gb;
    auto& pi_basis = ssCofs_[iCof].pi_basis;
    auto& pi_gen_Einf = ssCofs_[iCof].pi_gen_Einf;
    auto& pi_f_top_cell = ssCofs_[iCof].pi_f_top_cell;

    Mod tmp;
    for (int t = 0; t <= t_max; ++t) {
        for (int s = 0; s <= t; ++s) {
            AdamsDeg deg(s, t);
            algZ::MMod1d pi_basis_d = pi_gb.GenBasis(deg, ssS0_.pi_basis);
            if (basis.find(deg) != basis.end()) {
                /* Add new boundaries to ss */
                int2d boundaries = GetCofGbEinf(iCof, deg);
                const int r = deg.s + 1;
                const AdamsDeg deg_src = deg - AdamsDeg(r, r - 1);
                for (auto& boundary : boundaries) {
                    const int count1 = SetCofDiffLeibnizV2(iCof, deg_src, {}, boundary, r);
                    if (count1 > 0 && depth == 0)
                        std::cout << "Homotopy to SS:  iSS=" << iCof + 1 << "  " << deg.StrAdams() << "  " << boundary << " is a boundary           \n";
                    count_ss += count1;
                }

                /* Construct the projection map */
                auto& basis_d = basis.at(deg);
                int2d projs;
                {
                    auto& sc = GetRecentStaircase(basis_ss, deg);
                    size_t first_PC = GetFirstIndexOnLevel(sc, kLevelMax / 2);
                    for (auto& m : pi_basis_d) {
                        int1d proj = Mod2Indices(gb.Reduce(Proj(m, ssS0_.pi_gen_Einf, pi_gen_Einf)), basis_d);
                        proj = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_PC, proj);
                        projs.push_back(proj);
                    }
                }
                int2d image, kernel, g;
                lina::SetLinearMap(projs, image, kernel, g);

                /* Compute pi_basis */
                int1d non_leads = lina::AddVectors(ut::int_range((int)pi_basis_d.size()), lina::GetLeads(kernel));
                algZ::MMod1d pi_basis_d_v2;
                int2d pi_basis_d_Einf;
                for (int i : non_leads) {
                    pi_basis_d_v2.push_back(pi_basis_d[i]);
                    pi_basis_d_Einf.push_back(projs[i]);
                }

                /* Add new relations in homotopy */
                algZ::Mod1d pi_rels;
                for (auto& k : kernel) {
                    pi_rels.push_back(algZ::Indices2Mod(k, pi_basis_d) + algZ::MMod::O(deg.s + 1));
                    if (depth == 0)
                        std::cout << "SS to Homotopy:  iSS=" << iCof + 1 << "  " << pi_rels.back() << "=0          \n";
                }
                pi_gb.AddRels(pi_rels, t_max);
                count_homotopy += (int)pi_rels.size();

                /* Add new generators in homotopy */
                int2d Einf;
                {
                    auto& sc = GetRecentStaircase(basis_ss, deg);
                    size_t first_PC = GetFirstIndexOnLevel(sc, kLevelMax / 2);
                    size_t last_PC = GetFirstIndexOnLevel(sc, kLevelPC + 1);
                    for (size_t i = first_PC; i < last_PC; ++i)
                        Einf.push_back(sc.basis_ind[i]);
                }
                int2d new_generators = lina::QuotientSpace(Einf, image); /* image is supposed to be a subspace of Einf */
                size_t gen_size_old = pi_gb.v_degs().size();
                for (size_t i = 0; i < new_generators.size(); ++i) {
                    pi_gb.AddGen(deg);
                    pi_gen_Einf.push_back(Indices2Mod(new_generators[i], basis_d));
                    pi_basis_d_v2.push_back(algZ::MMod({}, (uint32_t)(gen_size_old + i), deg.s));
                    pi_basis_d_Einf.push_back(new_generators[i]);
                    if (depth == 0)
                        std::cout << "SS to Homotopy:  iSS=" << iCof + 1 << "  v_{" << std::to_string(pi_gb.v_degs().size() - 1) << "} detected by " << pi_gen_Einf.back() << "          \n";

                    Mod x = Indices2Mod(new_generators[i], basis.at(deg));
                    Poly fx = ssS0_.gb.Reduce(subsMod(x, ssCofs_[iCof].f_top_cell));
                    AdamsDeg deg_S0 = deg - ssCofs_[iCof].deg_f_top_cell;
                    if (fx) {
                        int1d ifx = Poly2Indices(fx, ssS0_.basis.at(deg_S0));
                        auto& sc = GetRecentStaircase(ssS0_.basis_ss, deg_S0);
                        size_t first = GetFirstIndexOnLevel(sc, kLevelMax / 2);
                        ifx = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first, std::move(ifx));

                        int2d image, g;
                        lina::GetInvMap(ssS0_.pi_basis_Einf.at(deg_S0), image, g);
                        int1d pi_ifx = lina::GetImage(image, g, ifx);
                        algZ::Poly pi_fx = algZ::Indices2Poly(pi_ifx, ssS0_.pi_basis.at(deg_S0)) + algZ::Mon::O(deg_S0.s + 1);
                        pi_f_top_cell.back().push_back(std::move(pi_fx));
                    }
                    else
                        pi_f_top_cell.back().push_back(algZ::Mon::O(deg_S0.s + 1));
                }
                count_homotopy += (int)new_generators.size();

                pi_basis[deg] = pi_basis_d_v2;
                ssCofs_[iCof].pi_basis_Einf[deg] = pi_basis_d_Einf;
            }
            else {
                /* Add new relations in homotopy */
                algZ::Mod1d pi_rels;
                for (auto& m : pi_basis_d)
                    pi_rels.push_back(m + algZ::MMod::O(deg.s + 1));
                pi_gb.AddRels(pi_rels, t_max);
            }
        }
    }
}

int Diagram::DeduceTrivialExtensions(int depth)
{
    int count_homotopy = 0;

    /* top cell maps */
    {
        algZ::Poly1d new_rels_S0; /* By exactness fi * h = 0 */
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
            auto& pi_f = ssCofs_[iCof].pi_f_top_cell.back();
            for (size_t i = 0; i < pi_f.size(); ++i) {
                auto& fi = pi_f[i];
                AdamsDeg deg = ssCofs_[iCof].pi_gb.v_degs()[i] - ssCofs_[iCof].deg_f_top_cell;
                algZ::Poly fi_extended;
                if (fi && ExtendRelS0(deg.stem(), fi, fi_extended)) {
                    if (depth == 0)
                        std::cout << "For degree reason:  iSS=" << iCof + 1 << "  f" << std::to_string(i) << "=" << fi << " --> " << fi_extended << '\n';
                    fi = std::move(fi_extended);
                    ++count_homotopy;

                    algZ::Poly h = ssS0_.pi_gb.Gen((uint32_t)iCof);
                    algZ::Poly relS0 = ssS0_.pi_gb.Reduce(fi * h);
                    if (algZ::IsValidRel(relS0)) {
                        if (depth == 0)
                            std::cout << "    -> S0 rel: " << relS0 << '\n';
                        new_rels_S0.push_back(std::move(relS0));
                        ++count_homotopy;
                    }
                }
            }
        }
        AddPiRelsS0(std::move(new_rels_S0));
    }

    /* multiplicative structures */
    int old_count_homotopy = count_homotopy;
    std::vector<size_t> indices_start(all_basis_ss_.size(), 0);
    while (true) {
        { /* sphere */
            auto& pi_gb = ssS0_.pi_gb;
            size_t iSS = 0;

            algZ::Poly1d new_rels;
            size_t indices_end = pi_gb.data().size();
            for (size_t i = indices_start[iSS]; i < indices_end; ++i) {
                auto& rel = pi_gb.data()[i];
                AdamsDeg deg = GetDeg(rel.GetLead(), pi_gb.gen_degs());
                algZ::Poly rel_extended;
                if (ExtendRelS0(deg.stem(), rel, rel_extended)) {
                    rel_extended = pi_gb.Reduce(std::move(rel_extended));
                    if (algZ::IsValidRel(rel_extended)) {
                        if (depth == 0)
                            std::cout << "For degree reason:  iSS=" << iSS << "  " << rel << " --> " << rel_extended << '\n';
                        new_rels.push_back(std::move(rel_extended));
                        ++count_homotopy;
                    }
                }
            }
            indices_start[iSS] = indices_end;
            AddPiRelsS0(std::move(new_rels));
        }
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) { /* module */
            auto& pi_gb = ssCofs_[iCof].pi_gb;
            auto& basis_ss = ssCofs_[iCof].basis_ss;
            size_t iSS = iCof + 1;

            algZ::Mod1d new_rels;
            size_t indices_end = pi_gb.data().size();
            for (size_t i = indices_start[iSS]; i < indices_end; ++i) {
                auto& rel = pi_gb.data()[i];
                AdamsDeg deg = GetDeg(rel.GetLead(), ssS0_.pi_gb.gen_degs(), pi_gb.v_degs());
                algZ::Mod rel_extended;
                if (ExtendRelCof(iCof, deg.stem(), rel, rel_extended)) {
                    rel_extended = pi_gb.Reduce(std::move(rel_extended));
                    if (algZ::IsValidRel(rel_extended)) {
                        if (depth == 0)
                            std::cout << "For degree reason:  iSS=" << iSS << "  " << rel << " --> " << rel_extended << '\n';
                        new_rels.push_back(std::move(rel_extended));
                        ++count_homotopy;
                    }
                }
            }
            indices_start[iSS] = indices_end;
            AddPiRelsCof(iCof, std::move(new_rels));
        }

        if (old_count_homotopy != count_homotopy) {
            old_count_homotopy = count_homotopy;
            continue;
        }

        {
            auto& pi_basis = ssS0_.pi_basis;
            auto& pi_gb = ssS0_.pi_gb;
            const int t_max = ssS0_.t_max;
            algZ::Poly1d new_rels_S0;
            algZ::Mod2d new_rels_Cofs(ssCofs_.size());
            algZ::Poly tmp;
            algZ::Mod tmpm;
            for (auto& [deg, basis_d] : pi_basis) {
                for (auto& m : basis_d) {
                    for (uint32_t gen_id = 0; gen_id < 4; ++gen_id) {
                        int t_prod = deg.t + pi_gb.gen_degs()[gen_id].t;
                        if (t_prod < t_max) {
                            int stem_prod = deg.stem() + pi_gb.gen_degs()[gen_id].stem();
                            algZ::Poly prod;
                            auto h = pi_gb.Gen(gen_id);
                            prod.iaddmulP(h, m, tmp, pi_gb.gen_2tor_degs());
                            auto prod_reduced = pi_gb.ReduceV2(prod);
                            algZ::Poly prod_extended;
                            if (prod_reduced && ExtendRelS0(stem_prod, prod_reduced, prod_extended)) {
                                if (depth == 0)
                                    std::cout << "For degree reason *h" << gen_id << ":  iSS=" << 0 << "  " << h << '*' << m << '=' << prod_reduced << " --> " << prod_extended << '\n';
                                prod_extended.isubP(prod, tmp, pi_gb.gen_2tor_degs());
                                new_rels_S0.push_back(std::move(prod_extended));
                                ++count_homotopy;
                            }
                        }
                    }
                    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
                        if (deg.t < ssCofs_[iCof].t_max) {
                            algZ::Mod x(m, 0, 0);
                            auto x_reduced = ssCofs_[iCof].pi_gb.ReduceV2(x);
                            algZ::Mod x_extended;
                            if (x_reduced && ExtendRelCof(iCof, deg.stem(), x_reduced, x_extended)) {
                                if (depth == 0)
                                    std::cout << "For degree reason *v0:  iSS=" << iCof + 1 << "  " << x << '=' << x_reduced << " --> " << x_extended << '\n';
                                x_extended.isubP(x, tmpm, pi_gb.gen_2tor_degs());
                                new_rels_Cofs[iCof].push_back(std::move(x_extended));
                                ++count_homotopy;
                            }
                        }
                    }
                }
            }
            AddPiRelsS0(std::move(new_rels_S0));
            for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
                AddPiRelsCof(iCof, std::move(new_rels_Cofs[iCof]));
        }
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
            auto& pi_basis = ssCofs_[iCof].pi_basis;
            auto& pi_gb = ssCofs_[iCof].pi_gb;
            const int t_max = ssCofs_[iCof].t_max;
            algZ::Mod1d new_rels_cof;
            algZ::Mod tmp;
            for (auto& [deg, basis_d] : pi_basis) {
                for (auto& m : basis_d) {
                    for (uint32_t gen_id = 0; gen_id < 4; ++gen_id) {
                        int t_prod = deg.t + ssS0_.pi_gb.gen_degs()[gen_id].t;
                        if (t_prod < t_max) {
                            int stem_prod = deg.stem() + ssS0_.pi_gb.gen_degs()[gen_id].stem();
                            auto h = ssS0_.pi_gb.Gen(gen_id);
                            auto prod = h * m;
                            auto prod_reduced = pi_gb.ReduceV2(prod);
                            algZ::Mod prod_extended;
                            if (prod_reduced && ExtendRelCof(iCof, stem_prod, prod_reduced, prod_extended)) {
                                if (depth == 0)
                                    std::cout << "For degree reason v2:  iSS=" << iCof + 1 << "  " << h << '*' << m << "=" << prod_reduced << " --> " << prod_extended << '\n';
                                prod_extended.isubP(prod, tmp, ssS0_.pi_gb.gen_2tor_degs());
                                new_rels_cof.push_back(std::move(prod_extended));
                                ++count_homotopy;
                            }
                        }
                    }
                }
            }
            AddPiRelsCof(iCof, std::move(new_rels_cof));
        }

        if (old_count_homotopy != count_homotopy)
            old_count_homotopy = count_homotopy;
        else
            break;
    }
    return count_homotopy;
}

int Diagram::DeduceExtensionsByExactness(int depth)
{
    int count_homotopy = 0;
    algZ::Mon2d basis_S0(size_t(ssS0_.t_max + 1));
    std::vector<algZ::MMod2d> basis_Cofs(ssCofs_.size());
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
        basis_Cofs[iCof].resize(size_t(ssCofs_[iCof].t_max + 1));
    for (auto& [deg, basis_d] : ssS0_.pi_basis)
        for (auto& m : basis_d)
            basis_S0[deg.stem()].push_back(m);
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
        for (auto& [deg, basis_d] : ssCofs_[iCof].pi_basis)
            for (auto& m : basis_d)
                basis_Cofs[iCof][deg.stem()].push_back(m);
    int1d sPossMoreEinfS0 = PossMoreEinfFirstS_S0();
    algZ::Poly tmp;
    algZ::Mod tmpm;

    algZ::Poly1d new_rels_S0;
    algZ::Mod2d new_rels_Cofs(ssCofs_.size());
    int1d f_changed(ssCofs_.size(), 0);
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
        auto& ssCof = ssCofs_[iCof];
        auto& pi_f = ssCof.pi_f_top_cell;
        int d_f = ssCof.deg_f_top_cell.t;
        auto& h = ssS0_.pi_gb.Gen((uint32_t)iCof);

        int1d sPossMoreEinfCof = PossMoreEinfFirstS_Cof(iCof);

        /* exactness at h * f */
        int stem_max = std::min(ssS0_.t_max, ssCof.t_max - d_f);
        for (int stem = 0; stem <= stem_max; ++stem) {
            /* kernel of h */
            auto& basis_stem = basis_S0[stem];
            algZ::Poly1d x_h, hx;
            for (size_t i = 0; i < basis_stem.size(); ++i) {
                x_h.push_back(basis_stem[i]);
                hx.push_back(ssS0_.pi_gb.ReduceV2(h * basis_stem[i]));
            }
            auto num_leads = GetImageLeads(x_h, hx, ssS0_.pi_gb, ssS0_.pi_gb);
            for (size_t i = 0; i < basis_stem.size(); ++i) {
                auto hxi = ssS0_.pi_gb.ReduceV2(h * basis_stem[i]);
                bool use_hxi = false;
                if (hxi.UnknownFil() < algZ::FIL_MAX) {
                    ExtendRelS0V2(stem + d_f - 1, hxi, num_leads);
                    if (hxi.UnknownFil() > algZ::FIL_MAX) {
                        hx[i] = std::move(hxi);
                        use_hxi = true;
                    }
                }
                if (!use_hxi)
                    hx[i] = ssS0_.pi_gb.Reduce(h * basis_stem[i]);
            }
            algZ::Poly1d kernel_h;
            GetKernel(x_h, hx, ssS0_.pi_gb, ssS0_.pi_gb, kernel_h);

            /* image of f */
            algZ::Mod1d x_f;
            algZ::Poly1d fx;
            int stem1 = stem + d_f;
            auto& basis_cof_stem1 = basis_Cofs[iCof][stem1];
            for (size_t i = 0; i < basis_cof_stem1.size(); ++i) {
                x_f.push_back(basis_cof_stem1[i]);
                auto fxi = algZ::subsMod(basis_cof_stem1[i], ssCofs_[iCof].pi_f_top_cell.back(), ssCofs_[iCof].pi_gb.v_degs());
                fxi = ssS0_.pi_gb.ReduceV2(std::move(fxi));
                fx.push_back(std::move(fxi));
            }
            algZ::Poly1d image1_f;
            int O1 = sPossMoreEinfCof[stem1];
            int O2 = O1;
            algZ::Mod gO = algZ::Mod::O(0);
            GetImage(x_f, fx, ssCof.pi_gb, ssS0_.pi_gb, image1_f, O1, O2, gO);

            /* Check exactness */
            for (size_t i = 0; i < kernel_h.size(); ++i) {
                auto k = Residue(image1_f, std::move(kernel_h[i]), ssS0_.pi_gb);
                if (k && !k.GetLead().IsUnKnown()) {
                    int s = k.GetLead().fil();
                    AdamsDeg deg(s, stem + s);
                    if (!IsPossTgt(ssS0_.basis_ss, deg, kRPC)) {
                        if (s < O1) {
                            throw SSPiException(0x91866326U, "Top cell map can not hit the kernel of h");
                        }
                        if (s < O2) {
                            if (gO.data.size() == 1 && !gO.GetLead().IsUnKnown() && !gO.GetLead().m) {
                                size_t v = (size_t)gO.GetLead().v;
                                auto& f = ssCof.pi_f_top_cell.back()[v];
                                if (depth == 0)
                                    std::cout << "Exactness h * f:  iSS=" << iCof + 1 << "  stem=" << stem << "  f" << std::to_string(v) << "=" << f << " --> ";
                                f.data.pop_back();
                                f.iaddP(k.LF(), tmp);
                                f.data.push_back(algZ::Mon::O(s + 1));
                                ExtendRelS0(stem, f);
                                if (depth == 0)
                                    std::cout << f << '\n';
                                ++count_homotopy;
                                f_changed[iCof] = 1;
                            }
                        }
                    }
                    else if (s < O1) {
                        new_rels_S0.push_back(k);
                        if (depth == 0)
                            std::cout << "Exactness h * f:  iSS=" << 0 << "  stem=" << stem << k << "=0";
                        ++count_homotopy;
                    }

                    // std::cout << "stem=" << stem << "  " << k << "  O_k=" << k.GetLead().fil() << "  O1=" << O1 << "  O2=" << O2 << "  gO=" << gO << '\n';
                }
            }
        }

        /* exactness at i * h */
        for (int stem = d_f - 1; stem <= ssS0_.t_max; ++stem) {
            /* kernel of i */
            auto& basis_stem = basis_S0[stem];
            algZ::Poly1d x_i;
            algZ::Mod1d ix;
            for (size_t i = 0; i < basis_stem.size(); ++i) {
                x_i.push_back(basis_stem[i]);
                ix.push_back(ssCof.pi_gb.ReduceV2(algZ::Mod(basis_stem[i], 0, 0)));
            }
            auto num_leads = GetImageLeads(x_i, ix, ssS0_.pi_gb, ssCof.pi_gb);
            for (size_t i = 0; i < basis_stem.size(); ++i) {
                auto ixi = ssCof.pi_gb.ReduceV2(algZ::Mod(basis_stem[i], 0, 0));
                bool use_ixi = false;
                if (ixi.UnknownFil() < algZ::FIL_MAX) {
                    ExtendRelCofV2(iCof, stem, ixi, num_leads);
                    if (ixi.UnknownFil() > algZ::FIL_MAX) {
                        ix[i] = std::move(ixi);
                        use_ixi = true;
                    }
                }
                if (!use_ixi)
                    ix[i] = ssCof.pi_gb.Reduce(algZ::Mod(basis_stem[i], 0, 0));
            }
            algZ::Poly1d kernel_i;
            GetKernel(x_i, ix, ssS0_.pi_gb, ssCofs_[iCof].pi_gb, kernel_i);

            /* image of h */
            algZ::Poly1d x_h;
            algZ::Poly1d image_h;
            int stem1 = stem - d_f + 1;
            auto& basis_S0_stem1 = basis_S0[stem1];
            for (size_t i = 0; i < basis_S0_stem1.size(); ++i) {
                x_h.push_back(basis_S0_stem1[i]);
                auto fxi = ssS0_.pi_gb.ReduceV2(h * basis_S0_stem1[i]);
                image_h.push_back(std::move(fxi));
            }
            algZ::Poly1d image1_h;
            int O1 = sPossMoreEinfS0[stem] + 1;
            int O2 = O1;
            algZ::Poly gO = algZ::Mon::O(0);
            GetImage(x_h, image_h, ssS0_.pi_gb, ssS0_.pi_gb, image1_h, O1, O2, gO);

            /* Check exactness */
            for (size_t i = 0; i < kernel_i.size(); ++i) {
                auto k = Residue(image1_h, std::move(kernel_i[i]), ssS0_.pi_gb);
                if (k && !k.GetLead().IsUnKnown()) {
                    int s = k.GetLead().fil();
                    AdamsDeg deg(s, stem + s);
                    if (!IsPossTgt(ssS0_.basis_ss, deg, kRPC)) {
                        if (s < O1)
                            throw SSPiException(0x91866326U, "h multiples can not hit the kernel of i");
                        if (s < O2) {
                            if (gO && !gO.GetLead().IsUnKnown()) {
                                auto& rel = algZ::Poly(h) * gO + k.LF() + algZ::Poly::O(s + 1);
                                if (depth == 0)
                                    std::cout << "Exactness i * h:  iSS=" << iCof + 1 << "  stem=" << stem << "  " << rel << "=0\n";
                                new_rels_S0.push_back(std::move(rel));
                                ++count_homotopy;
                            }
                        }
                    }
                    else if (s < O1) {
                        new_rels_S0.push_back(k);
                        if (depth == 0)
                            std::cout << "Exactness i * h:  iSS=" << 0 << "  stem=" << stem << "  " << k << "=0";
                        ++count_homotopy;
                    }

                    // std::cout << "stem=" << stem << "  " << k << "  O_k=" << k.GetLead().fil() << "  O1=" << O1 << "  O2=" << O2 << "  gO=" << gO << '\n';
                }
            }
        }

        /* exactness at f * i */
        stem_max = std::min(ssS0_.t_max, ssCof.t_max);
        for (int stem = 0; stem <= stem_max; ++stem) {
            /* kernel of f */
            auto& basis_stem = basis_Cofs[iCof][stem];
            algZ::Mod1d x_f;
            algZ::Poly1d fx;
            for (size_t i = 0; i < basis_stem.size(); ++i) {
                x_f.push_back(basis_stem[i]);
                auto fxi = algZ::subsMod(basis_stem[i], ssCofs_[iCof].pi_f_top_cell.back(), ssCofs_[iCof].pi_gb.v_degs());
                fxi = ssS0_.pi_gb.ReduceV2(std::move(fxi));
                fx.push_back(fxi);
            }
            auto num_leads = GetImageLeads(x_f, fx, ssCof.pi_gb, ssS0_.pi_gb);
            x_f.clear();
            fx.clear();
            for (size_t i = 0; i < basis_stem.size(); ++i) {
                auto fxi = algZ::subsMod(basis_stem[i], ssCofs_[iCof].pi_f_top_cell.back(), ssCofs_[iCof].pi_gb.v_degs());
                auto fxi_v2 = ssS0_.pi_gb.ReduceV2(fxi);
                bool use_fxi_v2 = false;
                if (fxi_v2.UnknownFil() < algZ::FIL_MAX) {
                    ExtendRelS0V2(stem - d_f, fxi_v2, num_leads);
                    if (fxi_v2.UnknownFil() > algZ::FIL_MAX) {
                        x_f.push_back(basis_stem[i]);
                        fx.push_back(std::move(fxi_v2));
                        use_fxi_v2 = true;
                    }
                }
                if (!use_fxi_v2) {
                    fxi = ssS0_.pi_gb.Reduce(std::move(fxi));
                    ExtendRelS0V2(stem - d_f, fxi, num_leads);
                    if (fxi.UnknownFil() > algZ::FIL_MAX) {
                        x_f.push_back(basis_stem[i]);
                        fx.push_back(std::move(fxi));
                    }
                }
            }

            algZ::Mod1d kernel_f;
            GetKernel(x_f, fx, ssCofs_[iCof].pi_gb, ssS0_.pi_gb, kernel_f);
            ut::RemoveIf(kernel_f, [](const algZ::Mod& k) { return k.GetLead().v == 0; });

            /* image of i */
            algZ::Poly1d x_i;
            algZ::Mod1d image_i;
            int stem1 = stem;
            auto& basis_S0_stem1 = basis_S0[stem1];
            for (size_t i = 0; i < basis_S0_stem1.size(); ++i) {
                x_i.push_back(basis_S0_stem1[i]);
                auto fxi = ssCof.pi_gb.ReduceV2(algZ::Mod(basis_S0_stem1[i], 0, 0));
                image_i.push_back(std::move(fxi));
            }
            auto indices = ut::size_t_range(image_i.size());
            std::stable_sort(indices.begin(), indices.end(), [&image_i](size_t i, size_t j) { return image_i[i].UnknownFil() > image_i[j].UnknownFil(); });
            algZ::Poly1d x_i_sorted;
            algZ::Mod1d image_i_sorted;
            for (size_t i : indices) {
                x_i_sorted.push_back(std::move(x_i[i]));
                image_i_sorted.push_back(std::move(image_i[i]));
            }

            algZ::Mod1d image1_i;
            int O1 = sPossMoreEinfS0[stem];
            int O2 = O1;
            algZ::Poly gO = algZ::Mon::O(0);
            GetImage(x_i_sorted, image_i_sorted, ssS0_.pi_gb, ssCof.pi_gb, image1_i, O1, O2, gO);

            /* Check exactness */
            for (size_t i = 0; i < kernel_f.size(); ++i) {
                auto k = Residue(image1_i, std::move(kernel_f[i]), ssCof.pi_gb);
                if (k && !k.GetLead().IsUnKnown()) {
                    int s = k.GetLead().fil();
                    AdamsDeg deg(s, stem + s);
                    if (!IsPossTgt(ssCof.basis_ss, deg, kRPC)) {
                        if (s < O1) {
                            throw SSPiException(0x91866326U, "i can not hit the kernel of f");
                        }
                        if (s < O2) {
                            if (gO && !gO.GetLead().IsUnKnown()) {
                                auto& rel = algZ::Mod(gO, 0, 0) + k.LF() + algZ::Mod::O(s + 1);
                                if (depth == 0)
                                    std::cout << "Exactness f * i:  iSS=" << iCof + 1 << "  stem=" << stem << "  " << rel << "=0\n";
                                new_rels_Cofs[iCof].push_back(std::move(rel));
                                ++count_homotopy;
                            }
                        }
                    }
                    else if (s < O1) {
                        new_rels_Cofs[iCof].push_back(k);
                        if (depth == 0)
                            std::cout << "Exactness f * i:  iSS=" << iCof + 1 << "  stem=" << stem << "  " << k << "=0\n";
                        ++count_homotopy;
                    }
                }
            }
        }
    }
    AddPiRelsS0(std::move(new_rels_S0));
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
        AddPiRelsCof(iCof, std::move(new_rels_Cofs[iCof]));
        if (f_changed[iCof])
            AddPiRelsCof2S0(iCof);
    }

    return 0;
}

/*
 * Stage 0: Fast nd cache, Try{Leibniz}, maxPoss = 4
 */
void Diagram::DeduceExtensions(int& count_ss, int& count_homotopy, int depth)
{
    SyncHomotopy(count_ss, count_homotopy, depth);
    count_homotopy += DeduceTrivialExtensions(depth);
    count_homotopy += DeduceExtensionsByExactness(depth);
    std::string color, color_end = "\033[0m";

    /* top cell maps */
    {
        algZ::Poly1d new_rels_S0; /* By exactness h * fi = 0 */
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
            size_t f_size = ssCofs_[iCof].pi_f_top_cell.back().size();
            for (size_t i = 0; i < f_size; ++i) {
                int s = ssCofs_[iCof].pi_f_top_cell.back()[i].UnknownFil();
                if (s > algZ::FIL_MAX)
                    continue;
                int stem = ssCofs_[iCof].pi_gb.v_degs()[i].stem() - ssCofs_[iCof].deg_f_top_cell.stem();
                AdamsDeg deg(s, stem + s);
                if (PossMoreEinf(ssS0_.basis_ss, deg) > 0)
                    continue;
                algZ::Poly O1 = algZ::Poly::O(s + 1);
                ExtendRelS0(stem, O1);
                int ne_cout = (int)ssS0_.pi_basis[deg].size();
                int count_pass = 0;
                unsigned i_max = 1 << ne_cout;
                algZ::Poly summand_s, summand_s_pass;
                for (unsigned b = 0; b < i_max; ++b) {
                    summand_s = O1;
                    for (int j : two_expansion(b))
                        summand_s += ssS0_.pi_basis[deg][j];

                    AddNode();
                    bool bException = false;
                    try {
                        auto& f = ssCofs_[iCof].pi_f_top_cell.back()[i];
                        f.data.pop_back();
                        f += summand_s;
                        AddPiRelsCof2S0(iCof);
                        int count_ss1 = 0, count_homotopy1 = 0;
                        SyncHomotopy(count_ss1, count_homotopy1, depth + 1);
                        DeduceTrivialExtensions(depth + 1);
                        DeduceExtensionsByExactness(depth + 1);
                    }
                    catch (SSException&) {
                        bException = true;
                    }
                    PopNode();

                    if (!bException) {
                        ++count_pass;
                        if (count_pass > 1)
                            break;
                        summand_s_pass = std::move(summand_s);
                    }
                }

                if (count_pass == 0)
                    throw SSPiException(0x577b825eU, "No compatible extensions");
                else if (count_pass == 1) {
                    ++count_homotopy;
                    auto& f = ssCofs_[iCof].pi_f_top_cell.back()[i];
                    if (depth == 0)
                        std::cout << "DeduceExtension:  iSS=" << iCof + 1 << "  f" << std::to_string(i) << "=" << f << " --> ";
                    f.data.pop_back();
                    f += summand_s_pass;
                    if (depth == 0)
                        std::cout << f << '\n';

                    AddPiRelsCof2S0(iCof);
                    SyncHomotopy(count_ss, count_homotopy, depth);
                    count_homotopy += DeduceTrivialExtensions(depth);
                    count_homotopy += DeduceExtensionsByExactness(depth);

                    algZ::Poly h = ssS0_.pi_gb.Gen((uint32_t)iCof);
                    algZ::Poly relS0 = ssS0_.pi_gb.Reduce(h * f);
                    if (algZ::IsValidRel(relS0)) {
                        std::cout << "    -> S0 rel: " << relS0 << '\n';
                        new_rels_S0.push_back(std::move(relS0));
                        ++count_homotopy;
                    }
                }
            }
        }
        AddPiRelsS0(std::move(new_rels_S0));
    }
     SimplifyPiRels();

    algZ::Poly tmp;
    algZ::Mod tmpm;

    /* multiplicative structures */
    int old_count_homotopy = count_homotopy;
    std::vector<size_t> indices_start(all_basis_ss_.size(), 0);
    while (true) {
        { /* sphere */
            auto& pi_gb = ssS0_.pi_gb;
            size_t iSS = 0;

            size_t indices_end = pi_gb.data().size();
            for (size_t i = indices_start[iSS]; i < indices_end; ++i) {
                algZ::Poly rel = pi_gb.data()[i];
                int s = rel.UnknownFil();
                if (s > algZ::FIL_MAX)
                    continue;
                int stem = GetDeg(rel.GetLead(), pi_gb.gen_degs()).stem();
                AdamsDeg deg(s, stem + s);
                if (deg.t > ssS0_.t_max || PossMoreEinf(ssS0_.basis_ss, deg) > 0)
                    continue;
                algZ::Poly O1 = algZ::Poly::O(s + 1);
                ExtendRelS0(stem, O1);
                int ne_cout = (int)ssS0_.pi_basis[deg].size();
                int count_pass = 0;
                unsigned i_max = 1 << ne_cout;
                algZ::Poly summand_s, summand_s_pass, rel1;
                for (unsigned b = 0; b < i_max; ++b) {
                    summand_s = O1;
                    for (int j : two_expansion(b))
                        summand_s.iaddP(ssS0_.pi_basis[deg][j], tmp);

                    AddNode();
                    bool bException = false;
                    try {
                        rel1 = rel;
                        rel1.data.pop_back();
                        rel1.iaddP(summand_s, tmp);

                        AddPiRelsS0({std::move(rel1)});
                        int count_ss1 = 0, count_homotopy1 = 0;
                        SyncHomotopy(count_ss1, count_homotopy1, depth + 1);
                        DeduceTrivialExtensions(depth + 1);
                        DeduceExtensionsByExactness(depth + 1);
                    }
                    catch (SSException&) {
                        bException = true;
                    }
                    PopNode();

                    if (!bException) {
                        ++count_pass;
                        if (count_pass > 1)
                            break;
                        summand_s_pass = std::move(summand_s);
                    }
                }

                if (count_pass == 0)
                    throw SSPiException(0x577b825eU, "No compatible extensions");
                else if (count_pass == 1) {
                    ++count_homotopy;
                    rel1 = rel;
                    rel1.data.pop_back();
                    rel1.iaddP(summand_s_pass, tmp);
                    if (depth == 0) {
                        if (algZ::IsValidRel(summand_s_pass))
                            color = "\033[38;2;255;128;128m";
                        else
                            color = "";
                        std::cout << "DeduceExtension:  iSS=" << 0 << "  " << i << '/' << pi_gb.data().size() << "  " << rel << " --> " << rel1 << '\n';
                    }
                    AddPiRelsS0({std::move(rel1)});
                    SyncHomotopy(count_ss, count_homotopy, depth);
                    count_homotopy += DeduceTrivialExtensions(depth);
                    count_homotopy += DeduceExtensionsByExactness(depth);
                }
            }
            indices_start[iSS] = indices_end;
        }

        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) { /* module */
            auto& basis_ss = ssCofs_[iCof].basis_ss;
            auto& pi_gb = ssCofs_[iCof].pi_gb;
            auto& pi_basis = ssCofs_[iCof].pi_basis;
            size_t iSS = iCof + 1;

            size_t indices_end = pi_gb.data().size();
            for (size_t i = indices_start[iSS]; i < indices_end; ++i) {
                algZ::Mod rel = pi_gb.data()[i];
                int s = rel.UnknownFil();
                if (s > algZ::FIL_MAX)
                    continue;
                int stem = GetDeg(rel.GetLead(), ssS0_.pi_gb.gen_degs(), pi_gb.v_degs()).stem();
                AdamsDeg deg(s, stem + s);
                if (deg.t > ssCofs_[iCof].t_max || PossMoreEinf(basis_ss, deg) > 0)
                    continue;

                algZ::Mod O1 = algZ::Mod::O(s + 1);
                ExtendRelCof(iCof, stem, O1);
                int ne_cout = (int)pi_basis[deg].size();
                int count_pass = 0;
                unsigned i_max = 1 << ne_cout;
                algZ::Mod summand_s, summand_s_pass, rel1;
                for (unsigned b = 0; b < i_max; ++b) {
                    summand_s = O1;
                    for (int j : two_expansion(b))
                        summand_s.iaddP(pi_basis[deg][j], tmpm);

                    AddNode();
                    bool bException = false;
                    try {
                        rel1 = rel;
                        rel1.data.pop_back();
                        rel1.iaddP(summand_s, tmpm);

                        AddPiRelsCof(iCof, {std::move(rel1)});
                        int count_ss1 = 0, count_homotopy1 = 0;
                        SyncHomotopy(count_ss1, count_homotopy1, depth + 1);
                        DeduceTrivialExtensions(depth + 1);
                        DeduceExtensionsByExactness(depth + 1);
                    }
                    catch (SSException&) {
                        bException = true;
                    }
                    PopNode();

                    if (!bException) {
                        ++count_pass;
                        if (count_pass > 1)
                            break;
                        summand_s_pass = std::move(summand_s);
                    }
                }

                if (count_pass == 0)
                    throw SSPiException(0x577b825eU, "No compatible extensions");
                else if (count_pass == 1) {
                    ++count_homotopy;
                    rel1 = rel;
                    rel1.data.pop_back();
                    rel1.iaddP(summand_s_pass, tmpm);
                    if (depth == 0) {
                        if (algZ::IsValidRel(summand_s_pass))
                            color = "\033[38;2;255;128;128m";
                        else
                            color = "";
                        std::cout << color << "DeduceExtension:  iSS=" << iCof + 1 << "  " << i << '/' << pi_gb.data().size() << "  " << rel << " --> " << rel1 << color_end << '\n';
                    }

                    AddPiRelsCof(iCof, {std::move(rel1)});
                    SyncHomotopy(count_ss, count_homotopy, depth);
                    count_homotopy += DeduceTrivialExtensions(depth);
                    count_homotopy += DeduceExtensionsByExactness(depth);
                }
            }
            indices_start[iSS] = indices_end;
        }

        if (old_count_homotopy != count_homotopy)
            old_count_homotopy = count_homotopy;
        else
            break;
    }
}
