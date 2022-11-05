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

bool Diagram::PossEinf(const Staircases1d& basis_ss, AdamsDeg deg) const
{
    if (basis_ss.front().find(deg) != basis_ss.front().end()) {
        const Staircase& sc = GetRecentStaircase(basis_ss, deg);
        size_t i_start_perm = GetFirstIndexOnLevel(sc, kLevelMax / 2);
        size_t i_stable = GetFirstIndexOfFixedLevels(basis_ss, deg, kLevelPC + 1);
        return i_start_perm < i_stable;
    }
    return false;
}

int Diagram::NextSExt(const Staircases1d& basis_ss, int t_max, int stem, int s_min) const
{
    int s_max = (stem + 3) / 2;
    for (int s = s_min; s <= s_max; ++s) {
        if (stem + s > t_max || PossEinf(basis_ss, AdamsDeg(s, stem + s)))
            return s;
    }
    return algZ::FIL_MAX + 1;
}

int2d Diagram::GetS0GbEinf(AdamsDeg deg) const
{
    auto& gb = ssS0_.gb;
    auto& pi_gb = ssS0_.pi_gb;
    auto& pi_gen_Einf = ssS0_.pi_gen_Einf;

    int2d result;
    Poly tmp;
    algZ::Poly1d gb_Einf = pi_gb.RelsLF(deg);
    Poly1d new_boundaries;
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
    Poly1d new_boundaries;
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

void Diagram::SyncS0Homotopy(int& count_ss, int& count_homotopy)
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
                    if (count1 > 0)
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

void Diagram::SyncCofHomotopy(int iCof, int& count_ss, int& count_homotopy)
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
                    if (count1 > 0)
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
                        pi_f_top_cell.push_back(std::move(pi_fx));
                    }
                    else
                        pi_f_top_cell.push_back(algZ::Mon::O(deg_S0.s + 1));
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

int Diagram::DeduceZeroExtensions()
{
    int count_homotopy = 0;

    /* homomorphisms */
    algZ::Poly1d new_rels_S0_exact; /* By exactness fi * h = 0 */
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
        auto& pi_f = ssCofs_[iCof].pi_f_top_cell;
        for (size_t i = 0; i < pi_f.size(); ++i) {
            auto& fi = pi_f[i];
            if (fi) {
                auto& m = fi.data.back();
                if (m.IsUnKnown()) {
                    AdamsDeg deg = ssCofs_[iCof].pi_gb.v_degs()[i] - ssCofs_[iCof].deg_f_top_cell;
                    int s = NextSExt(ssS0_.basis_ss, ssS0_.t_max, deg.stem(), m.fil());
                    if (s > m.fil()) {
                        std::cout << "For degree reason:  iSS=" << iCof + 1 << "  f" << std::to_string(i) << "=" << fi;
                        if (s <= algZ::FIL_MAX)
                            fi.data.back() = algZ::Mon::O(s);
                        else
                            fi.data.pop_back();
                        std::cout << " --> " << fi << '\n';
                        ++count_homotopy;

                        algZ::Poly h = ssS0_.pi_gb.Gen((uint32_t)iCof);
                        algZ::Poly relS0 = ssS0_.pi_gb.Reduce(fi * h);
                        if (algZ::IsValidRel(relS0)) {
                            std::cout << "    -> S0 rel: " << relS0 << '\n';
                            new_rels_S0_exact.push_back(std::move(relS0));
                            ++count_homotopy;
                        }
                    }
                }
            }
        }
    }
    ssS0_.pi_gb.AddRels(new_rels_S0_exact, ssS0_.t_max);
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
        ssCofs_[iCof].pi_gb.AddRels({}, ssCofs_[iCof].t_max);

    /* module structures */
    int old_count_homotopy = count_homotopy;
    while (true) {
        for (size_t iSS = 0; iSS < all_basis_ss_.size(); ++iSS) {
            auto& basis_ss = *all_basis_ss_[iSS];
            int t_max = all_t_max_[iSS];

            if (iSS == 0) { /* sphere */
                algZ::Poly1d new_rels;
                for (auto& rel : ssS0_.pi_gb.data()) {
                    auto& m = rel.data.back();
                    if (m.IsUnKnown()) {
                        AdamsDeg deg = GetDeg(rel.GetLead(), ssS0_.pi_gb.gen_degs());
                        int s = NextSExt(basis_ss, t_max, deg.stem(), m.fil());
                        if (s > m.fil()) {
                            algZ::Poly rel1 = rel;
                            if (s <= algZ::FIL_MAX)
                                rel1.data.back() = algZ::Mon::O(s);
                            else
                                rel1.data.pop_back();
                            rel1 = ssS0_.pi_gb.Reduce(std::move(rel1));
                            if (algZ::IsValidRel(rel1)) {
                                std::cout << "For degree reason:  iSS=" << iSS << "  " << rel << " --> " << rel1 << '\n';
                                new_rels.push_back(std::move(rel1));
                                ++count_homotopy;
                            }
                        }
                    }
                }
                ssS0_.pi_gb.AddRels(new_rels, t_max);
                for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
                    ssCofs_[iCof].pi_gb.AddRels({}, ssCofs_[iCof].t_max);
            }
            else { /* module */
                int iCof = int(iSS - 1);

                algZ::Mod1d new_rels;
                algZ::Poly1d new_rels_S0;
                for (auto& rel : ssCofs_[iCof].pi_gb.data()) {
                    auto& m = rel.data.back();
                    if (m.IsUnKnown()) {
                        AdamsDeg deg = GetDeg(rel.GetLead(), ssS0_.pi_gb.gen_degs(), ssCofs_[iCof].pi_gb.v_degs());
                        int s = NextSExt(basis_ss, t_max, deg.stem(), m.fil());
                        if (s > m.fil()) {
                            algZ::Mod rel1 = rel;
                            if (s <= algZ::FIL_MAX)
                                rel1.data.back() = algZ::MMod::O(s);
                            else
                                rel1.data.pop_back();
                            rel1 = ssCofs_[iCof].pi_gb.Reduce(std::move(rel1));
                            algZ::Poly relS0 = algZ::subsMod(rel1, ssCofs_[iCof].pi_f_top_cell, ssCofs_[iCof].pi_gb.v_degs());
                            relS0 = ssS0_.pi_gb.Reduce(std::move(relS0));
                            if (algZ::IsValidRel(rel1)) {
                                std::cout << "For degree reason:  iSS=" << iSS << "  " << rel << " --> " << rel1 << '\n';
                                if (algZ::IsValidRel(relS0)) {
                                    std::cout << "    -> S0 rel: " << relS0 << '\n';
                                    new_rels_S0.push_back(std::move(relS0));
                                    ++count_homotopy;
                                }
                                new_rels.push_back(std::move(rel1));
                                ++count_homotopy;
                            }
                        }
                    }
                }
                ssS0_.pi_gb.AddRels(new_rels_S0, ssS0_.t_max);
                ssCofs_[iCof].pi_gb.AddRels(new_rels, t_max);
            }
        }
        if (old_count_homotopy != count_homotopy)
            old_count_homotopy = count_homotopy;
        else
            break;
    }

    return count_homotopy;
}

int Diagram::DeduceExtensionsByExactness()
{
    int count_homotopy = 0;

    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
        auto& pi_f = ssCofs_[iCof].pi_f_top_cell;
        
         
    }
    return 0;
}
