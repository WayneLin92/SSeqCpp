#include "mylog.h"
#include "main.h"
#include "algebras/linalg.h"

bool EndWithIndAndO(const algZ::Poly& p, uint32_t& gen_id)
{
    if (p.data.size() >= 2 && p.data.back().IsUnKnown()) {
        auto& m = p.data[p.data.size() - 2];
        if (m.c() == 0 && m.m().size() == 1) {
            if (m.m().begin()->e_masked() == 1) {
                gen_id = m.m().begin()->g();
                return true;
            }
        }
    }
    return false;
}

bool EndWithIndAndO(const algZ::Mod& p, uint32_t& v_id)
{
    if (p.data.size() >= 2 && p.data.back().IsUnKnown()) {
        auto& m = p.data[p.data.size() - 2];
        if (m.c() == 0 && !m.m) {
            v_id = m.v;
            return true;
        }
    }
    return false;
}

algZ::Poly DefRel(algZ::Poly p, const algZ::Groebner& gb)
{
    algZ::Poly tmp;
    p.data.pop_back();
    algZ::Poly gen = p.data.back();
    p.data.pop_back();
    p.isubP(gen, tmp, gb.gen_2tor_degs());
    p = gb.Reduce(std::move(p));
    return p;
}

algZ::Mod DefRel(algZ::Mod p, const algZ::GroebnerMod& gb)
{
    algZ::Mod tmp;
    p.data.pop_back();
    algZ::Mod gen = p.data.back();
    p.data.pop_back();
    p.isubP(gen, tmp, gb.gen_2tor_degs());
    p = gb.Reduce(std::move(p));
    return p;
}

algZ::Poly DefPoly(algZ::Poly p, const algZ::Groebner& gb)
{
    algZ::Poly tmp;
    p.data.pop_back();
    p = gb.Reduce(std::move(p));
    return p;
}

/* Get a Z2-basis of the kernel.
 * fx should all be certain.
 */
template <typename T1, typename T2, typename GB1, typename GB2>
void GetKernel(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2, const int O, std::vector<T1>& kernel)
{
    /* Sort by certainty */
    auto indices = ut::size_t_range(x.size());
    std::stable_sort(indices.begin(), indices.end(), [&x](size_t i, size_t j) { return x[i].UnknownFil() > x[j].UnknownFil(); });

    std::vector<T2> image;
    std::vector<T1> g;
    T1 tmp1{};
    T2 tmp2{};
    /* f(g[i]) = image[i] */
    for (size_t i : indices) {
        T1 src = x[i];
        T2 tgt = fx[i];
        if (O <= algZ::FIL_MAX) {
            tgt.iaddP(T2::O(O), tmp2);
            tgt.data.pop_back();
        }
        {
            size_t index = 0;
            while (index < tgt.data.size()) {
                size_t d_index = 1;
                for (size_t j = 0; j < image.size(); j++) {
                    if (tgt.data[index] == image[j].GetLead()) {
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
            while (index < src.data.size() && !src.data[index].IsUnKnown()) {
                size_t d_index = 1;
                for (size_t j = 0; j < kernel.size(); j++) {
                    if (src.data[index] == kernel[j].GetLead()) {
                        int c = src.data[index].c() - kernel[j].GetLead().c();
                        src.isubmulP(algZ::Mon::twoTo(c), kernel[j], tmp1, gb1.gen_2tor_degs());
                        src = gb1.Reduce(std::move(src));
                        d_index = 0;
                        break;
                    }
                }
                index += d_index;
            }
            if (algZ::IsValidRel(src))
                kernel.push_back(std::move(src));
        }
    }
}

/* Group the leads by s and return the numbers of leads for each s */
template <typename T1, typename T2, typename GB1, typename GB2>
void GetImageLeads(const std::vector<T1>& x, const std::vector<T2>& fx, const GB1& gb1, const GB2& gb2, ut::map_seq<int, 0>& num_leads)
{
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
                ++num_leads[tgt.GetLead().fil()];
                image.push_back(std::move(tgt));
                g.push_back(std::move(src));
            }
        }
    }
}

int Diagram::DefineDependenceInExtensions(int depth)
{
    int count_homotopy = 0;

    /* multiplicative structures */
    { /* sphere */
        auto& pi_gb = ssS0_.pi_gb;
        if (ssS0_.pi_gen_defs.size() < ssS0_.pi_gb.gen_degs().size())
            ssS0_.pi_gen_defs.resize(ssS0_.pi_gb.gen_degs().size(), DefFlag::no_def);
        if (ssS0_.pi_gen_def_mons.size() < ssS0_.pi_gb.gen_degs().size())
            ssS0_.pi_gen_def_mons.resize(ssS0_.pi_gb.gen_degs().size());

        algZ::Poly1d new_rels;
        for (size_t i = 0; i < pi_gb.data().size(); ++i) {
            auto& rel = pi_gb.data()[i];
            uint32_t gen_id = -1;
            if (EndWithIndAndO(rel, gen_id)) {
                AdamsDeg deg = GetDeg(rel.GetLead(), pi_gb.gen_degs());
                if (ssS0_.pi_gen_defs[gen_id] == DefFlag::no_def) {
                    ssS0_.pi_gen_defs[gen_id] = DefFlag::dec;
                    algZ::Poly rel1 = DefRel(rel, pi_gb);

                    Logger::LogHtpyRel(depth, enumReason::def, "S0", deg, rel, rel1);
                    new_rels.push_back(std::move(rel1));
                    ++count_homotopy;
                }
            }
        }
        AddPiRelsS0(std::move(new_rels));
    }
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) { /* module */
        auto& ssCof = ssCofs_[iCof];
        auto& pi_gb = ssCof.pi_gb;
        auto& nodes_ss = ssCof.nodes_ss;

        if (ssCof.pi_gen_defs.size() < ssCof.pi_gb.v_degs().size())
            ssCof.pi_gen_defs.resize(ssCof.pi_gb.v_degs().size(), DefFlag::no_def);
        if (ssCof.pi_gen_def_mons.size() < ssCof.pi_gb.v_degs().size())
            ssCof.pi_gen_def_mons.resize(ssCof.pi_gb.v_degs().size());

        algZ::Mod1d new_rels;
        for (size_t i = 0; i < pi_gb.data().size(); ++i) {
            auto& rel = pi_gb.data()[i];
            uint32_t v_id = -1;
            if (EndWithIndAndO(rel, v_id)) {
                AdamsDeg deg = GetDeg(rel.GetLead(), ssS0_.pi_gb.gen_degs(), pi_gb.v_degs());
                if (ssCof.pi_gen_defs[v_id] == DefFlag::no_def) {
                    ssCof.pi_gen_defs[v_id] = DefFlag::dec;
                    algZ::Mod rel1 = DefRel(rel, pi_gb);
                    Logger::LogHtpyRel(depth, enumReason::def, ssCof.name, deg, rel, rel1);
                    new_rels.push_back(std::move(rel1));
                    ++count_homotopy;
                }
            }
        }
        AddPiRelsCof(iCof, std::move(new_rels));
    }

    /* top cell maps */
    //{
    //    algZ::Poly1d new_rels_S0; /* By exactness qi * h = 0 */
    //    int old_count_homotopy = count_homotopy;
    //    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
    //        auto& ssCof = ssCofs_[iCof];
    //        auto& q = ssCof.pi_qt.back();
    //        for (size_t i = 0; i < q.size(); ++i) {
    //            auto& qi = q[i];
    //            qi = ssS0_.pi_gb.Reduce(std::move(qi));
    //            uint32_t gen_id = -1;
    //            if (EndWithIndAndO(qi, gen_id)) {
    //                AdamsDeg deg = GetDeg(qi.GetLead(), ssS0_.pi_gb.gen_degs());
    //                if (ssS0_.pi_gen_defs[gen_id] == DefFlag::no_def) {
    //                    ssS0_.pi_gen_defs[gen_id] = DefFlag::dec;
    //                    algZ::Poly qi1 = DefPoly(qi, ssS0_.pi_gb);

    //                    if (depth == 0)
    //                        myio::Logger::cout_fout2() << "By definition:  " << ssCof.name << "  stem=" << deg.stem() << "  q" << std::to_string(i) << "=" << qi << " --> " << qi1 << '\n';
    //                    qi = std::move(qi1);
    //                    ++count_homotopy;

    //                    algZ::Poly h = ssS0_.pi_gb.Gen((uint32_t)iCof);
    //                    algZ::Poly relS0 = ssS0_.pi_gb.Reduce(qi * h);
    //                    if (algZ::IsValidRel(relS0)) {
    //                        if (depth == 0)
    //                            myio::Logger::cout_fout2() << "    -> S0 rel: " << relS0 << '\n';
    //                        new_rels_S0.push_back(std::move(relS0));
    //                        ++count_homotopy;
    //                    }
    //                }
    //            }
    //        }
    //        if (count_homotopy > old_count_homotopy) {
    //            AddPiRelsCof2S0(iCof);
    //            old_count_homotopy = count_homotopy;
    //        }
    //    }
    //    AddPiRelsS0(std::move(new_rels_S0));
    //}

    return count_homotopy;
}

void ReduceIndeterminancyId(algZ::Poly1d& indeterminancy, const GenConstraint& constraint, const algZ::Groebner& gb)
{
    algZ::Poly1d m_by_ind;
    for (auto& b : indeterminancy)
        m_by_ind.push_back(gb.Reduce(algZ::Poly(constraint.m) * b));
    algZ::Poly1d kernel;
    GetKernel(indeterminancy, m_by_ind, gb, gb, constraint.O, kernel);
    indeterminancy = std::move(kernel);
}

void ReduceIndeterminancy(algZ::Poly1d& indeterminancy, const std::vector<GenConstraint>& constraints, const algZ::Groebner& gb)
{
    for (auto& gc : constraints)
        ReduceIndeterminancyId(indeterminancy, gc, gb);
}

void ReduceIndeterminancyId(algZ::Mod1d& indeterminancy, const GenConstraint& constraint, const algZ::Groebner& gb, const algZ::GroebnerMod& gbm)
{
    algZ::Mod1d m_by_ind;
    for (auto& b : indeterminancy)
        m_by_ind.push_back(gbm.Reduce(algZ::Poly(constraint.m) * b));
    algZ::Mod1d kernel;
    GetKernel(indeterminancy, m_by_ind, gbm, gbm, constraint.O, kernel);
    indeterminancy = std::move(kernel);
}

void ReduceIndeterminancyQ(algZ::Mod1d& indeterminancy, const GenConstraint& constraint, const algZ::Groebner& gb, const SSMod& ssCof)
{
    algZ::Poly1d q_m_by_ind;
    for (auto& b : indeterminancy) {
        auto q = algZ::subsMod(ssCof.pi_gb.Reduce(algZ::Poly(constraint.m) * b), ssCof.pi_qt.back(), ssCof.pi_gb.v_degs());
        q = gb.Reduce(q);
        q_m_by_ind.push_back(std::move(q));
    }
    algZ::Mod1d kernel;
    GetKernel(indeterminancy, q_m_by_ind, ssCof.pi_gb, gb, constraint.O, kernel);
    indeterminancy = std::move(kernel);
}

void ReduceIndeterminancy(algZ::Mod1d& indeterminancy, const std::vector<GenConstraint>& constraints, const algZ::Groebner& gb, const SSMod& ssCof)
{
    for (auto& gc : constraints) {
        if (gc.map_id == 0)
            ReduceIndeterminancyId(indeterminancy, gc, gb, ssCof.pi_gb);
        else if (gc.map_id == 1)
            ReduceIndeterminancyQ(indeterminancy, gc, gb, ssCof);
        else
            throw MyException(0xf095ff26U, "Unknown map_id");
    }
}

int Diagram::DefineDependenceInExtensionsV2(int stem_max_mult, int depth)
{
    int count_homotopy = 0;

    algZ::Poly tmp;
    algZ::Mod tmpm;

    /* Populate `basis_S0` and `basis_Cofs` */
    algZ::Mon2d basis_S0(size_t(ssS0_.t_max + 1));
    std::vector<algZ::MMod2d> basis_Cofs(ssCofs_.size());
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
        basis_Cofs[iCof].resize(size_t(ssCofs_[iCof].t_max + 1));
    if (depth == 0) {
        for (auto& [deg, pi_basis_d] : ssS0_.pi_basis.front())
            for (auto& m : pi_basis_d.pi_basis)
                basis_S0[deg.stem()].push_back(m);
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
            for (auto& [deg, pi_basis_d] : ssCofs_[iCof].pi_basis.front())
                for (auto& m : pi_basis_d.pi_basis)
                    basis_Cofs[iCof][deg.stem()].push_back(m);
    }
    else {
        std::set<AdamsDeg> degs_S0;
        for (auto& pi_basis_i : ssS0_.pi_basis)
            for (auto& [deg, _] : pi_basis_i)
                degs_S0.insert(deg);
        for (AdamsDeg deg : degs_S0)
            for (auto& m : GetRecentPiBasis(ssS0_.pi_basis, deg)->pi_basis)
                basis_S0[deg.stem()].push_back(m);
        std::set<AdamsDeg> degs_Cof;
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
            degs_Cof.clear();
            for (auto& pi_basis_i : ssCofs_[iCof].pi_basis)
                for (auto& [deg, _] : pi_basis_i)
                    degs_Cof.insert(deg);
            for (AdamsDeg deg : degs_Cof)
                for (auto& m : GetRecentPiBasis(ssCofs_[iCof].pi_basis, deg)->pi_basis)
                    basis_Cofs[iCof][deg.stem()].push_back(m);
        }
    }

    algZ::Mon1d multipliers;
    auto degs_mult = OrderDegsByStem(ssS0_.pi_basis.front());
    for (auto& deg : degs_mult) {
        auto& bs = ssS0_.pi_basis.front().at(deg).pi_basis;
        if (deg.stem() <= stem_max_mult && deg.s <= 15)
            multipliers.insert(multipliers.end(), bs.begin(), bs.end());
    }

    for (auto& m : multipliers) {
        /* multiplicative structures */
        { /* sphere */
            auto& pi_gb = ssS0_.pi_gb;
            if (ssS0_.pi_gen_defs.size() < ssS0_.pi_gb.gen_degs().size())
                ssS0_.pi_gen_defs.resize(ssS0_.pi_gb.gen_degs().size(), DefFlag::no_def);
            if (ssS0_.pi_gen_def_mons.size() < ssS0_.pi_gb.gen_degs().size())
                ssS0_.pi_gen_def_mons.resize(ssS0_.pi_gb.gen_degs().size());

            algZ::Poly1d new_rels;
            for (size_t gen_id = 0; gen_id < pi_gb.gen_degs().size(); ++gen_id) {
                if (ssS0_.pi_gen_defs[gen_id] == DefFlag::no_def || ssS0_.pi_gen_defs[gen_id] == DefFlag::constraints) {
                    AdamsDeg d_g = ssS0_.pi_gb.gen_degs()[gen_id];
                    algZ::Poly g = ssS0_.pi_gb.ReduceV2(ssS0_.pi_gb.Gen((uint32_t)gen_id));
                    if (!g)
                        continue;
                    AdamsDeg d_m = GetDeg(m, ssS0_.pi_gb.gen_degs());
                    AdamsDeg d_prod = d_g + d_m;
                    algZ::Poly prod;
                    prod.iaddmulP(m, g, tmp, pi_gb.gen_2tor_degs());
                    auto prod_reduced = pi_gb.ReduceV2(prod);
                    if (prod_reduced.UnknownFil() > algZ::FIL_MAX)
                        continue;
                    ExtendRelS0(d_prod.stem(), prod_reduced);
                    if (d_prod.t + prod_reduced.UnknownFil() - d_prod.s > ssS0_.t_max)
                        continue;

                    algZ::Poly1d indeterminancy, m_by_ind;
                    auto& basis_stem = basis_S0[d_g.stem()];
                    for (auto& b : basis_stem)
                        if (b.fil() > d_g.s)
                            indeterminancy.push_back(b);
                    // if (gen_id == 38)
                    //     std::cout << "debug\n";
                    if (ssS0_.pi_gen_defs[gen_id] == DefFlag::constraints)
                        ReduceIndeterminancy(indeterminancy, ssS0_.pi_gen_def_mons[gen_id], ssS0_.pi_gb);
                    for (auto& b : indeterminancy)
                        m_by_ind.push_back(ssS0_.pi_gb.ReduceV2(algZ::Poly(m) * b));
                    ut::map_seq<int, 0> num_leads;
                    GetImageLeads(indeterminancy, m_by_ind, ssS0_.pi_gb, ssS0_.pi_gb, num_leads);
                    algZ::Poly prod_extended = prod_reduced;
                    if (ExtendRelS0V2(d_prod.stem(), prod_extended, num_leads)) {
                        Logger::LogHtpyProd(depth, enumReason::def, "S0", d_prod, m, g, prod_reduced, prod_extended);
                        prod_extended.isubP(prod, tmp, pi_gb.gen_2tor_degs());
                        ssS0_.pi_gen_defs[gen_id] = DefFlag::constraints;
                        ssS0_.pi_gen_def_mons[gen_id].push_back(GenConstraint{0, m, prod_extended.UnknownFil()});
                        new_rels.push_back(std::move(prod_extended));
                    }
                }
            }
            AddPiRelsS0(std::move(new_rels));
        }

        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) { /* module */
            auto& ssCof = ssCofs_[iCof];
            auto& pi_gb = ssCof.pi_gb;

            if (ssCof.pi_gen_defs.size() < ssCof.pi_gb.v_degs().size())
                ssCof.pi_gen_defs.resize(ssCof.pi_gb.v_degs().size(), DefFlag::no_def);
            if (ssCof.pi_gen_def_mons.size() < ssCof.pi_gb.v_degs().size())
                ssCof.pi_gen_def_mons.resize(ssCof.pi_gb.v_degs().size());

            algZ::Mod1d new_rels;
            for (size_t v_id = 0; v_id < pi_gb.v_degs().size(); ++v_id) {
                if (ssCof.pi_gen_defs[v_id] == DefFlag::no_def || ssCof.pi_gen_defs[v_id] == DefFlag::constraints) {
                    AdamsDeg d_g = ssCof.pi_gb.v_degs()[v_id];
                    algZ::Mod g = ssCof.pi_gb.ReduceV2(ssCof.pi_gb.Gen((uint32_t)v_id));
                    if (!g)
                        continue;
                    AdamsDeg d_m = GetDeg(m, ssS0_.pi_gb.gen_degs());
                    AdamsDeg d_prod = d_g + d_m;
                    algZ::Mod prod;
                    prod.iaddmulP(m, g, tmpm, pi_gb.gen_2tor_degs());
                    auto prod_reduced = pi_gb.ReduceV2(prod);
                    if (prod_reduced.UnknownFil() > algZ::FIL_MAX)
                        continue;
                    ExtendRelCof(iCof, d_prod.stem(), prod_reduced);
                    if (d_prod.t + prod_reduced.UnknownFil() - d_prod.s > ssCof.t_max)
                        continue;

                    /*if (iCof == 0 && v_id == 21)
                        std::cout << "debug\n";*/

                    algZ::Mod1d indeterminancy, m_by_ind;
                    auto& basis_stem = basis_Cofs[iCof][d_g.stem()];
                    for (auto& b : basis_stem)
                        if (b.fil() > d_g.s)
                            indeterminancy.push_back(b);
                    if (ssCof.pi_gen_defs[v_id] == DefFlag::constraints)
                        ReduceIndeterminancy(indeterminancy, ssCof.pi_gen_def_mons[v_id], ssS0_.pi_gb, ssCof);
                    for (auto& b : indeterminancy)
                        m_by_ind.push_back(ssCof.pi_gb.ReduceV2(algZ::Poly(m) * b));

                    ut::map_seq<int, 0> num_leads;
                    GetImageLeads(indeterminancy, m_by_ind, ssCof.pi_gb, ssCof.pi_gb, num_leads);
                    algZ::Mod prod_extended = prod_reduced;
                    if (ExtendRelCofV2(iCof, d_prod.stem(), prod_extended, num_leads)) {
                        Logger::LogHtpyProd(depth, enumReason::def, ssCof.name, d_prod, m, g, prod_reduced, prod_extended);
                        prod_extended.isubP(prod, tmpm, pi_gb.gen_2tor_degs());
                        ssCof.pi_gen_defs[v_id] = DefFlag::constraints;
                        ssCof.pi_gen_def_mons[v_id].push_back(GenConstraint{0, m, prod_extended.UnknownFil()});
                        new_rels.push_back(std::move(prod_extended));
                    }
                }
            }
            AddPiRelsCof(iCof, std::move(new_rels));
        }

        /* top cell maps */
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) { /* module */
            auto& ssCof = ssCofs_[iCof];
            auto& pi_gb = ssCof.pi_gb;
            int1d v_ids_changed;

            algZ::Poly1d new_rels_S0;
            for (size_t v_id = 0; v_id < pi_gb.v_degs().size(); ++v_id) {
                if (ssCof.pi_gen_defs[v_id] == DefFlag::no_def || ssCof.pi_gen_defs[v_id] == DefFlag::constraints) {
                    AdamsDeg d_g = ssCof.pi_gb.v_degs()[v_id];
                    algZ::Mod g = ssCof.pi_gb.ReduceV2(ssCof.pi_gb.Gen((uint32_t)v_id));
                    if (!g)
                        continue;
                    AdamsDeg d_m = GetDeg(m, ssS0_.pi_gb.gen_degs());
                    AdamsDeg d_q = d_g + d_m - ssCof.deg_qt;
                    algZ::Mod prod;
                    prod.iaddmulP(m, g, tmpm, ssCof.pi_gb.gen_2tor_degs());
                    auto q_prod = algZ::subsMod(prod, ssCof.pi_qt.back(), ssCof.pi_gb.v_degs());
                    q_prod = ssS0_.pi_gb.Reduce(q_prod);
                    if (!algZ::IsValidRel(q_prod))
                        continue;
                    auto q_prod_reduced = ssS0_.pi_gb.ReduceV2(q_prod);
                    if (q_prod_reduced.UnknownFil() > algZ::FIL_MAX)
                        continue;
                    ExtendRelS0(d_q.stem(), q_prod_reduced);
                    if (d_q.t + q_prod_reduced.UnknownFil() - d_q.s > ssS0_.t_max)
                        continue;

                    algZ::Mod1d indeterminancy;
                    algZ::Poly1d q_m_by_ind;
                    auto& basis_stem = basis_Cofs[iCof][d_g.stem()];
                    for (auto& b : basis_stem)
                        if (b.fil() > d_g.s)
                            indeterminancy.push_back(b);
                    if (ssCof.pi_gen_defs[v_id] == DefFlag::constraints)
                        ReduceIndeterminancy(indeterminancy, ssCof.pi_gen_def_mons[v_id], ssS0_.pi_gb, ssCof);
                    for (auto& b : indeterminancy) {
                        auto q = algZ::subsMod(ssCof.pi_gb.Reduce(algZ::Poly(m) * b), ssCof.pi_qt.back(), ssCof.pi_gb.v_degs());
                        q = ssS0_.pi_gb.Reduce(std::move(q));
                        q_m_by_ind.push_back(std::move(q));
                    }

                    ut::map_seq<int, 0> num_leads;
                    GetImageLeads(indeterminancy, q_m_by_ind, ssCof.pi_gb, ssS0_.pi_gb, num_leads);
                    algZ::Poly q_prod_extended = q_prod_reduced;
                    if (ExtendRelS0V2(d_q.stem(), q_prod_extended, num_leads)) {
                        if (m) {
                            std::string qm = fmt::format("q*{}", algZ::Poly(m));
                            Logger::LogHtpyMap(depth, enumReason::def, ssCof.name, d_q, qm, v_id, q_prod_reduced, q_prod_extended);
                            q_prod_extended.isubP(q_prod, tmp, pi_gb.gen_2tor_degs());
                            ssCof.pi_gen_defs[v_id] = DefFlag::constraints;
                            ssCof.pi_gen_def_mons[v_id].push_back(GenConstraint{1, m, q_prod_extended.UnknownFil()});
                            new_rels_S0.push_back(std::move(q_prod_extended));
                        }
                        else {
                            Logger::LogHtpyMap(depth, enumReason::def, ssCof.name, d_q, "q", v_id, q_prod_reduced, q_prod_extended);
                            ssCof.pi_gen_defs[v_id] = DefFlag::constraints;
                            ssCof.pi_gen_def_mons[v_id].push_back(GenConstraint{1, m, q_prod_extended.UnknownFil()});
                            ssCof.pi_qt.back()[v_id] = std::move(q_prod_extended);
                            v_ids_changed.push_back((int)v_id);
                        }
                    }
                }
            }
            AddPiRelsS0(std::move(new_rels_S0));
            if (!v_ids_changed.empty())
                AddPiRelsCof2S0(iCof);
        }
    }

    return count_homotopy;
}

int main_deduce_ext_def(int argc, char** argv, int index)
{
    std::string selector = "default";
    DeduceFlag flag = DeduceFlag::homotopy | DeduceFlag::homotopy_def;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss deduce ext_def [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);
    Diagram diagram(dbnames, flag);

    try {
        try {
            diagram.DeduceTrivialExtensions(0);
            diagram.DefineDependenceInExtensions(0);
        }
        catch (TerminationException&) {
        }
        diagram.SimplifyPiRels();
        diagram.save(dbnames, flag);
    }
#ifdef MYDEPLOY
    catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }
#endif
    catch (NoException&) {
        ;
    }

    // bench::Counter::print();
    return 0;
}

int main_deduce_ext_def2(int argc, char** argv, int index)
{
    int stem_max_mult = 9;
    std::string selector = "default";
    DeduceFlag flag = DeduceFlag::homotopy | DeduceFlag::homotopy_def;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss deduce ext_def2 [stem_max_mult] [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  stem_max_mult = " << stem_max_mult << "\n\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "stem_max_mult", stem_max_mult))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);
    Diagram diagram(dbnames, flag);

    try {
        try {
            int count_ss = 0, count_homotopy = 0;
            diagram.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_homotopy, 0);
            diagram.DefineDependenceInExtensionsV2(stem_max_mult, 0);
        }
        catch (TerminationException&) {
        }
        diagram.SimplifyPiRels();
        diagram.save(dbnames, flag);
    }
#ifdef MYDEPLOY
    catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }
#endif
    catch (NoException&) {
        ;
    }

    // bench::Counter::print();
    return 0;
}