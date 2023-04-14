#include "main.h"
#include "mylog.h"


Diagram::Diagram(const std::vector<std::string>& dbnames, DeduceFlag flag)
{
    {
        DBSS db(dbnames[0]);
        std::string complexName = "S0";
        std::string table_S0 = "S0_AdamsE2";

        ssS0_.name = complexName;
        ssS0_.basis = db.load_basis(table_S0);
        ssS0_.degs_basis_order_by_stem = OrderDegsByStem(ssS0_.basis);
        ssS0_.t_max = ssS0_.basis.rbegin()->first.t;
        ssS0_.nodes_ss = {db.load_basis_ss(table_S0), {}};
        ssS0_.nodes_ss.reserve(MAX_NUM_NODES + 1);
        ssS0_.gb = Groebner(ssS0_.t_max, {}, db.load_gb(table_S0, DEG_MAX));

        if (flag & DeduceFlag::homotopy) {
            ssS0_.pi_gen_Einf = db.get_column_from_str<Poly>(complexName + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Poly>);
            ssS0_.pi_gb = algZ::Groebner(ssS0_.t_max, db.load_pi_gen_adamsdegs(complexName), db.load_pi_gb(complexName, DEG_MAX), true);
            ssS0_.pi_basis.reserve(MAX_NUM_NODES);
            if (flag & DeduceFlag::homotopy_def)
                db.load_pi_def(complexName, ssS0_.pi_gen_defs, ssS0_.pi_gen_def_mons);
        }

        all_names_.push_back(ssS0_.name);
        all_basis_ss_.push_back(&ssS0_.nodes_ss);
        all_t_max_.push_back(ssS0_.t_max);
    }

    for (size_t i = 1; i < dbnames.size(); ++i) {
        DBSS dbCof(dbnames[i]);
        std::string table_CW = GetE2TablePrefix(dbnames[i]);
        std::string complexName = GetComplexName(dbnames[i]);
        SSMod ssCof;
        ssCof.name = complexName;
        ssCof.basis = dbCof.load_basis_mod(table_CW);
        ssCof.degs_basis_order_by_stem = OrderDegsByStem(ssCof.basis);
        ssCof.t_max = ssCof.basis.rbegin()->first.t;
        ssCof.nodes_ss = {dbCof.load_basis_ss(table_CW), {}};
        ssCof.nodes_ss.reserve(MAX_NUM_NODES + 1);
        Mod1d xs = dbCof.load_gb_mod(table_CW, DEG_MAX);
        ssCof.gb = GroebnerMod(&ssS0_.gb, ssCof.t_max, {}, std::move(xs));
        ssCof.qt = dbCof.get_column_from_str<Poly>(table_CW + "_generators", "to_S0", "", myio::Deserialize<Poly>);
        ssCof.deg_qt = AdamsDeg(0, GetTopCellT(dbnames[i]));

        if (flag & DeduceFlag::homotopy) {
            ssCof.pi_gen_Einf = dbCof.get_column_from_str<Mod>(complexName + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Mod>);
            ssCof.pi_gb = algZ::GroebnerMod(&ssS0_.pi_gb, ssCof.t_max, dbCof.load_pi_gen_adamsdegs(complexName), dbCof.load_pi_gb_mod(complexName, DEG_MAX), true);
            ssCof.pi_basis.reserve(MAX_NUM_NODES);
            if (flag & DeduceFlag::homotopy_def)
                dbCof.load_pi_def(complexName, ssCof.pi_gen_defs, ssCof.pi_gen_def_mons);
            ssCof.pi_qt = {dbCof.get_column_from_str<algZ::Poly>(complexName + "_pi_generators", "to_S0", "ORDER BY id", myio::Deserialize<algZ::Poly>)};
            ssCof.pi_qt.reserve(MAX_NUM_NODES);
        }

        ssCofs_.push_back(std::move(ssCof));
    }
    for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
        all_names_.push_back(ssCofs_[iCof].name);
        all_basis_ss_.push_back(&ssCofs_[iCof].nodes_ss);
        all_t_max_.push_back(ssCofs_[iCof].t_max);
    }

    if (flag & DeduceFlag::homotopy)
        UpdateAllPossEinf();

    // VersionConvertReorderRels();
}

void Diagram::save(const std::vector<std::string>& dbnames, DeduceFlag flag)
{
    for (size_t k = 0; k < dbnames.size(); ++k) {
        DBSS db(dbnames[k]);
        auto pi_table = GetComplexName(dbnames[k]);
        db.begin_transaction();

        db.update_basis_ss(pi_table + "_AdamsE2", (*all_basis_ss_[k])[1]);
        if (flag & DeduceFlag::homotopy) {
            db.drop_and_create_pi_relations(pi_table);
            db.drop_and_create_pi_basis(pi_table);
            if (flag & DeduceFlag::homotopy_def)
                db.drop_and_create_pi_definitions(pi_table);
            if (k == 0) {
                db.drop_and_create_pi_generators(pi_table);
                db.save_pi_generators(pi_table, ssS0_.pi_gb.gen_degs(), ssS0_.pi_gen_Einf);
                db.save_pi_gb(pi_table, ssS0_.pi_gb.OutputForDatabase(), GetS0GbEinf());
                db.save_pi_basis(pi_table, ssS0_.pi_basis.front());
                if (flag & DeduceFlag::homotopy_def)
                    db.save_pi_def(pi_table, ssS0_.pi_gen_defs, ssS0_.pi_gen_def_mons);
            }
            else {
                db.drop_and_create_pi_generators_mod(pi_table);
                auto& Cof = ssCofs_[k - 1];
                if (Cof.pi_qt.size() != 1) {
                    Logger::LogException(0, 0x925afecU, "Not on the initial node\n");
                    throw MyException(0x925afecU, "Not on the initial node");
                }
                db.save_pi_generators_mod(pi_table, Cof.pi_gb.v_degs(), Cof.pi_gen_Einf, Cof.pi_qt.front());
                db.save_pi_gb_mod(pi_table, Cof.pi_gb.OutputForDatabase(), GetCofGbEinf(k - 1));
                db.save_pi_basis_mod(pi_table, Cof.pi_basis.front());
                if (flag & DeduceFlag::homotopy_def)
                    db.save_pi_def(pi_table, Cof.pi_gen_defs, Cof.pi_gen_def_mons);
            }
        }

        db.end_transaction();
    }
}

/* Add a node */
void Diagram::AddNode(DeduceFlag flag)
{
    for (size_t k = 0; k < all_basis_ss_.size(); ++k)
        all_basis_ss_[k]->push_back({});

    if (flag & DeduceFlag::homotopy) {
        ssS0_.pi_nodes_gen.push_back(ssS0_.pi_gb.gen_degs().size());
        ssS0_.pi_nodes_rel.push_back(ssS0_.pi_gb.data().size());
        ssS0_.pi_nodes_gen_2tor_degs.push_back(ssS0_.pi_gb.gen_2tor_degs());
        ssS0_.pi_basis.push_back({});
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
            auto& ssCof = ssCofs_[iCof];
            ssCof.pi_nodes_gen.push_back(ssCof.pi_gb.v_degs().size());
            ssCof.pi_nodes_rel.push_back(ssCof.pi_gb.data().size());
            ssCof.pi_qt.push_back(ssCof.pi_qt.back());
            ssCof.pi_basis.push_back({});
        }
    }
}

/* Pop the lastest node */
void Diagram::PopNode(DeduceFlag flag)
{
    for (size_t k = 0; k < all_basis_ss_.size(); ++k)
        all_basis_ss_[k]->pop_back();
    if (flag & DeduceFlag::homotopy) {
        ssS0_.pi_gb.Pop(ssS0_.pi_nodes_gen.back(), ssS0_.pi_nodes_rel.back());
        ssS0_.pi_gb.set_gen_2tor_degs(std::move(ssS0_.pi_nodes_gen_2tor_degs.back()));
        ssS0_.pi_gen_Einf.resize(ssS0_.pi_nodes_gen.back());
        ssS0_.pi_nodes_gen.pop_back();
        ssS0_.pi_nodes_rel.pop_back();
        ssS0_.pi_nodes_gen_2tor_degs.pop_back();
        ssS0_.pi_basis.pop_back();
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof) {
            auto& ssCof = ssCofs_[iCof];
            ssCof.pi_gb.Pop(ssCof.pi_nodes_gen.back(), ssCof.pi_nodes_rel.back());
            ssCof.pi_gen_Einf.resize(ssCof.pi_nodes_gen.back());
            ssCof.pi_qt.pop_back();
            ssCof.pi_nodes_gen.pop_back();
            ssCof.pi_nodes_rel.pop_back();
            ssCof.pi_basis.pop_back();
        }
        UpdateAllPossEinf();
    }
}



int Diagram::SetS0DiffGlobal(AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry)
{
    AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
    int result = 0;

    auto& basis = ssS0_.basis;
    auto& nodes_ss = ssS0_.nodes_ss;
    int t_max = ssS0_.t_max;

    if (x.empty()) {
        if (!dx.empty() && !IsZeroOnLevel(GetRecentStaircase(nodes_ss, deg_dx), dx, r))
            result += SetS0ImageLeibniz(deg_dx, dx, r - 1);
    }
    else if (IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        int r_min = kLevelMin;
        while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
            ++r_min;
        if (dx.empty()) {
            int r_max = NextRTgt(nodes_ss, t_max, deg_x, r + 1);
            if (r_max == -1)
                r = kRPC - 1;
            else
                r = r_max - 1;
        }
        result += SetS0DiffLeibniz(deg_x, x, dx, r, r_min, bFastTry);
    }
    return result;
}

int Diagram::SetCofDiffGlobal(size_t iCof, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry)
{
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    int count = 0;

    auto& basis = ssCofs_[iCof].basis;
    auto& nodes_ss = ssCofs_[iCof].nodes_ss;
    int t_max = ssCofs_[iCof].t_max;

    if (x.empty()) {
        if (!dx.empty() && !IsZeroOnLevel(GetRecentStaircase(nodes_ss, deg_dx), dx, r)) {
            count += SetCofImageLeibniz(iCof, deg_dx, dx, r - 1);

            Mod mod_dx = Indices2Mod(dx, basis.at(deg_dx));
            Poly poly_dx_S0 = ssS0_.gb.Reduce(subsMod(mod_dx, ssCofs_[iCof].qt));
            if (poly_dx_S0) {
                AdamsDeg deg_dx_S0 = deg_dx - ssCofs_[iCof].deg_qt;
                int1d dx_S0 = Poly2Indices(poly_dx_S0, ssS0_.basis.at(deg_dx_S0));
                count += SetS0ImageLeibniz(deg_dx_S0, dx_S0, r - 1);
            }
        }
    }
    else if (IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        int r_min = kLevelMin;
        while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
            ++r_min;
        if (dx.empty()) {
            int r_max = NextRTgt(nodes_ss, t_max, deg_x, r + 1);
            if (r_max == -1)
                r = kRPC - 1;
            else
                r = r_max - 1;
        }
        count += SetCofDiffLeibniz(iCof, deg_x, x, dx, r, r_min, bFastTry);

        Mod mod_x = Indices2Mod(x, basis.at(deg_x));
        Poly poly_x_S0 = ssS0_.gb.Reduce(subsMod(mod_x, ssCofs_[iCof].qt));
        Mod mod_dx = !dx.empty() ? Indices2Mod(dx, basis.at(deg_dx)) : Mod();
        Poly poly_dx_S0 = ssS0_.gb.Reduce(subsMod(mod_dx, ssCofs_[iCof].qt));
        if (poly_x_S0 || poly_dx_S0) {
            AdamsDeg deg_x_S0 = deg_x - ssCofs_[iCof].deg_qt;
            int1d x_S0 = poly_x_S0 ? Poly2Indices(poly_x_S0, ssS0_.basis.at(deg_x_S0)) : int1d{};
            AdamsDeg deg_dx_S0 = deg_dx - ssCofs_[iCof].deg_qt;
            int1d dx_S0 = poly_dx_S0 ? Poly2Indices(poly_dx_S0, ssS0_.basis.at(deg_dx_S0)) : int1d{};
            count += SetS0DiffGlobal(deg_x_S0, x_S0, dx_S0, r, bFastTry);
        }
    }
    return count;
}

int Diagram::SetDiffGlobal(size_t iSS, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry)
{
    if (iSS == 0)
        return SetS0DiffGlobal(deg_x, x, dx, r, bFastTry);
    else
        return SetCofDiffGlobal(iSS - 1, deg_x, x, dx, r, bFastTry);
}