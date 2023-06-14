#include "algebras/linalg.h"
#include "main.h"

class MyDB : public DBSS
{
    using Statement = myio::Statement;

public:
    MyDB() = default;
    explicit MyDB(const std::string& filename) : DBSS(filename) {}

    void create_basis_products(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis_products (id1 INTEGER, id2 INTEGER, prod TEXT, PRIMARY KEY(id1, id2));");
    }

    void drop_and_create_basis_products(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_basis_products");
        create_basis_products(table_prefix);
    }

    void create_pi_basis_products(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_basis_products (id1 INTEGER, id2 INTEGER, prod TEXT, O SMALLINT, PRIMARY KEY(id1, id2));");
    }

    void drop_and_create_pi_basis_products(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_pi_basis_products");
        create_pi_basis_products(table_prefix);
    }

    void create_pi_basis_maps_S0() const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS S0_pi_basis_maps (id INTEGER PRIMARY KEY, to_C2 TEXT, to_Ceta TEXT, to_Cnu TEXT, to_Csigma TEXT, O_C2 SMALLINT, O_Ceta SMALLINT, O_Cnu SMALLINT, O_Csigma SMALLINT);");
    }

    void drop_and_create_pi_basis_maps_S0() const
    {
        drop_table("S0_pi_basis_maps");
        create_pi_basis_maps_S0();
    }

    void create_pi_basis_maps_Cof(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_basis_maps (id INTEGER PRIMARY KEY, to_S0 TEXT, O_S0 SMALLINT);");
    }

    void drop_and_create_pi_basis_maps_Cof(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_pi_basis_maps");
        create_pi_basis_maps_Cof(table_prefix);
    }

    void create_ss_products(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_products (id1 INTEGER, id2 INTEGER, prod TEXT, PRIMARY KEY(id1, id2));");
    }

    void drop_and_create_ss_products(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss_products");
        create_ss_products(table_prefix);
    }

    void create_ss_diff(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_diffs (src INTEGER PRIMARY KEY, r SMALLINT, tgt TEXT);");
    }

    void drop_and_create_ss_diff(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss_diffs");
        create_ss_diff(table_prefix);
    }

    void create_ss_nd(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_nd (src INTEGER PRIMARY KEY, r SMALLINT, tgt TEXT);");
    }

    void drop_and_create_ss_nd(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss_nd");
        create_ss_nd(table_prefix);
    }

    void create_ss_stable_levels(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_stable_levels (s SMALLINT, t SMALLINT, l SMALLINT, PRIMARY KEY(s, t));");
    }

    void drop_and_create_ss_stable_levels(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss_stable_levels");
        create_ss_stable_levels(table_prefix);
    }

    void load_basis_v2(const std::string& table_prefix, Mon1d& basis, AdamsDeg1d& deg_basis) const
    {
        Statement stmt(*this, "SELECT mon, s, t FROM " + table_prefix + "_basis ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            basis.push_back(myio::Deserialize<Mon>(stmt.column_str(0)));
            deg_basis.push_back(AdamsDeg(stmt.column_int(1), stmt.column_int(2)));
        }
        // myio::Logger::out() << "basis loaded from " << table_prefix << "_basis, size=" << basis.size() << '\n';
    }

    void load_basis_mod_v2(const std::string& table_prefix, MMod1d& basis, AdamsDeg1d& deg_basis) const
    {
        Statement stmt(*this, "SELECT mon, s, t FROM " + table_prefix + "_basis ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            basis.push_back(myio::Deserialize<MMod>(stmt.column_str(0)));
            deg_basis.push_back(AdamsDeg(stmt.column_int(1), stmt.column_int(2)));
        }
        // myio::Logger::out() << "basis loaded from " << table_prefix << "_basis, size=" << basis.size() << '\n';
    }

    void load_basis_ss_v2(const std::string& table_prefix, int2d& nodes_ss, AdamsDeg1d& deg_basis) const
    {
        Statement stmt(*this, "SELECT base, s, t FROM " + table_prefix + "_ss ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            nodes_ss.push_back(myio::Deserialize<int1d>(stmt.column_str(0)));
            deg_basis.push_back(AdamsDeg(stmt.column_int(1), stmt.column_int(2)));
        }
        // myio::Logger::out() << "basis loaded from " << table_prefix << "_ss, size=" << ss.size() << '\n';
    }
};

void serialize_ss(const Staircases& ss, int2d& seq_ss, AdamsDeg1d& seq_deg)
{
    auto degs = OrderDegsV2(ss);
    for (auto& d : degs) {
        for (auto& b : ss.at(d).basis) {
            seq_ss.push_back(b);
            seq_deg.push_back(d);
        }
    }
}

int main_plot(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Diagram diagram(diagram_name, DeduceFlag::no_op);
    diagram.DeduceTrivialDiffs();
    auto& rings = diagram.GetRings();
    auto& mods = diagram.GetModules();

    std::vector<MyDB> dbPlots;
    std::vector<std::string> names, paths;
    std::vector<int> isRing;
    GetAllDbNames(diagram_name, names, paths, isRing);

    dbPlots.reserve(paths.size());
    std::vector<std::map<AdamsDeg, int>> deg2ids;
    for (size_t i = 0; i < paths.size(); ++i) {
        std::string dbname = paths[i];
        dbname.insert(dbname.size() - 3, "_plot");
        dbPlots.emplace_back(dbname);
        MyDB db(paths[i]);
        deg2ids.push_back(db.load_basis_indices(fmt::format("{}_AdamsE2", names[i])));
    }

    for (auto& db : dbPlots)
        db.begin_transaction();

    /* ss_products */
    {
        int3d seq_ss_rings(rings.size()), seq_ss_mods(mods.size());
        AdamsDeg2d seq_deg_rings(rings.size()), seq_deg_mods(mods.size());
        for (size_t iRing = 0; iRing < rings.size(); ++iRing)
            serialize_ss(rings[iRing].nodes_ss.front(), seq_ss_rings[iRing], seq_deg_rings[iRing]);
        for (size_t iMod = 0; iMod < mods.size(); ++iMod)
            serialize_ss(mods[iMod].nodes_ss.front(), seq_ss_mods[iMod], seq_deg_mods[iMod]);

        int1d arr_factors = {1, 3, 7, 15, 23, 29, 33, 42};  //// TODO: change
        for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
            auto& ring = rings[iRing];
            auto& seq_ss_ring = seq_ss_rings[iRing];
            auto& seq_deg_ring = seq_deg_rings[iRing];
            dbPlots[iRing].drop_and_create_ss_products(fmt::format("{}_AdamsE2", ring.name));
            myio::Statement stmt(dbPlots[iRing], fmt::format("INSERT INTO {}_AdamsE2_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);", ring.name));
            for (int i : arr_factors) {
                for (size_t j = 0; j < seq_ss_ring.size(); ++j) {
                    const AdamsDeg deg_prod = seq_deg_ring[i] + seq_deg_ring[j];
                    if (deg_prod.t <= ring.t_max) {
                        Poly bi = Indices2Poly(seq_ss_ring[i], ring.basis.at(seq_deg_ring[i]));
                        Poly bj = Indices2Poly(seq_ss_ring[j], ring.basis.at(seq_deg_ring[j]));
                        Poly poly_prod = ring.gb.Reduce(bi * bj);
                        if (poly_prod) {
                            int1d prod = Poly2Indices(poly_prod, ring.basis.at(deg_prod));
                            prod = lina::GetInvImage(ring.nodes_ss.front().at(deg_prod).basis, prod);
                            for (size_t k = 0; k < prod.size(); ++k)
                                prod[k] += deg2ids[iRing][deg_prod];

                            stmt.bind_and_step((int)i, (int)j, myio::Serialize(prod));
                        }
                    }
                }
            }
        }

        for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
            auto& mod = mods[iMod];
            auto& seq_ss_mod = seq_ss_mods[iMod];
            auto& seq_deg_mod = seq_deg_mods[iMod];
            auto& ring = rings[mod.iRing];
            auto& seq_ss_ring = seq_ss_rings[mod.iRing];
            auto& seq_deg_ring = seq_deg_rings[mod.iRing];
            auto& nodes_ss = mod.nodes_ss;
            auto& basis = mod.basis;
            auto& gb = mod.gb;
            int t_max = mod.t_max;
            dbPlots[iMod + rings.size()].drop_and_create_ss_products(fmt::format("{}_AdamsE2", mod.name));
            myio::Statement stmt(dbPlots[iMod + rings.size()], fmt::format("INSERT INTO {}_AdamsE2_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);", mod.name));
            for (int i : arr_factors) {
                for (size_t j = 0; j < seq_ss_mod.size(); ++j) {
                    const AdamsDeg deg_prod = seq_deg_ring[i] + seq_deg_mod[j];
                    if (deg_prod.t <= t_max) {
                        Poly ai = Indices2Poly(seq_ss_ring[i], ring.basis.at(seq_deg_ring[i]));
                        if (ai.data.size() != 1)
                            throw MyException(0x8f0340d6U, "the factor is supposed to be the only one in its degree.");
                        Mod bj = Indices2Mod(seq_ss_mod[j], basis.at(seq_deg_mod[j]));
                        Mod x_prod = gb.Reduce(ai.GetLead() * bj);
                        if (x_prod) {
                            int1d prod = Mod2Indices(x_prod, basis.at(deg_prod));
                            prod = lina::GetInvImage(nodes_ss.front().at(deg_prod).basis, prod);
                            for (size_t l = 0; l < prod.size(); ++l)
                                prod[l] += deg2ids[iMod + rings.size()][deg_prod];

                            stmt.bind_and_step((int)i, (int)j, myio::Serialize(prod));
                        }
                    }
                }
            }
        }
    }

    /* ss_diffs */
    for (size_t iCw = 0; iCw < dbPlots.size(); ++iCw) {
        auto& nodes_ss = iCw < rings.size() ? rings[iCw].nodes_ss : mods[iCw - rings.size()].nodes_ss;
        auto& name = iCw < rings.size() ? rings[iCw].name : mods[iCw - rings.size()].name;

        dbPlots[iCw].drop_and_create_ss_diff(fmt::format("{}_AdamsE2", name));
        myio::Statement stmt(dbPlots[iCw], fmt::format("INSERT INTO {}_AdamsE2_ss_diffs (src, r, tgt) VALUES (?1, ?2, ?3);", name));
        for (auto& [deg, basis_ss_d] : nodes_ss.front()) {
            for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
                int src = deg2ids[iCw][deg] + (int)i;
                if (basis_ss_d.levels[i] > 9800 && basis_ss_d.diffs[i] != NULL_DIFF) {
                    int r = LEVEL_MAX - basis_ss_d.levels[i];
                    AdamsDeg deg_tgt = deg + AdamsDeg(r, r - 1);
                    int1d tgt = lina::GetInvImage(nodes_ss.front().at(deg_tgt).basis, basis_ss_d.diffs[i]);
                    for (size_t l = 0; l < tgt.size(); ++l)
                        tgt[l] += deg2ids[iCw][deg_tgt];

                    stmt.bind_and_step(src, r, myio::Serialize(tgt));
                }
            }
        }
    }

    /* ss_nd */
    for (size_t iCw = 0; iCw < dbPlots.size(); ++iCw) {
        auto& name = iCw < rings.size() ? rings[iCw].name : mods[iCw - rings.size()].name;
        auto& nodes_ss = iCw < rings.size() ? rings[iCw].nodes_ss : mods[iCw - rings.size()].nodes_ss;
        int t_max = iCw < rings.size() ? rings[iCw].t_max : mods[iCw - rings.size()].t_max;
        dbPlots[iCw].drop_and_create_ss_nd(fmt::format("{}_AdamsE2", name));
        myio::Statement stmt(dbPlots[iCw], fmt::format("INSERT INTO {}_AdamsE2_ss_nd (src, r, tgt) VALUES (?1, ?2, ?3);", name));

        auto& degs = iCw < rings.size() ? rings[iCw].degs_basis_order_by_stem : mods[iCw - rings.size()].degs_basis_order_by_stem;

        for (AdamsDeg deg : degs) {
            NullDiff1d nds;
            diagram.CacheNullDiffs(nodes_ss, t_max, deg, DeduceFlag::no_op, nds);
            for (auto& nd : nds) {
                auto& basis_ss_d = nodes_ss.front().at(deg);
                size_t index = 0;
                while (index < basis_ss_d.basis.size() && basis_ss_d.basis[index] != nd.x)
                    ++index;
                if (index == basis_ss_d.basis.size())
                    throw MyException(0xef63215fU, "nd could be wrong");
                int src = deg2ids[iCw][deg] + (int)index;
                if (basis_ss_d.levels[index] > 9800) {
                    int1d tgt;
                    for (int j = nd.first; j < nd.first + nd.count; ++j)
                        tgt.push_back(j);
                    int r = LEVEL_MAX - basis_ss_d.levels[index];
                    AdamsDeg deg_tgt = deg + AdamsDeg(r, r - 1);
                    for (size_t j = 0; j < tgt.size(); ++j)
                        tgt[j] += deg2ids[iCw][deg_tgt];

                    stmt.bind_and_step(src, r, myio::Serialize(tgt));
                }
            }
        }
    }

    /* ss_stable_levels */
    for (size_t iCw = 0; iCw < dbPlots.size(); ++iCw) {
        auto& name = iCw < rings.size() ? rings[iCw].name : mods[iCw - rings.size()].name;
        auto& nodes_ss = iCw < rings.size() ? rings[iCw].nodes_ss : mods[iCw - rings.size()].nodes_ss;
        dbPlots[iCw].drop_and_create_ss_stable_levels(fmt::format("{}_AdamsE2", name));
        myio::Statement stmt(dbPlots[iCw], fmt::format("INSERT INTO {}_AdamsE2_ss_stable_levels (s, t, l) VALUES (?1, ?2, ?3);", name));

        for (auto& [deg, basis_ss_d] : nodes_ss.front()) {
            int level = diagram.GetFirstFixedLevelForPlot(nodes_ss, deg);
            stmt.bind_and_step(deg.s, deg.t, level);
        }
    }

    for (auto& db : dbPlots)
        db.end_transaction();
    return 0;
}

void ToIndices(const algZ::Poly& x, const PiBasis& basis, const std::map<AdamsDeg, int>& pi_deg2id, int stem, int t_max, int1d& result, int& O)
{
    for (auto& m : x.data) {
        AdamsDeg d = AdamsDeg(m.fil(), stem + m.fil());
        if (m.IsUnKnown() || d.t > t_max) {
            O = m.fil();
            break;
        }
        else {
            const auto& basis_d = basis.at(d).nodes_pi_basis;
            const auto& p = std::lower_bound(basis_d.begin(), basis_d.end(), m);
            int index = (int)(p - basis_d.begin());
            result.push_back(pi_deg2id.at(d) + index);
        }
    }
}

void ToIndices(const algZ::Mod& x, const PiBasisMod& basis, const std::map<AdamsDeg, int>& pi_deg2id, int stem, int t_max, int1d& result, int& O)
{
    for (auto& m : x.data) {
        AdamsDeg d = AdamsDeg(m.fil(), stem + m.fil());
        if (m.IsUnKnown() || d.t > t_max) {
            O = m.fil();
            break;
        }
        else {
            const auto& basis_d = basis.at(d).nodes_pi_basis;
            const auto& p = std::lower_bound(basis_d.begin(), basis_d.end(), m);
            int index = (int)(p - basis_d.begin());
            result.push_back(pi_deg2id.at(d) + index);
        }
    }
}

int main_plot_htpy(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Diagram diagram(diagram_name, DeduceFlag::homotopy);
    /* pi_basis_products */
    int count_ss = 0, count_homotopy = 0;
    diagram.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_homotopy, 0);
    diagram.DeduceTrivialExtensions(0);

    auto& ssS0 = diagram.GetRings();
    auto& ssCofs = diagram.GetModules();
    // auto& all_basis_ss = diagram.GetAllBasisSs();

    // std::vector<MyDB> dbPlots;
    // dbPlots.reserve(dbnames.size()); /* To avoid reallocation is critical */
    // std::vector<std::string> tablesE2, complexNames;
    // std::vector<std::map<AdamsDeg, int>> pi_deg2ids(dbnames.size());
    // for (size_t i = 0; i < dbnames.size(); ++i) {
    //     std::string dbname = dbnames[i];
    //     dbname.insert(dbname.size() - 3, "_plot");
    //     dbPlots.emplace_back(dbname);
    //     tablesE2.push_back(GetE2TablePrefix(dbname));
    //     complexNames.push_back(GetComplexName(dbname));
    //     MyDB db(dbnames[i]);
    // }

    // for (auto& db : dbPlots)
    //     db.begin_transaction();

    //{
    //    algZ::Mon1d S0_pi_basis;
    //    AdamsDeg1d S0_pi_basis_deg;

    //    int pi_index = 0;
    //    AdamsDeg1d S0_degs = OrderDegsV2(diagram.GetRings().nodes_pi_basis.front());
    //    for (AdamsDeg d : S0_degs) {
    //        auto& basis_d = diagram.GetRings().nodes_pi_basis.front().at(d).nodes_pi_basis;
    //        pi_deg2ids[0][d] = pi_index;
    //        pi_index += (int)basis_d.size();
    //        for (auto& b : basis_d) {
    //            S0_pi_basis.push_back(b);
    //            S0_pi_basis_deg.push_back(d);
    //        }
    //    }

    //    algZ::MMod2d Cofs_pi_basis(diagram.GetModules().size());
    //    AdamsDeg2d Cofs_pi_basis_deg(diagram.GetModules().size());

    //    for (size_t iCof = 0; iCof < diagram.GetModules().size(); ++iCof) {
    //        auto& Cof = diagram.GetModules()[iCof];
    //        int pi_index = 0;
    //        AdamsDeg1d Cof_degs = OrderDegsV2(Cof.nodes_pi_basis.front());
    //        for (AdamsDeg d : Cof_degs) {
    //            auto& basis_d = Cof.nodes_pi_basis.front().at(d).nodes_pi_basis;
    //            pi_deg2ids[iCof + 1][d] = pi_index;
    //            pi_index += (int)basis_d.size();
    //            for (auto& b : basis_d) {
    //                Cofs_pi_basis[iCof].push_back(b);
    //                Cofs_pi_basis_deg[iCof].push_back(d);
    //            }
    //        }
    //    }

    //    const int1d arr_factors = {1, 3, 7, 15, 23, 29, 33, 39, 40, 42, 47};
    //    {
    //        auto& pi_gb = diagram.GetRings().pi_gb;
    //        int t_max = diagram.GetRings().t_max;
    //        dbPlots[0].drop_and_create_pi_basis_products(complexNames[0]);
    //        myio::Statement stmt(dbPlots[0], "INSERT INTO " + complexNames[0] + "_pi_basis_products (id1, id2, prod, O) VALUES (?1, ?2, ?3, ?4);");
    //        for (int i : arr_factors) {
    //            for (size_t j = 0; j < S0_pi_basis.size(); ++j) {
    //                const AdamsDeg deg_prod = S0_pi_basis_deg[i] + S0_pi_basis_deg[j];
    //                if (deg_prod.t <= t_max) {
    //                    algZ::Poly poly_prod = pi_gb.ReduceV2(S0_pi_basis[i] * S0_pi_basis[j]);
    //                    int1d prod;
    //                    int O = -1;
    //                    ToIndices(poly_prod, diagram.GetRings().nodes_pi_basis.front(), pi_deg2ids[0], deg_prod.stem(), t_max, prod, O);

    //                    if (O != -1 || !prod.empty()) {
    //                        stmt.bind_and_step((int)i, (int)j, myio::Serialize(prod), O);
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    dbPlots[0].drop_and_create_pi_basis_maps_S0();
    //    for (size_t iSS = 1; iSS < dbPlots.size(); ++iSS) {
    //        size_t iCof = iSS - 1;
    //        auto& pi_gb = ssCofs[iCof].pi_gb;
    //        int t_max = ssCofs[iCof].t_max;

    //        /* module structure */
    //        {
    //            dbPlots[iSS].drop_and_create_pi_basis_products(complexNames[iSS]);
    //            myio::Statement stmt(dbPlots[iSS], "INSERT INTO " + complexNames[iSS] + "_pi_basis_products (id1, id2, prod, O) VALUES (?1, ?2, ?3, ?4);");
    //            for (int i : arr_factors) {
    //                for (size_t j = 0; j < Cofs_pi_basis[iCof].size(); ++j) {
    //                    const AdamsDeg deg_prod = S0_pi_basis_deg[i] + Cofs_pi_basis_deg[iCof][j];
    //                    if (deg_prod.t <= t_max) {
    //                        algZ::Mod x_prod = pi_gb.ReduceV2(S0_pi_basis[i] * Cofs_pi_basis[iCof][j]);
    //                        int1d prod;
    //                        int O = -1;
    //                        ToIndices(x_prod, diagram.GetModules()[iCof].nodes_pi_basis.front(), pi_deg2ids[iSS], deg_prod.stem(), t_max, prod, O);

    //                        if (O != -1 || !prod.empty()) {
    //                            stmt.bind_and_step((int)i, (int)j, myio::Serialize(prod), O);
    //                        }
    //                    }
    //                }
    //            }
    //        }

    //        /* S0->Cof */
    //        {
    //            myio::Statement stmt(dbPlots[0], "INSERT INTO S0_pi_basis_maps (id, to_" + complexNames[iSS] + ", O_" + complexNames[iSS] + ") VALUES (?1, ?2, ?3) ON CONFLICT DO UPDATE SET to_" + complexNames[iSS] + "=excluded.to_"
    //                                                 + complexNames[iSS] + ", O_" + complexNames[iSS] + "=excluded.O_" + complexNames[iSS]);
    //            for (size_t i = 0; i < S0_pi_basis.size(); ++i) {
    //                const AdamsDeg deg = S0_pi_basis_deg[i];
    //                if (deg.t <= t_max) {
    //                    algZ::Mod x = pi_gb.ReduceV2(algZ::Mod(S0_pi_basis[i], 0, 0));
    //                    int1d image;
    //                    int O = -1;
    //                    ToIndices(x, diagram.GetModules()[iCof].nodes_pi_basis.front(), pi_deg2ids[iSS], deg.stem(), t_max, image, O);

    //                    if (O != -1 || !image.empty()) {
    //                        stmt.bind_and_step((int)i, myio::Serialize(image), O);
    //                    }
    //                }
    //            }
    //        }

    //        /* Cof->S0 */
    //        {
    //            dbPlots[iSS].drop_and_create_pi_basis_maps_Cof(complexNames[iSS]);
    //            myio::Statement stmt(dbPlots[iSS], "INSERT INTO " + complexNames[iSS] + "_pi_basis_maps (id, to_S0, O_S0) VALUES (?1, ?2, ?3)");
    //            for (size_t i = 0; i < Cofs_pi_basis[iCof].size(); ++i) {
    //                const AdamsDeg deg = Cofs_pi_basis_deg[iCof][i] - ssCofs[iCof].deg_qt;
    //                if (deg.t <= ssS0.t_max) {
    //                    algZ::Poly a = ssS0.pi_gb.ReduceV2(algZ::subs(Cofs_pi_basis[iCof][i], ssCofs[iCof].nodes_pi_qt.back(), ssCofs[iCof].pi_gb.v_degs()));
    //                    diagram.ExtendRelRing(deg.stem(), a);
    //                    int1d image;
    //                    int O = -1;
    //                    ToIndices(a, diagram.GetRings().nodes_pi_basis.front(), pi_deg2ids[0], deg.stem(), t_max, image, O);

    //                    if (O != -1 || !image.empty()) {
    //                        stmt.bind_and_step((int)i, myio::Serialize(image), O);
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}

    // for (auto& db : dbPlots)
    //     db.end_transaction();

    // diagram.save(dbnames, DeduceFlag::homotopy);
    return 0;
}
