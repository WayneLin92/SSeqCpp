#include "algebras/linalg.h"
#include "main.h"

using namespace alg2;

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
        // myio::Logger::out() << "basis loaded from " << table_prefix << "_ss, size=" << nodes_ss.size() << '\n';
    }
};

int main_basis_prod(int argc, char** argv, int index)
{
    std::string db_filename = "S0_AdamsSS_t259.db";
    std::string table_prefix = GetE2TablePrefix(db_filename);

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate the basis_prod table for S0\n";
        std::cout << "Usage:\n  ss basis_prod [db_filename] [table_prefix]\n\n";  // TODO: delete tables

        std::cout << "Default values:\n";
        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.1 (2022-7-16)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_prefix", table_prefix))
        return index;

    MyDB dbIn(db_filename);
    Poly1d polys = dbIn.load_gb(table_prefix, DEG_MAX);
    auto basis = dbIn.load_basis(table_prefix);
    auto ids = dbIn.load_basis_indices(table_prefix);
    Mon1d basis_v2;
    AdamsDeg1d deg_basis;
    dbIn.load_basis_v2(table_prefix, basis_v2, deg_basis);

    int t_max = basis.rbegin()->first.t;
    Groebner gb(t_max, {}, polys);

    std::string db_out = db_filename;
    auto p = db_out.insert(db_out.size() - 3, "_plot");
    MyDB dbOut(db_out);

    dbOut.begin_transaction();
    dbOut.drop_and_create_basis_products(table_prefix);
    myio::Statement stmt(dbOut, "INSERT INTO " + table_prefix + "_basis_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
    for (size_t i = 0; i < basis_v2.size(); ++i) {
        std::cout << "i=" << i << "          \r";
        for (size_t j = i; j < basis_v2.size(); ++j) {
            const AdamsDeg deg_prod = deg_basis[i] + deg_basis[j];
            if (deg_prod.t <= t_max) {
                Poly poly_prod = gb.Reduce(Poly(basis_v2[i]) * basis_v2[j]);
                int1d prod = Poly2Indices(poly_prod, basis[deg_prod]);
                for (size_t k = 0; k < prod.size(); ++k)
                    prod[k] += ids[deg_prod];

                stmt.bind_int(1, (int)i);
                stmt.bind_int(2, (int)j);
                stmt.bind_str(3, myio::Serialize(prod));
                stmt.step_and_reset();
            }
        }
    }
    dbOut.end_transaction();
    return 0;
}

int main_mod_basis_prod(int argc, char** argv, int index)
{
    std::string db_S0 = "S0_AdamsSS_t259.db";
    std::string table_S0 = GetE2TablePrefix(db_S0);
    std::string db_complex = "S0_AdamsSS_t259.db";
    std::string table_complex = GetE2TablePrefix(db_complex);

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate the basis_prod table for a complex\n";
        std::cout << "Usage:\n  ss mod basis_prod <db_complex> <table_complex> [db_S0] [table_S0]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  db_S0 = " << db_S0 << "\n";
        std::cout << "  table_S0 = " << table_S0 << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "db_complex", db_complex))
        return index;
    if (myio::load_arg(argc, argv, ++index, "table_complex", table_complex))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_S0", db_S0))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_S0", table_S0))
        return index;

    MyDB dbMod(db_complex);
    Mod1d xs_mod = dbMod.load_gb_mod(table_complex, DEG_MAX);
    auto basis_mod = dbMod.load_basis_mod(table_complex);
    auto ids_mod = dbMod.load_basis_indices(table_complex);
    MMod1d basis_mod_v2;
    AdamsDeg1d deg_basis_mod;
    dbMod.load_basis_mod_v2(table_complex, basis_mod_v2, deg_basis_mod);
    int t_max = basis_mod.rbegin()->first.t;

    MyDB dbS0(db_S0);
    Poly1d polys = dbS0.load_gb(table_S0, t_max);
    auto basis = dbS0.load_basis(table_S0);
    Mon1d basis_v2;
    AdamsDeg1d deg_basis;
    dbS0.load_basis_v2(table_S0, basis_v2, deg_basis);
    Groebner gb(t_max, {}, polys);

    GroebnerMod gbm(&gb, t_max, {}, xs_mod);

    for (size_t i = 0; i < 10 && i < gbm.data().size(); ++i)
        std::cout << gbm.data()[i].Str() << '\n';

    std::string db_out = db_complex;
    auto p = db_out.insert(db_out.size() - 3, "_plot");
    MyDB dbOut(db_out);

    dbOut.begin_transaction();
    dbOut.drop_and_create_basis_products(table_complex);
    myio::Statement stmt(dbOut, "INSERT INTO " + table_complex + "_basis_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
    for (size_t i = 0; i < basis_v2.size(); ++i) {
        std::cout << "i=" << i << "          \r";
        for (size_t j = 0; j < basis_mod_v2.size(); ++j) {
            const AdamsDeg deg_prod = deg_basis[i] + deg_basis_mod[j];
            if (deg_prod.t <= t_max) {
                Mod x_prod = gbm.Reduce(basis_v2[i] * basis_mod_v2[j]);
                int1d prod = Mod2Indices(x_prod, basis_mod[deg_prod]);
                for (size_t k = 0; k < prod.size(); ++k)
                    prod[k] += ids_mod[deg_prod];

                stmt.bind_int(1, (int)i);
                stmt.bind_int(2, (int)j);
                stmt.bind_str(3, myio::Serialize(prod));
                stmt.step_and_reset();
            }
        }
    }
    dbOut.end_transaction();
    return 0;
}

int main_plot(int argc, char** argv, int index)
{
    std::string selector = "default";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate tables: ss_prod, ss_diff, ss_nd, ss_stable_levels for plotting\n";
        std::cout << "Usage:\n  ss plot [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    Diagram diagram(dbnames, DeduceFlag::no_op);
    diagram.DeduceTrivialDiffs();
    auto& ssS0 = diagram.GetS0();
    auto& ssCof = diagram.GetCofs();
    auto& all_basis_ss = diagram.GetAllBasisSs();

    std::vector<MyDB> dbPlots;
    dbPlots.reserve(dbnames.size());
    std::vector<std::string> tablesE2, complexNames;
    std::vector<std::map<AdamsDeg, int>> deg2ids;
    for (size_t i = 0; i < dbnames.size(); ++i) {
        std::string dbname = dbnames[i];
        dbname.insert(dbname.size() - 3, "_plot");
        dbPlots.emplace_back(dbname);
        tablesE2.push_back(GetE2TablePrefix(dbname));
        complexNames.push_back(GetComplexName(dbname));
        MyDB db(dbnames[i]);
        deg2ids.push_back(db.load_basis_indices(tablesE2.back()));
    }

    for (auto& db : dbPlots)
        db.begin_transaction();

    /* ss_products */
    {
        int3d basis_sss(all_basis_ss.size());
        AdamsDeg2d deg_bases(all_basis_ss.size());
        for (size_t i = 0; i < all_basis_ss.size(); ++i) {
            auto& nodes_ss = all_basis_ss[i]->front();
            auto degs = OrderDegsV2(nodes_ss);
            for (auto& d : degs) {
                for (auto& b : nodes_ss.at(d).basis) {
                    basis_sss[i].push_back(b);
                    deg_bases[i].push_back(d);
                }
            }
        }

        int1d arr_factors = {1, 3, 7, 15, 23, 29, 33, 42};
        {
            dbPlots[0].drop_and_create_ss_products(tablesE2[0]);
            myio::Statement stmt(dbPlots[0], "INSERT INTO " + tablesE2[0] + "_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
            for (int i : arr_factors) {
                for (size_t j = 0; j < basis_sss[0].size(); ++j) {
                    const AdamsDeg deg_prod = deg_bases[0][i] + deg_bases[0][j];
                    if (deg_prod.t <= ssS0.t_max) {
                        Poly bi = Indices2Poly(basis_sss[0][i], ssS0.basis.at(deg_bases[0][i]));
                        Poly bj = Indices2Poly(basis_sss[0][j], ssS0.basis.at(deg_bases[0][j]));
                        Poly poly_prod = ssS0.gb.Reduce(bi * bj);
                        if (poly_prod) {
                            int1d prod = Poly2Indices(poly_prod, ssS0.basis.at(deg_prod));
                            prod = lina::GetInvImage(ssS0.nodes_ss.front().at(deg_prod).basis, prod);
                            for (size_t k = 0; k < prod.size(); ++k)
                                prod[k] += deg2ids[0][deg_prod];

                            stmt.bind_int(1, (int)i);
                            stmt.bind_int(2, (int)j);
                            stmt.bind_str(3, myio::Serialize(prod));
                            stmt.step_and_reset();
                        }
                    }
                }
            }
        }

        for (size_t k = 1; k < dbPlots.size(); ++k) {
            auto& nodes_ss = *all_basis_ss[k];
            auto& basis = ssCof[k - 1].basis;
            auto& gb = ssCof[k - 1].gb;
            int t_max = ssCof[k - 1].t_max;
            dbPlots[k].drop_and_create_ss_products(tablesE2[k]);
            myio::Statement stmt(dbPlots[k], "INSERT INTO " + tablesE2[k] + "_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
            for (int i : arr_factors) {
                for (size_t j = 0; j < basis_sss[k].size(); ++j) {
                    const AdamsDeg deg_prod = deg_bases[0][i] + deg_bases[k][j];
                    if (deg_prod.t <= t_max) {
                        Poly ai = Indices2Poly(basis_sss[0][i], ssS0.basis.at(deg_bases[0][i]));
                        if (ai.data.size() != 1)
                            throw MyException(0x8f0340d6U, "the factor is supposed to be the only one in its degree.");
                        Mod bj = Indices2Mod(basis_sss[k][j], basis.at(deg_bases[k][j]));
                        Mod x_prod = gb.Reduce(ai.GetLead() * bj);
                        if (x_prod) {
                            int1d prod = Mod2Indices(x_prod, basis.at(deg_prod));
                            prod = lina::GetInvImage(nodes_ss.front().at(deg_prod).basis, prod);
                            for (size_t l = 0; l < prod.size(); ++l)
                                prod[l] += deg2ids[k][deg_prod];

                            stmt.bind_int(1, (int)i);
                            stmt.bind_int(2, (int)j);
                            stmt.bind_str(3, myio::Serialize(prod));
                            stmt.step_and_reset();
                        }
                    }
                }
            }
        }
    }

    /* ss_diffs */
    for (size_t k = 0; k < dbPlots.size(); ++k) {
        dbPlots[k].drop_and_create_ss_diff(tablesE2[k]);
        auto& nodes_ss = *all_basis_ss[k];
        myio::Statement stmt(dbPlots[k], "INSERT INTO " + tablesE2[k] + "_ss_diffs (src, r, tgt) VALUES (?1, ?2, ?3);");
        for (auto& [deg, basis_ss_d] : nodes_ss.front()) {
            for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
                int src = deg2ids[k][deg] + (int)i;
                if (basis_ss_d.levels[i] > 9800 && basis_ss_d.diffs[i] != NULL_DIFF) {
                    int r = LEVEL_MAX - basis_ss_d.levels[i];
                    AdamsDeg deg_tgt = deg + AdamsDeg(r, r - 1);
                    int1d tgt = lina::GetInvImage(nodes_ss.front().at(deg_tgt).basis, basis_ss_d.diffs[i]);
                    for (size_t l = 0; l < tgt.size(); ++l)
                        tgt[l] += deg2ids[k][deg_tgt];

                    stmt.bind_int(1, src);
                    stmt.bind_int(2, r);
                    stmt.bind_str(3, myio::Serialize(tgt));
                    stmt.step_and_reset();
                }
            }
        }
    }

    /* ss_nd */
    for (size_t iSS = 0; iSS < dbPlots.size(); ++iSS) {
        dbPlots[iSS].drop_and_create_ss_nd(tablesE2[iSS]);
        myio::Statement stmt(dbPlots[iSS], "INSERT INTO " + tablesE2[iSS] + "_ss_nd (src, r, tgt) VALUES (?1, ?2, ?3);");

        auto& degs = [&]() {
            if (iSS == 0)
                return diagram.GetS0().degs_basis_order_by_stem;
            else
                return diagram.GetCofs()[iSS - 1].degs_basis_order_by_stem;
        }();

        for (AdamsDeg deg : degs) {
            NullDiff1d nds;
            diagram.CacheNullDiffs(iSS, deg, DeduceFlag::no_op, nds);
            for (auto& nd : nds) {
                auto& basis_ss_d = all_basis_ss[iSS]->front().at(deg);
                size_t index = 0;
                while (index < basis_ss_d.basis.size() && basis_ss_d.basis[index] != nd.x)
                    ++index;
                if (index == basis_ss_d.basis.size())
                    throw MyException(0xef63215fU, "nd could be wrong");
                int src = deg2ids[iSS][deg] + (int)index;
                if (basis_ss_d.levels[index] > 9800) {
                    int1d tgt;
                    for (int j = nd.first; j < nd.first + nd.count; ++j)
                        tgt.push_back(j);
                    int r = LEVEL_MAX - basis_ss_d.levels[index];
                    AdamsDeg deg_tgt = deg + AdamsDeg(r, r - 1);
                    for (size_t j = 0; j < tgt.size(); ++j)
                        tgt[j] += deg2ids[iSS][deg_tgt];

                    stmt.bind_int(1, src);
                    stmt.bind_int(2, r);
                    stmt.bind_str(3, myio::Serialize(tgt));
                    stmt.step_and_reset();
                }
            }
        }
    }

    /* ss_stable_levels */
    for (size_t k = 0; k < dbPlots.size(); ++k) {
        dbPlots[k].drop_and_create_ss_stable_levels(tablesE2[k]);
        myio::Statement stmt(dbPlots[k], "INSERT INTO " + tablesE2[k] + "_ss_stable_levels (s, t, l) VALUES (?1, ?2, ?3);");
        auto& nodes_ss = *all_basis_ss[k];

        for (auto& [deg, basis_ss_d] : nodes_ss.front()) {
            stmt.bind_int(1, deg.s);
            stmt.bind_int(2, deg.t);
            int level = diagram.GetFirstFixedLevelForPlot(nodes_ss, deg);
            stmt.bind_int(3, level);
            stmt.step_and_reset();
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
            auto& basis_d = basis.at(d).pi_basis;
            auto& p = std::lower_bound(basis_d.begin(), basis_d.end(), m);
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
            auto& basis_d = basis.at(d).pi_basis;
            auto& p = std::lower_bound(basis_d.begin(), basis_d.end(), m);
            int index = (int)(p - basis_d.begin());
            result.push_back(pi_deg2id.at(d) + index);
        }
    }
}

int main_plotpi(int argc, char** argv, int index)
{
    std::string selector = "default";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate tables: pi_basis_products pi_basis_maps for plotting\n";
        std::cout << "Usage:\n  ss plotpi [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    Diagram diagram(dbnames, DeduceFlag::homotopy);
    /* pi_basis_products */
    int count_ss = 0, count_homotopy = 0;
    diagram.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_homotopy, 0);
    diagram.DeduceTrivialExtensions(0);

    auto& ssS0 = diagram.GetS0();
    auto& ssCofs = diagram.GetCofs();
    auto& all_basis_ss = diagram.GetAllBasisSs();

    std::vector<MyDB> dbPlots;
    dbPlots.reserve(dbnames.size()); /* To avoid reallocation is critical */
    std::vector<std::string> tablesE2, complexNames;
    std::vector<std::map<AdamsDeg, int>> pi_deg2ids(dbnames.size());
    for (size_t i = 0; i < dbnames.size(); ++i) {
        std::string dbname = dbnames[i];
        dbname.insert(dbname.size() - 3, "_plot");
        dbPlots.emplace_back(dbname);
        tablesE2.push_back(GetE2TablePrefix(dbname));
        complexNames.push_back(GetComplexName(dbname));
        MyDB db(dbnames[i]);
    }

    for (auto& db : dbPlots)
        db.begin_transaction();

    {
        algZ::Mon1d S0_pi_basis;
        AdamsDeg1d S0_pi_basis_deg;

        int pi_index = 0;
        AdamsDeg1d S0_degs = OrderDegsV2(diagram.GetS0().pi_basis.front());
        for (AdamsDeg d : S0_degs) {
            auto& basis_d = diagram.GetS0().pi_basis.front().at(d).pi_basis;
            pi_deg2ids[0][d] = pi_index;
            pi_index += (int)basis_d.size();
            for (auto& b : basis_d) {
                S0_pi_basis.push_back(b);
                S0_pi_basis_deg.push_back(d);
            }
        }

        algZ::MMod2d Cofs_pi_basis(diagram.GetCofs().size());
        AdamsDeg2d Cofs_pi_basis_deg(diagram.GetCofs().size());

        for (size_t iCof = 0; iCof < diagram.GetCofs().size(); ++iCof) {
            auto& Cof = diagram.GetCofs()[iCof];
            int pi_index = 0;
            AdamsDeg1d Cof_degs = OrderDegsV2(Cof.pi_basis.front());
            for (AdamsDeg d : Cof_degs) {
                auto& basis_d = Cof.pi_basis.front().at(d).pi_basis;
                pi_deg2ids[iCof + 1][d] = pi_index;
                pi_index += (int)basis_d.size();
                for (auto& b : basis_d) {
                    Cofs_pi_basis[iCof].push_back(b);
                    Cofs_pi_basis_deg[iCof].push_back(d);
                }
            }
        }

        const int1d arr_factors = {1, 3, 7, 15, 23, 29, 33, 39, 40, 42, 47};
        {
            auto& pi_gb = diagram.GetS0().pi_gb;
            int t_max = diagram.GetS0().t_max;
            dbPlots[0].drop_and_create_pi_basis_products(complexNames[0]);
            myio::Statement stmt(dbPlots[0], "INSERT INTO " + complexNames[0] + "_pi_basis_products (id1, id2, prod, O) VALUES (?1, ?2, ?3, ?4);");
            for (int i : arr_factors) {
                for (size_t j = 0; j < S0_pi_basis.size(); ++j) {
                    const AdamsDeg deg_prod = S0_pi_basis_deg[i] + S0_pi_basis_deg[j];
                    if (deg_prod.t <= t_max) {
                        algZ::Poly poly_prod = pi_gb.ReduceV2(S0_pi_basis[i] * S0_pi_basis[j]);
                        int1d prod;
                        int O = -1;
                        ToIndices(poly_prod, diagram.GetS0().pi_basis.front(), pi_deg2ids[0], deg_prod.stem(), t_max, prod, O);

                        if (O != -1 || !prod.empty()) {
                            stmt.bind_int(1, (int)i);
                            stmt.bind_int(2, (int)j);
                            stmt.bind_str(3, myio::Serialize(prod));
                            stmt.bind_int(4, O);
                            stmt.step_and_reset();
                        }
                    }
                }
            }
        }

        dbPlots[0].drop_and_create_pi_basis_maps_S0();
        for (size_t iSS = 1; iSS < dbPlots.size(); ++iSS) {
            size_t iCof = iSS - 1;
            auto& pi_gb = ssCofs[iCof].pi_gb;
            int t_max = ssCofs[iCof].t_max;

            /* module structure */
            {
                dbPlots[iSS].drop_and_create_pi_basis_products(complexNames[iSS]);
                myio::Statement stmt(dbPlots[iSS], "INSERT INTO " + complexNames[iSS] + "_pi_basis_products (id1, id2, prod, O) VALUES (?1, ?2, ?3, ?4);");
                for (int i : arr_factors) {
                    for (size_t j = 0; j < Cofs_pi_basis[iCof].size(); ++j) {
                        const AdamsDeg deg_prod = S0_pi_basis_deg[i] + Cofs_pi_basis_deg[iCof][j];
                        if (deg_prod.t <= t_max) {
                            algZ::Mod x_prod = pi_gb.ReduceV2(S0_pi_basis[i] * Cofs_pi_basis[iCof][j]);
                            int1d prod;
                            int O = -1;
                            ToIndices(x_prod, diagram.GetCofs()[iCof].pi_basis.front(), pi_deg2ids[iSS], deg_prod.stem(), t_max, prod, O);

                            if (O != -1 || !prod.empty()) {
                                stmt.bind_int(1, (int)i);
                                stmt.bind_int(2, (int)j);
                                stmt.bind_str(3, myio::Serialize(prod));
                                stmt.bind_int(4, O);
                                stmt.step_and_reset();
                            }
                        }
                    }
                }
            }

            /* S0->Cof */
            {
                myio::Statement stmt(dbPlots[0], "INSERT INTO S0_pi_basis_maps (id, to_" + complexNames[iSS] + ", O_" + complexNames[iSS] + ") VALUES (?1, ?2, ?3) ON CONFLICT DO UPDATE SET to_" + complexNames[iSS] + "=excluded.to_"
                                                     + complexNames[iSS] + ", O_" + complexNames[iSS] + "=excluded.O_" + complexNames[iSS]);
                for (size_t i = 0; i < S0_pi_basis.size(); ++i) {
                    const AdamsDeg deg = S0_pi_basis_deg[i];
                    if (deg.t <= t_max) {
                        algZ::Mod x = pi_gb.ReduceV2(algZ::Mod(S0_pi_basis[i], 0, 0));
                        int1d image;
                        int O = -1;
                        ToIndices(x, diagram.GetCofs()[iCof].pi_basis.front(), pi_deg2ids[iSS], deg.stem(), t_max, image, O);

                        if (O != -1 || !image.empty()) {
                            stmt.bind_int(1, (int)i);
                            stmt.bind_str(2, myio::Serialize(image));
                            stmt.bind_int(3, O);
                            stmt.step_and_reset();
                        }
                    }
                }
            }

            /* Cof->S0 */
            {
                dbPlots[iSS].drop_and_create_pi_basis_maps_Cof(complexNames[iSS]);
                myio::Statement stmt(dbPlots[iSS], "INSERT INTO " + complexNames[iSS] + "_pi_basis_maps (id, to_S0, O_S0) VALUES (?1, ?2, ?3)");
                for (size_t i = 0; i < Cofs_pi_basis[iCof].size(); ++i) {
                    const AdamsDeg deg = Cofs_pi_basis_deg[iCof][i] - ssCofs[iCof].deg_qt;
                    if (deg.t <= ssS0.t_max) {
                        algZ::Poly a = ssS0.pi_gb.ReduceV2(algZ::subsMod(Cofs_pi_basis[iCof][i], ssCofs[iCof].pi_qt.back(), ssCofs[iCof].pi_gb.v_degs()));
                        diagram.ExtendRelS0(deg.stem(), a);
                        int1d image;
                        int O = -1;
                        ToIndices(a, diagram.GetS0().pi_basis.front(), pi_deg2ids[0], deg.stem(), t_max, image, O);

                        if (O != -1 || !image.empty()) {
                            stmt.bind_int(1, (int)i);
                            stmt.bind_str(2, myio::Serialize(image));
                            stmt.bind_int(3, O);
                            stmt.step_and_reset();
                        }
                    }
                }
            }
        }
    }

    for (auto& db : dbPlots)
        db.end_transaction();

    diagram.save(dbnames, DeduceFlag::homotopy);
    return 0;
}

int main_mod(int argc, char** argv, int index)
{
    std::string cmd;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  ss mod <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  basis_prod: Generate the basis_prod table for a complex\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "basis_prod")
        return main_mod_basis_prod(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << cmd << '\n';
    return 0;
}
