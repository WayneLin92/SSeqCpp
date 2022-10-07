#include "algebras/linalg.h"
#include "main.h"

using namespace alg;

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
        std::clog << "basis loaded from " << table_prefix << "_basis, size=" << basis.size() << '\n';
    }

    void load_basis_mod_v2(const std::string& table_prefix, MMod1d& basis, AdamsDeg1d& deg_basis) const
    {
        Statement stmt(*this, "SELECT mon, s, t FROM " + table_prefix + "_basis ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            basis.push_back(myio::Deserialize<MMod>(stmt.column_str(0)));
            deg_basis.push_back(AdamsDeg(stmt.column_int(1), stmt.column_int(2)));
        }
        std::clog << "basis loaded from " << table_prefix << "_basis, size=" << basis.size() << '\n';
    }

    void load_basis_ss_v2(const std::string& table_prefix, int2d& basis_ss, AdamsDeg1d& deg_basis) const
    {
        Statement stmt(*this, "SELECT base, s, t FROM " + table_prefix + "_ss ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            basis_ss.push_back(myio::Deserialize<int1d>(stmt.column_str(0)));
            deg_basis.push_back(AdamsDeg(stmt.column_int(1), stmt.column_int(2)));
        }
        std::clog << "basis loaded from " << table_prefix << "_ss, size=" << basis_ss.size() << '\n';
    }
};

int main_basis_prod(int argc, char** argv, int index)
{
    std::string db_filename = DB_DEFAULT;
    std::string table_prefix = GetTablePrefix(db_filename);

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate the basis_prod table for S0\n";
        std::cout << "Usage:\n  ss basis_prod [db_filename] [table_prefix]\n\n"; //TODO: delete tables

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
    auto ids = dbIn.load_indices(table_prefix);
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
    std::string db_S0 = DB_DEFAULT;
    std::string table_S0 = GetTablePrefix(db_S0);
    std::string db_complex = DB_DEFAULT;
    std::string table_complex = GetTablePrefix(db_complex);

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate the basis_prod table for a complex\n";
        std::cout << "Usage:\n  ss mod basis_prod <db_complex> <table_complex> [db_S0] [table_S0]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  db_S0 = " << db_S0 << "\n";
        std::cout << "  table_S0 = " << table_S0 << "\n\n";

        std::cout << "Version:\n  2.1 (2022-9-15)" << std::endl;
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
    auto ids_mod = dbMod.load_indices(table_complex);
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
    std::string db_S0 = DB_DEFAULT;
    std::vector<std::string> dbnames = {
        "C2_AdamsSS_t221.db",
        "Ceta_AdamsSS_t200.db",
        "Cnu_AdamsSS_t200.db",
        "Csigma_AdamsSS_t200.db",
    };

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate tables: ss_prod, ss_diff, ss_nd, ss_stable_levels for plotting\n";
        std::cout << "Usage:\n  ss plot [db_S0] [db_Cofib ...]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  db_S0 = " << db_S0 << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_S0", db_S0))
        return index;
    if (myio::load_args(argc, argv, ++index, "db_Cofib", dbnames))
        return index;
    dbnames.insert(dbnames.begin(), db_S0);


    Diagram diagram(dbnames);
    auto& ssS0 = diagram.GetS0();
    auto& ssCof = diagram.GetCofs();
    auto& all_basis_ss = diagram.GetAllBasisSs();
    auto& all_nd = diagram.GetAllNd();

    std::vector<MyDB> dbPlots;
    dbPlots.reserve(dbnames.size());
    std::vector<std::string> tables;
    std::vector<std::map<AdamsDeg, int>> deg2ids;
    for (size_t i = 0; i < dbnames.size(); ++i) {
        std::string dbname = dbnames[i];
        dbname.insert(dbname.size() - 3, "_plot");
        dbPlots.emplace_back(dbname);
        tables.push_back(GetTablePrefix(dbname));
        deg2ids.push_back(MyDB(dbnames[i]).load_indices(tables.back()));
    }

    for (auto& db : dbPlots)
        db.begin_transaction();

    /* ss_products */
    {
        int3d basis_sss(all_basis_ss.size());
        AdamsDeg2d deg_bases(all_basis_ss.size());
        for (size_t i = 0; i < all_basis_ss.size(); ++i) {
            auto& basis_ss = all_basis_ss[i]->front();
            auto degs = OrderDegsV2(basis_ss);
            for (auto& d : degs) {
                for (auto& b : basis_ss.at(d).basis_ind) {
                    basis_sss[i].push_back(b);
                    deg_bases[i].push_back(d);
                }
            }
        }

        int1d arr_factors = {1, 3, 7, 15, 23, 29, 33, 42};
        {
            dbPlots[0].drop_and_create_ss_products(tables[0]);
            myio::Statement stmt(dbPlots[0], "INSERT INTO " + tables[0] + "_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
            for (int i : arr_factors) {
                for (size_t j = 0; j < basis_sss[0].size(); ++j) {
                    const AdamsDeg deg_prod = deg_bases[0][i] + deg_bases[0][j];
                    if (deg_prod.t <= ssS0.t_max) {
                        Poly bi = Indices2Poly(basis_sss[0][i], ssS0.basis.at(deg_bases[0][i]));
                        Poly bj = Indices2Poly(basis_sss[0][j], ssS0.basis.at(deg_bases[0][j]));
                        Poly poly_prod = ssS0.gb.Reduce(bi * bj);
                        if (poly_prod) {
                            int1d prod = Poly2Indices(poly_prod, ssS0.basis.at(deg_prod));
                            prod = lina::GetInvImage(ssS0.basis_ss.front().at(deg_prod).basis_ind, prod);
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
            auto& basis_ss = *all_basis_ss[k];
            auto& basis = ssCof[k - 1].basis;
            auto& gb = ssCof[k - 1].gb;
            int t_max = ssCof[k - 1].t_max;
            dbPlots[k].drop_and_create_ss_products(tables[k]);
            myio::Statement stmt(dbPlots[k], "INSERT INTO " + tables[k] + "_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
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
                            prod = lina::GetInvImage(basis_ss.front().at(deg_prod).basis_ind, prod);
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
        dbPlots[k].drop_and_create_ss_diff(tables[k]);
        auto& basis_ss = *all_basis_ss[k];
        myio::Statement stmt(dbPlots[k], "INSERT INTO " + tables[k] + "_ss_diffs (src, r, tgt) VALUES (?1, ?2, ?3);");
        for (auto& [deg, basis_ss_d] : basis_ss.front()) {
            for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
                int src = deg2ids[k][deg] + (int)i;
                if (basis_ss_d.levels[i] > 9800 && basis_ss_d.diffs_ind[i] != int1d{-1}) {
                    int r = kLevelMax - basis_ss_d.levels[i];
                    AdamsDeg deg_tgt = deg + AdamsDeg(r, r - 1);
                    int1d tgt = lina::GetInvImage(basis_ss.front().at(deg_tgt).basis_ind, basis_ss_d.diffs_ind[i]);
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
    for (size_t k = 0; k < dbPlots.size(); ++k) {
        dbPlots[k].drop_and_create_ss_nd(tables[k]);
        myio::Statement stmt(dbPlots[k], "INSERT INTO " + tables[k] + "_ss_nd (src, r, tgt) VALUES (?1, ?2, ?3);");
        diagram.CacheNullDiffs(4, 127, 0);

        for (size_t i = 0; i < all_nd[k]->back().size(); ++i) {
            auto& nd = all_nd[k]->back()[i];
            auto& basis_ss_d = all_basis_ss[k]->front().at(nd.deg);
            size_t index = 0;
            while (index < basis_ss_d.basis_ind.size() && basis_ss_d.basis_ind[index] != nd.x)
                ++index;
            if (index == basis_ss_d.basis_ind.size())
                throw MyException(0xef63215fU, "nd could be wrong");
            int src = deg2ids[k][nd.deg] + (int)index;
            if (basis_ss_d.levels[index] > 9800) {
                int r = kLevelMax - basis_ss_d.levels[index];
                AdamsDeg deg_tgt = nd.deg + AdamsDeg(r, r - 1);
                int1d tgt;
                for (int j = nd.first; j < nd.first + nd.count; ++j)
                    tgt.push_back(j);
                for (size_t j = 0; j < tgt.size(); ++j)
                    tgt[j] += deg2ids[k][deg_tgt];

                stmt.bind_int(1, src);
                stmt.bind_int(2, r);
                stmt.bind_str(3, myio::Serialize(tgt));
                stmt.step_and_reset();
            }
        }
    }

    /* ss_stable_levels */
    for (size_t k = 0; k < dbPlots.size(); ++k) {
        dbPlots[k].drop_and_create_ss_stable_levels(tables[k]);
        myio::Statement stmt(dbPlots[k], "INSERT INTO " + tables[k] + "_ss_stable_levels (s, t, l) VALUES (?1, ?2, ?3);");
        auto& basis_ss = *all_basis_ss[k];

        for (auto& [deg, basis_ss_d] : basis_ss.front()) {
            stmt.bind_int(1, deg.s);
            stmt.bind_int(2, deg.t);
            int level = diagram.GetFirstFixedLevelForPlot(basis_ss, deg);
            stmt.bind_int(3, level);
            stmt.step_and_reset();
        }
    }

    for (auto& db : dbPlots)
        db.end_transaction();
    std::cout << "Done" << std::endl;
    return 0;
}

int main_mod(int argc, char** argv, int index)
{
    std::string cmd;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  ss mod <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";

        std::cout << "  basis_prod: Generate the basis_prod table for a complex\n\n";

        std::cout << "Version:\n  2.1 (2022-9-15)" << std::endl;
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