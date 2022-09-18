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
    std::string db_filename = db_ss_default;
    std::string table_prefix = table_ss_default;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate the basis_prod table for S0\n";
        std::cout << "Usage:\n  ss basis_prod [db_filename] [table_prefix]\n\n";

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
    std::string db_S0 = db_ss_default;
    std::string table_S0 = table_ss_default;
    std::string db_complex = db_ss_default;
    std::string table_complex = table_ss_default;

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
    std::string db_filename = db_ss_default;
    std::string table_prefix = table_ss_default;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate tables: ss_prod,ss_diff,ss_nd,ss_stable_levels for plotting\n";
        std::cout << "Usage:\n  ss plot [db_filename] [table_prefix]\n\n";

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
    const Staircases basis_ss = dbIn.load_basis_ss(table_prefix);
    auto ids = dbIn.load_indices(table_prefix);
    int2d basis_ss_v2;
    AdamsDeg1d deg_basis;
    dbIn.load_basis_ss_v2(table_prefix, basis_ss_v2, deg_basis);

    int t_max = deg_basis.back().t;
    Groebner gb(t_max, {}, polys);

    std::string db_out = db_filename;
    auto p = db_out.insert(db_out.size() - 3, "_plot");
    MyDB dbOut(db_out);

    dbOut.begin_transaction();

    /* ss_products */
    {
        dbOut.drop_and_create_ss_products(table_prefix);
        myio::Statement stmt(dbOut, "INSERT INTO " + table_prefix + "_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
        for (size_t i = 0; i < basis_ss_v2.size(); ++i) {
            std::cout << "i=" << i << "          \r";
            for (size_t j = i; j < basis_ss_v2.size(); ++j) {
                const AdamsDeg deg_prod = deg_basis[i] + deg_basis[j];
                if (deg_prod.t <= t_max) {
                    Poly bi = Indices2Poly(basis_ss_v2[i], basis[deg_basis[i]]);
                    Poly bj = Indices2Poly(basis_ss_v2[j], basis[deg_basis[j]]);
                    Poly poly_prod = gb.Reduce(bi * bj);
                    int1d prod = Poly2Indices(poly_prod, basis[deg_prod]);
                    if (!prod.empty())
                        prod = lina::GetInvImage(basis_ss.at(deg_prod).basis_ind, prod);
                    for (size_t k = 0; k < prod.size(); ++k)
                        prod[k] += ids[deg_prod];

                    stmt.bind_int(1, (int)i);
                    stmt.bind_int(2, (int)j);
                    stmt.bind_str(3, myio::Serialize(prod));
                    stmt.step_and_reset();
                }
            }
        }
    }

    /* ss_diffs */
    {
        dbOut.drop_and_create_ss_diff(table_prefix);
        myio::Statement stmt(dbOut, "INSERT INTO " + table_prefix + "_ss_diffs (src, r, tgt) VALUES (?1, ?2, ?3);");
        for (auto& [deg, basis_ss_d] : basis_ss) {
            for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
                int src = ids[deg] + (int)i;
                if (basis_ss_d.levels[i] > 9800 && basis_ss_d.diffs_ind[i] != int1d{-1}) {
                    int r = kLevelMax - basis_ss_d.levels[i];
                    AdamsDeg deg_tgt = deg + AdamsDeg(r, r - 1);
                    int1d tgt = lina::GetInvImage(basis_ss.at(deg_tgt).basis_ind, basis_ss_d.diffs_ind[i]);
                    for (size_t k = 0; k < tgt.size(); ++k)
                        tgt[k] += ids[deg_tgt];

                    stmt.bind_int(1, src);
                    stmt.bind_int(2, r);
                    stmt.bind_str(3, myio::Serialize(tgt));
                    stmt.step_and_reset();
                }
            }
        }
    }

    /* ss_nd */
    {
        dbOut.drop_and_create_ss_nd(table_prefix);
        myio::Statement stmt(dbOut, "INSERT INTO " + table_prefix + "_ss_nd (src, r, tgt) VALUES (?1, ?2, ?3);");
        SS ss(gb, basis, basis_ss);
        ss.CacheNullDiffs(5);
        auto& nds = ss.GetND();

        for (size_t i = 0; i < nds.back().size(); ++i) {
            auto& nd = nds.back()[i];
            int src = ids[nd.deg] + (int)nd.index;
            auto& basis_ss_d = basis_ss.at(nd.deg);
            if (basis_ss_d.levels[nd.index] > 9800) {
                int r = kLevelMax - basis_ss_d.levels[nd.index];
                AdamsDeg deg_tgt = nd.deg + AdamsDeg(r, r - 1);
                int1d tgt;
                for (int j = nd.first; j < nd.first + nd.count; ++j)
                    tgt.push_back(j);
                for (size_t k = 0; k < tgt.size(); ++k)
                    tgt[k] += ids[deg_tgt];

                stmt.bind_int(1, src);
                stmt.bind_int(2, r);
                stmt.bind_str(3, myio::Serialize(tgt));
                stmt.step_and_reset();
            }
        }
    }

    /* ss_stable_levels */
    {
        dbOut.drop_and_create_ss_stable_levels(table_prefix);
        myio::Statement stmt(dbOut, "INSERT INTO " + table_prefix + "_ss_stable_levels (s, t, l) VALUES (?1, ?2, ?3);");
        SS ss(gb, basis, basis_ss);

        for (auto& [deg, basis_ss_d] : basis_ss) {
            stmt.bind_int(1, deg.s);
            stmt.bind_int(2, deg.t);
            int level = ss.GetFirstFixedLevelForPlot(deg);
            stmt.bind_int(3, level);
            stmt.step_and_reset();
        }
    }

    dbOut.end_transaction();
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