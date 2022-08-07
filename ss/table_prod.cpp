#include "algebras/linalg.h"
#include "main.h"

using namespace alg;

class MyDB : public SSDB
{
    using Statement = myio::Statement;

public:
    MyDB() = default;
    explicit MyDB(const std::string& filename) : SSDB(filename) {}

    void create_basis_products(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis_products (id1 INTEGER, id2 INTEGER, prod TEXT, PRIMARY KEY(id1, id2));");
    }

    void create_basis_products_and_delete(const std::string& table_prefix) const
    {
        create_basis_products(table_prefix);
        delete_from(table_prefix + "_basis_products");
    }

    void create_ss_products(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_products (id1 INTEGER, id2 INTEGER, prod TEXT, PRIMARY KEY(id1, id2));");
    }

    void create_ss_products_and_delete(const std::string& table_prefix) const
    {
        create_ss_products(table_prefix);
        delete_from(table_prefix + "_ss_products");
    }

    void create_ss_diff(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_diffs (src INTEGER PRIMARY KEY, r SMALLINT, tgt TEXT);");
    }

    void create_ss_diff_and_delete(const std::string& table_prefix) const
    {
        create_ss_diff(table_prefix);
        delete_from(table_prefix + "_ss_diffs");
    }

    void create_ss_nd(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_nd (src INTEGER PRIMARY KEY, r SMALLINT, tgt TEXT);");
    }

    void create_ss_nd_and_delete(const std::string& table_prefix) const
    {
        create_ss_nd(table_prefix);
        delete_from(table_prefix + "_ss_nd");
    }

    void create_ss_stable_levels(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_stable_levels (s SMALLINT, t SMALLINT, l SMALLINT, PRIMARY KEY(s, t));");
    }

    void create_ss_stable_levels_and_delete(const std::string& table_prefix) const
    {
        create_ss_stable_levels(table_prefix);
        delete_from(table_prefix + "_ss_stable_levels");
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
    std::string db_filename = "AdamsE2Export_t220.db";
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate the basis_prod table\n";
        std::cout << "Usage:\n  ss basis_prod [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "Version:\n  1.1 (2022-7-16)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_prefix", table_prefix))
        return index;

    MyDB db(db_filename);
    Poly1d polys = db.load_gb(table_prefix, DEG_MAX);
    auto basis = db.load_basis(table_prefix);
    auto ids = db.load_indices(table_prefix);
    Mon1d basis_v2;
    AdamsDeg1d deg_basis;
    db.load_basis_v2(table_prefix, basis_v2, deg_basis);

    int t_max = basis.rbegin()->first.t;
    Groebner gb(t_max, {}, polys);

    db.begin_transaction();
    db.create_basis_products_and_delete(table_prefix);
    myio::Statement stmt(db, "INSERT INTO " + table_prefix + "_basis_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
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
    db.end_transaction();
    return 0;
}

int main_plot(int argc, char** argv, int index)
{
    std::string db_filename = "AdamsE2Export_t220.db";
    std::string table_prefix = "AdamsE2";

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

    MyDB db(db_filename);
    Poly1d polys = db.load_gb(table_prefix, DEG_MAX);
    auto basis = db.load_basis(table_prefix);
    const Staircases basis_ss = db.load_basis_ss(table_prefix);
    auto ids = db.load_indices(table_prefix);
    int2d basis_ss_v2;
    AdamsDeg1d deg_basis;
    db.load_basis_ss_v2(table_prefix, basis_ss_v2, deg_basis);

    int t_max = deg_basis.back().t;
    Groebner gb(t_max, {}, polys);

    db.begin_transaction();

    /* ss_products */
    {
        db.create_ss_products_and_delete(table_prefix);
        myio::Statement stmt(db, "INSERT INTO " + table_prefix + "_ss_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
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
        db.create_ss_diff_and_delete(table_prefix);
        myio::Statement stmt(db, "INSERT INTO " + table_prefix + "_ss_diffs (src, r, tgt) VALUES (?1, ?2, ?3);");
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
        db.create_ss_nd_and_delete(table_prefix);
        myio::Statement stmt(db, "INSERT INTO " + table_prefix + "_ss_nd (src, r, tgt) VALUES (?1, ?2, ?3);");
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
        db.create_ss_stable_levels_and_delete(table_prefix);
        myio::Statement stmt(db, "INSERT INTO " + table_prefix + "_ss_stable_levels (s, t, l) VALUES (?1, ?2, ?3);");
        SS ss(gb, basis, basis_ss);

        for (auto& [deg, basis_ss_d] : basis_ss) {
            stmt.bind_int(1, deg.s);
            stmt.bind_int(2, deg.t);
            int level = ss.GetFirstFixedLevelForPlot(deg);
            stmt.bind_int(3, level);
            stmt.step_and_reset();
        }
    }

    db.end_transaction();
    return 0;
}