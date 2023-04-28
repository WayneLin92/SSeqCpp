#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/linalg.h"
#include "algebras/utility.h"
#include "groebner_res_const.h"
#include "main.h"
#include <cstring>

int1d HomToK(const Mod& x)
{
    int1d result;
    for (MMod m : x.data)
        if (m.deg_m() == 0)
            result.push_back((int)m.v());
    return result;
}

int1d HomToMSq(const Mod& x, int t_cell);

class DbAdamsResProd : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsResProd(const std::string& filename) : Database(filename) {}

public:
    void create_products(const std::string& table_prefix)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, indecomposable TINYINT, s SMALLINT, t SMALLINT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_products (id INTEGER, id_ind INTEGER, prod BLOB, prod_h BLOB, prod_h_glo TEXT, PRIMARY KEY (id, id_ind));");
    }

    void create_time(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_products_time (s SMALLINT, t SMALLINT, time REAL, PRIMARY KEY (s, t));");
    }

    void save_time(const std::string& table_prefix, int s, int t, double time)
    {
        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_products_time (s, t, time) VALUES (?1, ?2, ?3);");

        stmt.bind_int(1, s);
        stmt.bind_int(2, t);
        stmt.bind_double(3, time);
        stmt.step_and_reset();
    }

    int1d load_old_ids(std::string_view table_prefix) const
    {
        int1d result;
        Statement stmt(*this, fmt::format("SELECT DISTINCT id FROM {}_products ORDER BY id", table_prefix));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            result.push_back(id);
        }
        return result;
    }

    /* result[g][index] = image(v_index) in degree s */
    std::map<int, Mod1d> load_products(std::string_view table_prefix, int s, const LocId1d& glo2loc, const int1d& gs_exclude) const
    {
        std::map<int, Mod1d> result;
        Statement stmt(*this, fmt::format("SELECT id, id_ind, prod FROM {}_products LEFT JOIN {}_generators USING(id) WHERE s={} ORDER BY id;", table_prefix, table_prefix, s));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            int g = stmt.column_int(1);
            if (!ut::has(gs_exclude, g)) {
                Mod prod;
                prod.data = stmt.column_blob_tpl<MMod>(2);

                int index = glo2loc[id].v;
                if (result[g].size() <= index)
                    result.at(g).resize(size_t(index + 1));
                result[g][index] = std::move(prod);
            }
        }
        return result;
    }
};

AdamsResConst AdamsResConst::load(const DbAdamsResLoader& db, const std::string& table, int t_trunc)
{
    DataMResConst2d data = db.load_data(table, t_trunc);
    int2d basis_degrees = db.load_basis_degrees(table, t_trunc);
    return AdamsResConst(std::move(data), std::move(basis_degrees));
}

/*  F_s ----f----> F_{s-g}
 *   |                |
 *   d                d
 *   |                |
 *   V                V
 *  F_{s-1} --f--> F_{s-1-g}
 */
void compute_products_by_t(int t_trunc, int stem_trunc, const std::string& db_res, const std::string& table_res, const std::string& db_out, const std::string& table_out)
{
    DbAdamsResProd dbProd(db_out);
    dbProd.create_products(table_out);
    dbProd.create_time(table_out);
    myio::Statement stmt_gen(dbProd, fmt::format("INSERT INTO {}_generators (id, indecomposable, s, t) values (?1, ?2, ?3, ?4);", table_out));   /* (id, is_indecomposable, s, t) */
    myio::Statement stmt_set_ind(dbProd, fmt::format("UPDATE {}_generators SET indecomposable=1 WHERE id=?1 and indecomposable=0;", table_out)); /* id */
    myio::Statement stmt_prod(dbProd, fmt::format("INSERT INTO {}_products (id, id_ind, prod, prod_h) VALUES (?1, ?2, ?3, ?4);", table_out));    /* (id, id_ind, prod, prod_h) */

    DbAdamsResLoader dbRes(db_res);
    auto gb = AdamsResConst::load(dbRes, table_res, t_trunc);
    std::map<int, AdamsDegV2> id_deg;  /* pairs (id, deg) where `id` is the first id in deg */
    int2d vid_num;                     /* vid_num[s][stem] is the number of generators in (<=stem, s) */
    std::map<AdamsDegV2, Mod1d> diffs; /* diffs[deg] is the list of differentials of v in deg */
    dbRes.load_generators(table_res, id_deg, vid_num, diffs, t_trunc, stem_trunc);
    int2d loc2glo; /* loc2glo[s][vid]=id */
    LocId1d glo2loc; /* loc2glo[id]=(s, vid) */
    dbRes.load_id_converter(table_res, loc2glo, glo2loc);
    int1d gs_hopf, t_gs_hopf;
    if (table_res == "S0_Adams_res" || table_res == "tmf_Adams_res") {
        gs_hopf = dbRes.get_column_int(fmt::format("{}_generators", table_res), "id", "WHERE s=1 ORDER BY id");
        t_gs_hopf = dbRes.get_column_int(fmt::format("{}_generators", table_res), "t", "WHERE s=1 ORDER BY id");
    }
    dbRes.disconnect();

    const Mod one = MMod(MMilnor(), 0);
    const int1d one_h = {0};

    /* Remove computed range */
    {
        int1d ids_old = dbProd.load_old_ids(table_out);
        for (auto p = id_deg.begin(); p != id_deg.end();) {
            if (ut::has(ids_old, p->first)) {
                diffs.erase(p->second);
                id_deg.erase(p++);
            }
            else
                ++p;
        }
    }

    bench::Timer timer;
    timer.SuppressPrint();

    for (const auto& [id, deg] : id_deg) {
        const auto& diffs_d = diffs.at(deg);
        const size_t diffs_d_size = diffs_d.size();

        dbProd.begin_transaction();

        /* save generators to database */
        for (size_t i = 0; i < diffs_d_size; ++i) {
            stmt_gen.bind(id + (int)i, 0, deg.s, deg.t);
            stmt_gen.step_and_reset();
        }

        if (deg.t == 0) {
            diffs.erase(deg);
            continue;
        }

        /* f_{s-1}[g] is the map F_{s-1} -> F_{s-1-deg(g)} dual to the multiplication of g */
        std::map<int, Mod1d> f_sm1 = dbProd.load_products(table_out, deg.s - 1, glo2loc, gs_hopf);
        int1d gs = ut::get_keys(f_sm1); /* indecomposables id's */

        /*# compute fd */
        std::map<int, Mod1d> fd;
        size_t vid_num_sm1 = deg.s > 0 ? (size_t)vid_num[size_t(deg.s - 1)][deg.stem()] : 0;
        for (auto& [g, f_sm1_g] : f_sm1) {
            f_sm1_g.resize(vid_num_sm1);
            fd[g].resize(diffs_d_size);
        }
        ut::for_each_par128(diffs_d_size * gs.size(), [&gs, &fd, &diffs_d, &f_sm1, diffs_d_size](size_t i) {
            int g = gs[i / diffs_d_size];
            size_t j = i % diffs_d_size;
            fd.at(g)[j] = subs(diffs_d[j], f_sm1.at(g));
        });

        /*# compute f */
        std::map<int, Mod1d> f;
        int1d s1;
        for (auto& [g, _] : fd) {
            s1.push_back(deg.s - 1 - glo2loc.at(g).s);
            f[g].resize(diffs_d_size);
        }
        ut::for_each_par128(gs.size(), [&gs, &gb, &fd, &f, &s1](size_t i) { gb.DiffInvBatch(fd[gs[i]], f[gs[i]], s1[i]); });

        /*# compute fh */
        std::map<int, int2d> fh;
        for (auto& [g, f_g] : f)
            for (size_t i = 0; i < diffs_d_size; ++i)
                fh[g].push_back(HomToK(f_g[i]));
        if (deg.s > 1) {
            for (size_t i_g = 0; i_g < gs_hopf.size(); ++i_g) {
                int g = gs_hopf[i_g];
                int t = t_gs_hopf[i_g];
                for (size_t i = 0; i < diffs_d_size; ++i)
                    fh[g].push_back(HomToMSq(diffs_d[i], t));
            }
        }

        /*# save products to database */
        for (auto& [g, f_g] : f) {
            for (size_t i = 0; i < diffs_d_size; ++i) {
                if (f_g[i]) {
                    stmt_prod.bind(id + (int)i, g, f_g[i].data, fh.at(g)[i]);
                    stmt_prod.step_and_reset();
                }
            }
        }
        if (deg.s > 1) {
            for (size_t i_g = 0; i_g < gs_hopf.size(); ++i_g) {
                int g = gs_hopf[i_g];
                int t = t_gs_hopf[i_g];
                for (size_t i = 0; i < diffs_d_size; ++i) {
                    if (!fh.at(g)[i].empty()) {
                        stmt_prod.bind(id + (int)i, g, myio::SQL_NULL(), fh.at(g)[i]);
                        stmt_prod.step_and_reset();
                    }
                }
            }
        }

        /*# find indecomposables */
        int2d fx;
        for (const auto& [_, fh_g] : fh) {
            size_t offset = fx.size();
            for (size_t i = 0; i < fh_g.size(); ++i) {
                for (int k : fh_g[i]) {
                    if (fx.size() <= offset + (size_t)k)
                        fx.resize(offset + (size_t)k + 1);
                    fx[offset + (size_t)k].push_back((int)i);
                }
            }
        }
        int1d lead_image = lina::GetLeads(lina::GetSpace(fx));
        int1d indices = lina::AddVectors(ut::int_range((int)diffs_d_size), lead_image);

        /*# mark indecomposables in database */
        for (int i : indices) {
            stmt_set_ind.bind_int(1, id + i);
            stmt_set_ind.step_and_reset();
        }

        /*# indecomposable comultiply with itself */
        for (int i : indices) {
            stmt_prod.bind(id + (int)i, id + (int)i, one.data, one_h);
            stmt_prod.step_and_reset();
        }

        double time = timer.Elapsed();
        timer.Reset();
        fmt::print("t={} s={} time={}\n", deg.t, deg.s, time);
        dbProd.save_time(table_out, deg.s, deg.t, time);

        dbProd.end_transaction();
        diffs.erase(deg);
    }
}

/*  F_s ----f----> F_{s-g}
 *   |                |
 *   d                d
 *   |                |
 *   V                V
 *  F_{s-1} --f--> F_{s-1-g}
 */
void compute_mod_products_by_t(int t_trunc, int stem_trunc, const std::string& db_mod, const std::string& table_mod, const std::string& db_ring, const std::string& table_ring, const std::string& db_out, const std::string& table_out)
{
    DbAdamsResProd dbProd(db_out);
    dbProd.create_products(table_out);
    dbProd.create_time(table_out);
    myio::Statement stmt_gen(dbProd, fmt::format("INSERT INTO {}_generators (id, indecomposable, s, t) values (?1, ?2, ?3, ?4);", table_out));   /* (id, is_indecomposable, s, t) */
    myio::Statement stmt_set_ind(dbProd, fmt::format("UPDATE {}_generators SET indecomposable=1 WHERE id=?1 and indecomposable=0;", table_out)); /* id */
    myio::Statement stmt_prod(dbProd, fmt::format("INSERT INTO {}_products (id, id_ind, prod, prod_h) VALUES (?1, ?2, ?3, ?4);", table_out));    /* (id, id_ind, prod, prod_h) */
    int1d ids_old = dbProd.load_old_ids(table_out);

    DbAdamsResLoader dbResRing(db_ring);
    auto gbRing = AdamsResConst::load(dbResRing, table_ring, t_trunc);
    dbResRing.disconnect();

    std::map<int, AdamsDegV2> id_deg;  /* pairs (id, deg) where `id` is the first id in deg */
    int2d vid_num;                     /* vid_num[s][stem] is the number of generators in (<=stem, s) */
    std::map<AdamsDegV2, Mod1d> diffs; /* diffs[deg] is the list of differentials of v in deg */
    int2d loc2glo;               /* loc2glo[s][vid]=id */
    LocId1d glo2loc;               /* loc2glo[id]=(s, vid) */
    {
        DbAdamsResLoader dbRes(db_mod);
        dbRes.load_generators(table_mod, id_deg, vid_num, diffs, t_trunc, stem_trunc);
        dbRes.load_id_converter(table_mod, loc2glo, glo2loc);
    }

    const Mod one = MMod(MMilnor(), 0);
    const int1d one_h = {0};

    /* Remove computed range */
    {
        int1d ids_old = dbProd.load_old_ids(table_out);
        for (auto p = id_deg.begin(); p != id_deg.end();) {
            if (ut::has(ids_old, p->first)) {
                diffs.erase(p->second);
                id_deg.erase(p++);
            }
            else
                ++p;
        }
    }

    bench::Timer timer;
    timer.SuppressPrint();

    for (const auto& [id, deg] : id_deg) {
        const auto& diffs_d = diffs.at(deg);
        const size_t diffs_d_size = diffs_d.size();

        dbProd.begin_transaction();

        /* save generators to database */
        for (size_t i = 0; i < diffs_d_size; ++i) {
            stmt_gen.bind(id + (int)i, 0, deg.s, deg.t);
            stmt_gen.step_and_reset();
        }

        /* f_{s-1}[g] is the map F_{s-1} -> F_{s-1-deg(g)} dual to the multiplication of g */
        const int1d empty;
        std::map<int, Mod1d> f_sm1 = dbProd.load_products(table_out, deg.s - 1, glo2loc, empty);
        int1d gs = ut::get_keys(f_sm1); /* indecomposables id's */

        /*# compute fd */
        std::map<int, Mod1d> fd;
        size_t vid_num_sm1 = deg.s > 0 ? (size_t)vid_num[size_t(deg.s - 1)][deg.stem()] : 0;
        for (auto& [g, f_sm1_g] : f_sm1) {
            f_sm1_g.resize(vid_num_sm1);
            fd[g].resize(diffs_d_size);
        }
        ut::for_each_par128(diffs_d_size * gs.size(), [&gs, &fd, &diffs_d, &f_sm1, diffs_d_size](size_t i) {
            int g = gs[i / diffs_d_size];
            size_t j = i % diffs_d_size;
            fd.at(g)[j] = subs(diffs_d[j], f_sm1.at(g));
        });

        /*# compute f */
        std::map<int, Mod1d> f;
        int1d s1;
        for (auto& [g, _] : fd) {
            s1.push_back(deg.s - 1 - glo2loc.at(g).s);
            f[g].resize(diffs_d_size);
        }
        ut::for_each_par128(gs.size(), [&gs, &gbRing, &fd, &f, &s1](size_t i) { gbRing.DiffInvBatch(fd[gs[i]], f[gs[i]], s1[i]); });

        /*# compute fh */
        std::map<int, int2d> fh;
        for (auto& [g, f_g] : f)
            for (size_t i = 0; i < diffs_d_size; ++i)
                fh[g].push_back(HomToK(f_g[i]));

        /*# save products to database */
        for (auto& [g, f_g] : f) {
            for (size_t i = 0; i < diffs_d_size; ++i) {
                if (f_g[i]) {
                    stmt_prod.bind(id + (int)i, g, f_g[i].data, fh.at(g)[i]);
                    stmt_prod.step_and_reset();
                }
            }
        }

        /*# find indecomposables */
        int2d fx;
        for (const auto& [_, fh_g] : fh) {
            size_t offset = fx.size();
            for (size_t i = 0; i < fh_g.size(); ++i) {
                for (int k : fh_g[i]) {
                    if (fx.size() <= offset + (size_t)k)
                        fx.resize(offset + (size_t)k + 1);
                    fx[offset + (size_t)k].push_back((int)i);
                }
            }
        }
        int1d lead_image = lina::GetLeads(lina::GetSpace(fx));
        int1d indices = lina::AddVectors(ut::int_range((int)diffs_d_size), lead_image);

        /*# mark indecomposables in database */
        for (int i : indices) {
            stmt_set_ind.bind_int(1, id + i);
            stmt_set_ind.step_and_reset();
        }

        /*# indecomposable comultiply with itself */
        for (int i : indices) {
            stmt_prod.bind(id + (int)i, id + (int)i, one.data, one_h);
            stmt_prod.step_and_reset();
        }

        double time = timer.Elapsed();
        timer.Reset();
        fmt::print("t={} s={} time={}\n", deg.t, deg.s, time);
        dbProd.save_time(table_out, deg.s, deg.t, time);

        dbProd.end_transaction();
        diffs.erase(deg);
    }
}

void compute_products_with_hi(const std::string& db_res_S0, const std::string& db_res_mod, const std::string& table_res_mod, const std::string& db_out)
{
    int1d ids_h;
    {
        DbAdamsResLoader dbResS0(db_res_S0);
        ids_h = dbResS0.get_column_int("S0_Adams_res_generators", "id", "WHERE s=1 ORDER BY id");
    }

    std::map<int, AdamsDegV2> id_deg;
    int2d loc2glo;
    std::map<AdamsDegV2, Mod1d> diffs;
    {
        DbAdamsResLoader dbResMod(db_res_mod);
        int2d vid_num;
        dbResMod.load_generators(table_res_mod, id_deg, vid_num, diffs, 500, 500);
        LocId1d glo2loc;
        dbResMod.load_id_converter(table_res_mod, loc2glo, glo2loc);
    }

    DbAdamsResProd dbProd(db_out);
    dbProd.create_products(table_res_mod);
    int id_start = dbProd.get_int("SELECT COALESCE(MAX(id), -1) FROM " + table_res_mod + "_products") + 1;

    myio::Statement stmt_prod(dbProd, "INSERT INTO " + table_res_mod + "_products (id, id_ind, prod_h, prod_h_glo) VALUES (?1, ?2, ?3, ?4);");

    for (const auto& [id, deg] : id_deg) {
        if (id < id_start) {
            diffs.erase(deg);
            continue;
        }

        auto& diffs_d = diffs.at(deg);

        dbProd.begin_transaction();

        /* Save the hj products */
        for (size_t i = 0; i < diffs_d.size(); ++i) {
            for (size_t j = 0; j < ids_h.size(); ++j) {
                int1d prod_hi = HomToMSq(diffs_d[i], 1 << j);
                if (!prod_hi.empty()) {
                    int1d prod_hi_glo;
                    for (int k : prod_hi)
                        prod_hi_glo.push_back(loc2glo[size_t(deg.s - 1)][k]);

                    stmt_prod.bind_int(1, id + (int)i);
                    stmt_prod.bind_int(2, ids_h[j]);
                    stmt_prod.bind_blob(3, prod_hi);
                    stmt_prod.bind_str(4, myio::Serialize(prod_hi_glo));
                    stmt_prod.step_and_reset();
                }
            }
        }

        dbProd.end_transaction();
        diffs.erase(deg);
    }
}

/* Compute the product with hi */
int main_prod_hi(int argc, char** argv, int index)
{
    std::string complex = "S0";
    std::string db_S0 = "S0_Adams_res.db";
    std::string db_mod = "<cw>_Adams_res.db";
    std::string db_out = "<cw>_Adams_res_prod.db";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  Adams prod_hi <cw> [db_S0] [db_mod] [db_out]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  db_S0 = " << db_S0 << "\n";
        std::cout << "  db_mod = " << db_mod << "\n";
        std::cout << "  db_out = " << db_out << "\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "complex", complex))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_mod", db_mod))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_out", db_out))
        return index;

    if (db_mod == "<cw>_Adams_res.db")
        db_mod = complex + "_Adams_res.db";
    std::string table_mod = complex + "_Adams_res";
    if (db_out == "<cw>_Adams_res_prod.db")
        db_out = complex + "_Adams_res_prod.db";

    compute_products_with_hi(db_S0, db_mod, table_mod, db_out);

    return 0;
}

int main_prod(int argc, char** argv, int index)
{
    std::string cw = "S0";
    int t_max = 152, stem_max = 100;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Usage:\n  Adams prod <cw> [t_max] [stem_max]\n\n");

        fmt::print("Default values:\n");
        fmt::print("  t_max = {}\n", t_max);
        fmt::print("  stem_max = {}\n", stem_max);

        fmt::print("{}\n", VERSION);
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cw", cw))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return index;

    std::string db_res = cw + "_Adams_res.db";
    std::string table_res = cw + "_Adams_res";
    std::string db_out = cw + "_Adams_res_prod.db";
    std::string table_out = cw + "_Adams_res";

    myio::AssertFileExists(db_res);
    compute_products_by_t(t_max, stem_max, db_res, table_res, db_out, table_out);

    return 0;
}

int main_prod_mod(int argc, char** argv, int index)
{
    std::string cw;
    std::string ring;
    int t_max = 152, stem_max = 100;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Usage:\n  Adams prod_mod <cw> <ring> [t_max] [stem_max]\n\n");

        fmt::print("Default values:\n");
        fmt::print("  t_max = {}\n", t_max);
        fmt::print("  stem_max = {}\n", stem_max);

        fmt::print("{}\n", VERSION);
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cw", cw))
        return index;
    if (myio::load_arg(argc, argv, ++index, "ring", ring))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return index;

    std::string db_mod = cw + "_Adams_res.db";
    std::string table_mod = cw + "_Adams_res";
    std::string db_ring = ring + "_Adams_res.db";
    std::string table_ring = ring + "_Adams_res";
    std::string db_out = cw + "_Adams_res_prod.db";
    std::string table_out = cw + "_Adams_res";

    myio::AssertFileExists(db_mod);
    myio::AssertFileExists(db_ring);
    compute_mod_products_by_t(t_max, stem_max, db_mod, table_mod, db_ring, table_ring, db_out, table_out);

    return 0;
}
