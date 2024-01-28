#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/myio.h"
#include "groebner_res_const.h"
#include "main.h"
#include "steenrod_sec.h"

using namespace steenrod;

namespace steenrod {
void SortMod2(MMod1d& data);
void SortMod4(MMilnorSec1d& data);
void MulMilnorMod4(MMilnorSec lhs, MMilnorSec rhs, MilnorSec& result_app, Milnor& tmp1, Milnor& tmp2);
}

class DbAdamsd2Map : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsd2Map(const std::string& filename) : Database(filename)
    {
        if (newFile_)
            SetVersion();
    }

    void SetVersion()
    {
        create_db_version(*this);
        Statement stmt(*this, "INSERT INTO version (id, name, value) VALUES (?1, ?2, ?3) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
        stmt.bind_and_step(0, std::string("version"), DB_ADAMS_VERSION);
        stmt.bind_and_step(1, std::string("change notes"), std::string(DB_VERSION_NOTES));
        stmt.bind_and_step(817812698, std::string("t_max"), -1);
    }

public:
    void create_tables(const std::string& table_prefix)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + " (id INTEGER PRIMARY KEY, d2 BLOB, d2_h TEXT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_time (s SMALLINT, t SMALLINT, time REAL, PRIMARY KEY (s, t));");
    }

    void save_time(const std::string& table_prefix, int s, int t, double time)
    {
        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_time (s, t, time) VALUES (?1, ?2, ?3);");
        stmt.bind_and_step(s, t, time);
    }

    int1d load_old_ids(std::string_view table_prefix) const
    {
        int1d result;
        Statement stmt(*this, fmt::format("SELECT DISTINCT id FROM {} ORDER BY id", table_prefix));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            result.push_back(id);
        }
        return result;
    }

    /* result[index] = image(v_index) in degree s */
    Mod1d load_d2(std::string_view table_prefix, int s) const
    {
        Mod1d result;
        Statement stmt(*this, fmt::format("SELECT id, d2 FROM {} WHERE (id>>19)={} ORDER BY id;", table_prefix, s));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            Mod d2;
            d2.data = stmt.column_blob_tpl<MMod>(1);
            int v = LocId(id).v;
            ut::get(result, v) = std::move(d2);
        }
        return result;
    }
};

Mod A_dd(const Mod& dg, const Mod2d& diffs, int s)
{
    Mod result;
    Milnor tmp, tmp1, tmp2;
    std::vector<MilnorSec> dd;
    //std::vector<Milnor> dd_check;  //////
    for (auto& m : dg.data) {
        auto& dv = diffs[size_t(s - 1)][m.v()];
        dd.clear();
        for (auto& mdv : dv.data) {
            auto& dv_mdv = diffs[size_t(s - 2)][mdv.v()];
            /* ut::get(dd, mdv_mdv.v()) += MilnorSec(mdv.m()) * MilnorSec(mdv_mdv.m()) */
            for (auto& mdv_mdv : dv_mdv.data)
                MulMilnorMod4(mdv.m(), mdv_mdv.m(), ut::get(dd, mdv_mdv.v()), tmp1, tmp2); 
        }
        for (size_t v = 0; v < dd.size(); ++v) {
            SortMod4(dd[v].data);
            for (auto& mdd : dd[v].data) {
                MulSec(m.m(), mdd, tmp, tmp1, tmp2);
                for (auto& m_tmp : tmp.data)
                    result.data.push_back(MMod(m_tmp, v));
            }
        }
    }
    SortMod2(result.data);
    return result;
}

/*
 * Ext^2(cw) --> H^*(cw)
 *
 *  F_s -----f-----> F_{s-2}
 *   |   \            |
 *   d       A        d
 *   |            \   |
 *   V                V
 *  F_{s-1} --f--> F_{s-3}
 */
void compute_d2(const std::string& cw, int t_trunc, int stem_trunc)
{
    std::string db_d2 = fmt::format("{}_Adams_d2.db", cw);
    std::string table_d2 = fmt::format("{}_Adams_d2", cw);
    DbAdamsd2Map dbD2(db_d2);
    int old_t_max_d2 = get_db_t_max(dbD2);

    std::string db_cw = cw + "_Adams_res.db";
    std::string table_cw = cw + "_Adams_res";
    myio::AssertFileExists(db_cw);

    dbD2.create_tables(table_d2);
    myio::Statement stmt_map(dbD2, fmt::format("INSERT OR IGNORE INTO {} (id, d2, d2_h) VALUES (?1, ?2, ?3);", table_d2)); /* (id, map, map_h) */
    myio::Statement stmt_t_max(dbD2, "INSERT INTO version (id, name, value) VALUES (817812698, \"t_max\", ?1) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
    myio::Statement stmt_time(dbD2, "INSERT INTO version (id, name, value) VALUES (1954841564, \"timestamp\", unixepoch()) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");

    DbAdamsResLoader dbRes(db_cw);
    int t_max_cw = get_db_t_max(dbRes);

    std::vector<std::pair<int, AdamsDegV2>> id_deg; /* pairs (id, deg) where `id` is the first id in deg */
    int2d vid_num;                                  /* vid_num[s][stem] is the number of generators in (<=stem, s) */
    Mod2d diffs;                                    /* diffs[s][i] is dv_i in filtration s */
    std::map<AdamsDegV2, size_t> num_diffs;         /* num_diffs[deg] is the number of differentials in deg */
    if (t_trunc > t_max_cw) {
        t_trunc = t_max_cw;
        fmt::print("t_max is truncated to {}\n", t_max_cw);
    }
    dbRes.load_generators(table_cw, id_deg, vid_num, diffs, num_diffs, t_trunc, stem_trunc);
    auto gb = AdamsResConst::load(dbRes, table_cw, t_trunc);

    /* Remove computed range */
    int1d ids_old = dbD2.load_old_ids(table_d2);
    ut::RemoveIf(id_deg, [&ids_old](const std::pair<int, AdamsDegV2>& p) { return ut::has(ids_old, p.first); });

    bench::Timer timer;
    timer.SuppressPrint();

    int t_prev = -1;
    AdamsDegV2 deg1_old(-1, -1);
    AdamsDegV2 deg2_old(-1, -1);
    for (const auto& [id, deg] : id_deg) {
        int s = deg.s;
        if (s < 2)
            continue;
        const auto& diffs_d = &diffs[s][LocId(id).v];
        const size_t diffs_d_size = num_diffs.at(deg);

        Mod1d f;
        f.resize(diffs_d_size);
        int2d fh;
        fh.resize(diffs_d_size);

        if (s == 2) { /* fill in the first map */
        }
        else {
            /* f_{s-1} is the map F_{s-1} -> F_{s-3} */
            Mod1d f_sm1 = dbD2.load_d2(table_d2, deg.s - 1);
            size_t vid_num_sm1 = deg.s > 0 ? (size_t)vid_num[size_t(deg.s - 1)][deg.stem()] : 0;
            f_sm1.resize(vid_num_sm1);

            /*# compute fd+A */
            Mod1d fd;
            fd.resize(diffs_d_size);
            ut::for_each_par128(diffs_d_size, [&fd, &diffs_d, &diffs, &f_sm1, s](size_t i) { fd[i] = subs(diffs_d[i], f_sm1) + A_dd(diffs_d[i], diffs, s); });

            //for (size_t i = 0; i < diffs_d_size; ++i)  ////
            //    fmt::print("fd[{}] = {}\n", i, fd[i].StrXi());
            //fmt::print("\n");

            /*# compute f */
            gb.DiffInvBatch(fd, f, size_t(deg.s - 3));

            
            //for (size_t i = 0; i < diffs_d_size; ++i)  ////
            //    fmt::print("f[{}] = {}\n", i, f[i].StrXi());
            //fmt::print("\n");

            /*# compute fh */
            for (size_t i = 0; i < diffs_d_size; ++i)
                fh[i] = HomToK(f[i]);
        }

        dbD2.begin_transaction();
        if (t_prev != deg.t) {
            if (t_prev != -1) {
                stmt_t_max.bind_and_step(t_prev);
                stmt_time.step_and_reset();
            }
            t_prev = deg.t;
        }
        /*# save products to database */
        for (size_t i = 0; i < diffs_d_size; ++i)
            stmt_map.bind_and_step(id + (int)i, f[i].data, myio::Serialize(fh[i]));

        double time = timer.Elapsed();
        timer.Reset();
        fmt::print("t={} s={} time={}\n", deg.t, deg.s, time);
        std::fflush(stdout);
        dbD2.save_time(table_d2, deg.s, deg.t, time);

        dbD2.end_transaction();
    }

    stmt_t_max.bind_and_step(std::max({t_prev, old_t_max_d2, t_trunc}));
    stmt_time.step_and_reset();
}

// TEST case
// const std::array<uint32_t, XI_MAX> R1 = {0, 2, 0, 0, 2};
// const std::array<uint32_t, XI_MAX> R2 = {1, 1, 4, 4};
// const std::array<uint32_t, XI_MAX> R3 = {0, 2, 1, 8};

int main_d2(int argc, char** argv, int& index, const char* desc)
{
    std::string cw;
    int t_max = 0;

    myio::CmdArg1d args = {{"cw", &cw}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    compute_d2(cw, t_max, DEG_MAX);
    return 0;
}