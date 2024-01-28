#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/linalg.h"
#include "algebras/utility.h"
#include "groebner_res_const.h"
#include "main.h"
#include <cstring>

int get_db_t_verified(const myio::Database& db)
{
    try {
        if (db.has_table("version")) {
            try {
                return db.get_int("select value from version where id=1121750147");
            }
            catch (MyException&) {
                return -1;
            }
        }
    }
    catch (MyException&) {
        return -3;
    }
    return -2;
}

class DbAdamsVerifyLoader : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsVerifyLoader(const std::string& filename) : Database(filename) {}

public:
    int1d load_ids(std::string_view table_prefix) const
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
    Mod1d load_map(std::string_view table_prefix, int s) const
    {
        Mod1d result;
        Statement stmt(*this, fmt::format("SELECT id, map FROM {} WHERE (id>>19)={} ORDER BY id;", table_prefix, s));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            Mod map_;
            map_.data = stmt.column_blob_tpl<MMod>(1);
            int v = LocId(id).v;
            ut::get(result, v) = std::move(map_);
        }
        return result;
    }

    /* vid_num[s][stem] is the number of generators in (<=stem, s) */
    Mod1d load_generators(const std::string& table_prefix, int s, int t_max) const
    {
        Mod1d result;
        Statement stmt(*this, fmt::format("SELECT diff FROM {}_generators WHERE t<={} AND s={} ORDER BY id;", table_prefix, t_max, s));
        while (stmt.step() == MYSQLITE_ROW) {
            Mod diff;
            diff.data = stmt.column_blob_tpl<MMod>(0);
            result.push_back(std::move(diff));
        }
        return result;
    }
};

/* cw1 --> cw2
 * Ext^fil(cw2) --> H^*(cw1)
 *
 *  F_s -----f-----> F_{s-fil}
 *   |                |
 *   d                d
 *   |                |
 *   V                V
 *  F_{s-1} --f--> F_{s-1-fil}
 */
void verify_map(const std::string& cw1, const std::string& cw2)
{
    std::string db_map = fmt::format("map_Adams_res_{}__{}.db", cw1, cw2);
    std::string table_map = fmt::format("map_Adams_res_{}__{}", cw1, cw2);
    myio::AssertFileExists(db_map);
    DbAdamsVerifyLoader dbMap(db_map);
    int t_max_map = get_db_t_max(dbMap);

    auto from = dbMap.get_str("select value from version where id=446174262");
    auto to = dbMap.get_str("select value from version where id=1713085477");

    std::string db_cw1 = from + "_Adams_res.db";
    std::string table_cw1 = from + "_Adams_res";
    std::string db_cw2 = to + "_Adams_res.db";
    std::string table_cw2 = to + "_Adams_res";
    myio::AssertFileExists(db_cw1);
    myio::AssertFileExists(db_cw2);

    int fil = 0;
    try { /* For compatibility */
        fil = dbMap.get_int("select value from version where id=651971502");
    }
    catch (MyException&) {
    }
    const int sus = dbMap.get_int("select value from version where id=1585932889");
    myio::Statement stmt_verify(dbMap, "INSERT INTO version (id, name, value) VALUES (1121750147, \"verified\", ?1) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");

    std::vector<std::pair<int, AdamsDegV2>> id_deg; /* pairs (id, deg) where `id` is the first id in deg */
    int2d vid_num;                                  /* vid_num[s][stem] is the number of generators in (<=stem, s) */
    std::map<AdamsDegV2, Mod1d> diffs_cw2;          /* diffs[deg] is the list of differentials of v in deg */
    int t_max_cw2;
    {
        DbAdamsResLoader dbResCw2(db_cw2);
        t_max_cw2 = get_db_t_max(dbResCw2);
        dbResCw2.load_generators(table_cw2, id_deg, vid_num, diffs_cw2, std::max(t_max_map + fil - sus, 0), 10000);
    }

    DbAdamsVerifyLoader dbResCw1(db_cw1);
    dbResCw1.execute_cmd(fmt::format("CREATE INDEX IF NOT EXISTS index_s ON {}_generators (s)", table_cw1));
    dbResCw1.execute_cmd(fmt::format("CREATE INDEX IF NOT EXISTS index_t ON {}_generators (t)", table_cw1));

    /* Remove computed range */
    int1d ids = dbMap.load_ids(table_map);

    int t_prev_cw2 = -1;
    int t_verified = get_db_t_verified(dbMap);
    for (const auto& [id, deg] : id_deg) {
        if (deg.s <= fil || deg.t - fil + sus <= t_verified)
            continue;
        const auto& diffs_cw2_d = diffs_cw2.at(deg);
        const size_t diffs_cw2_d_size = diffs_cw2_d.size();

        /* f_s: F_s -> F_{s-fil}
         * f_{s-1}: F_{s-1} -> F_{s-1-fil}
         */
        Mod1d f_sm1 = dbMap.load_map(table_map, deg.s - 1);
        size_t vid_num_sm1 = deg.s > 0 ? (size_t)vid_num[size_t(deg.s - 1)][deg.stem()] : 0;
        f_sm1.resize(vid_num_sm1);

        /*# compute fd */
        Mod1d fd;
        fd.resize(diffs_cw2_d_size);
        ut::for_each_par128(diffs_cw2_d_size, [&fd, &diffs_cw2_d, &f_sm1](size_t i) { fd[i] = subs(diffs_cw2_d[i], f_sm1); });

        /*# compute df */
        Mod1d f_s = dbMap.load_map(table_map, deg.s);
        f_s.resize(vid_num[deg.s][deg.stem()]);
        Mod1d f;
        int v = LocId(id).v;
        for (size_t i = 0; i < diffs_cw2_d_size; ++i)
            f.push_back(f_s[v + i]);

        Mod1d df;
        df.resize(diffs_cw2_d_size);
        auto diffs_cw1 = dbResCw1.load_generators(table_cw1, deg.s - fil, deg.t - fil + sus);
        ut::for_each_par128(diffs_cw2_d_size, [&df, &diffs_cw1, &f](size_t i) { df[i] = subs(f[i], diffs_cw1); });

        if (fd != df) {
            fmt::print("Error! cw1={} cw2={} (s, t)=({}, {})\n", cw1, cw2, deg.s, deg.t);
            dbMap.begin_transaction();
            stmt_verify.bind_and_step(deg.t - fil + sus + 10000);
            dbMap.end_transaction();
            return;
        }
        fmt::print("t={} s={}\n", deg.t - fil + sus, deg.s - fil);
        std::fflush(stdout);
        if (t_prev_cw2 != deg.t) {
            if (t_prev_cw2 != -1 && t_prev_cw2 - fil + sus >= 0) {
                dbMap.begin_transaction();
                stmt_verify.bind_and_step(t_prev_cw2 - fil + sus);
                dbMap.end_transaction();
            }
            t_prev_cw2 = deg.t;
        }
        diffs_cw2.erase(deg);
    }

    stmt_verify.bind_and_step(std::max(t_prev_cw2 - fil + sus, t_verified));
}

int main_verify_map(int argc, char** argv, int& index, const char* desc)
{
    std::string cw1, cw2;

    myio::CmdArg1d args = {{"cw1", &cw1}, {"cw2", &cw2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    verify_map(cw1, cw2);
    return 0;
}