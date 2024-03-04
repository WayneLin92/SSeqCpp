#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/groebner_steenrod.h"
#include "algebras/myio.h"
#include "algebras/steenrod.h"
#include "groebner_res_const.h"
#include "main.h"

using namespace steenrod;

namespace steenrod {
void SortMod2(MMod1d& data);
void MulMilnor(const std::array<uint32_t, XI_MAX>& R, const std::array<uint32_t, XI_MAX>& S, MMilnor1d& result);

/* Compute the contraction C(\xi_i^{2^k}, Sq(R)) inplace */
int Contr(uint32_t i, uint32_t k, std::array<uint32_t, XI_MAX>& R)
{
    if (i == 0)
        return 1;
    --i;
    if (R[i] >> k) {
        R[i] -= 1 << k;
        return 1;
    }
    return 0;
}

/* Compute the contraction C(\xi_i^{2^k}, Sq(R)) */
int Contr(uint32_t i, uint32_t k, const std::array<uint32_t, XI_MAX>& R, std::array<uint32_t, XI_MAX>& result)
{
    result = R;
    return Contr(i, k, result);
}

/* Compute the contraction C(\xi_i^{2^k}\xi_j^{2^l}, m)
 */
int Contr(uint32_t i, uint32_t k, uint32_t j, uint32_t l, const std::array<uint32_t, XI_MAX>& R, std::array<uint32_t, XI_MAX>& result)
{
    if (Contr(i, k, R, result))
        return Contr(j, l, result);
    return 0;
}

/**
 * Sort the sequence and each time remove a pair of identical elements
 */
void SortMod4(std::vector<uint64_t>& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i) {
        uint64_t c = data[i];
        if ((data[i] | 3) == (data[i + 1] | 3)) {
            size_t j = i + 1;
            c += data[j];
            data[j] = MMILNOR_NULL;
            for (++j; j < data.size() && (data[i] | 3) == (data[j] | 3); ++j) {
                c += data[j];
                data[j] = MMILNOR_NULL;
            }
            if ((c & 3) == 0)
                data[i] = MMILNOR_NULL;
            else
                data[i] = (data[i] & (~(uint64_t)3)) | (c & 3);
            i = j - 1;
        }
    }
    ut::RemoveIf(data, [](uint64_t m) { return m == MMILNOR_NULL; });
}

/* Binomial coefficient modulo 4 */
uint32_t BinomMod4(uint32_t n, uint32_t m)
{
    if ((n & m) == m) {
        int count3 = 0;
        while (n) {
            if (n % 4 == 3 && (m % 4 == 1 || m % 4 == 2))
                ++count3;
            n >>= 1;
            m >>= 1;
        }
        return count3 % 2 ? 3 : 1;
    }
    if (n > m && ut::popcount(m) + ut::popcount(n - m) - ut::popcount(n) == 1)
        return 2;
    return 0;
}

void MulMilnorMod4(const std::array<uint32_t, XI_MAX>& R, const std::array<uint32_t, XI_MAX>& S, std::vector<uint64_t>& result_app, uint64_t v_raw_shifted)
{
    constexpr size_t N = XI_MAX_MULT;
    constexpr size_t N1 = N + 1;

    /* XR[i,j], XS[i,j] controls the upper bound of X[i,j] when we consider R(X)=R and S(X)=S.
     * XT[i,j] is the partial sum of diagonals of X[i,j].
     */
    std::array<uint32_t, N1 * N1 - N> X, XT;  // TODO: use less memory
    std::array<uint32_t, N1 * N1 - N - N1> Xb, XR, XS;
    XT[N * N1] = S[N - 1] + R[N - 1];
    Xb[N * N1 - N1] = BinomMod4(XT[N * N1], R[N - 1]);
    if (Xb[N * N1 - N1] == 0)
        return;
    X[N] = R[N - 1];

    for (size_t i = 1; i <= N - 1; ++i)
        XR[(N - i) * N1 + i - N1] = R[i - 1];
    for (size_t i = 1; i <= N - 1; ++i)
        XS[i * N1 + (N - i) - N1] = S[i - 1];

    size_t i = N - 1, j = 1;
    bool decrease = false;
    while (true) {
        bool move_right = false;
        if (j) {
            const size_t index = i * N1 + j;
            const size_t index1 = index - N1;
            const size_t index_up = (i - 1) * N1 + j;
            const size_t index_up_right = (i - 1) * N1 + j + 1;
            const size_t index_left = i * N1 + j - 1;
            const size_t index_prev_b = i + j == N ? (i + 1) * N1 - N1 : i * N1 + j + 1 - N1;
            const size_t index_bottom_left = (i + 1) * N1 + j - 1;

            if (decrease)
                --X[index];
            else
                X[index] = std::min(XR[index - N1] >> i, XS[index - N1]);
            /* Row 1 is special because X[index_up_right] is determined before X[index] is determined */
            XT[index] = XT[(i > 1 || j == N - 1) ? index_bottom_left : index_up_right] + X[index];

            while (X[index]) {
                uint32_t b = (Xb[index_prev_b] * BinomMod4(XT[index], X[index])) % 4;
                if (b) {
                    Xb[index1] = b;
                    break;
                }
                --X[index];
                --XT[index];
            }
            if (X[index] == 0)
                Xb[index1] = Xb[index_prev_b];
            decrease = false;

            if (i == 1) {
                X[index_up] = XR[index - N1] - (X[index] << 1);
                if (j > 1) {
                    XT[index_up] = XT[(i + 1) * N1 + j - 2] + X[index_up];
                    Xb[index1] = (Xb[index1] * BinomMod4(XT[index_up], X[index_up])) % 4;
                    if (Xb[index1] == 0) {
                        if (X[index])
                            decrease = true;
                        else {
                            // fmt::print("X=\n{}\nXb=\n{}\n", X, Xb);
                            move_right = true;
                        }
                    }
                    else {
                        XS[index_left - N1] = XS[index - N1] - X[index];
                        --j;
                    }
                }
                else {
                    XS[index_left - N1] = XS[index - N1] - X[index];
                    --j;
                }
            }
            else {
                XS[index_left - N1] = XS[index - N1] - X[index];
                XR[index_up - N1] = XR[index - N1] - (X[index] << i);
                --j;
            }
        }
        else {
            if (i == 1) {
                XT[N + 1] = XS[N1 - N1] + X[1];
                uint32_t b = (Xb[N1 + 1 - N1] * BinomMod4(XT[N + 1], X[1])) % 4;
                // fmt::print("push: X=\n{}\nXS=\n{}\nXT=\n{}\nXb=\n{}\n", X, XS, XT, Xb);
                if (b) {
                    result_app.push_back((MMilnor::Xi(XT.data() + N + 1).data() << 2) | b | v_raw_shifted);
                }
                move_right = true;
            }
            else {
                const size_t index = i * N1;
                const size_t index_b = index - N1;
                XT[index] = XS[index - N1];
                Xb[index_b] = Xb[index_b + 1];
                j = N - (--i);
            }
        }
        if (move_right) {
            /* Find the next nonzero element. */
            size_t index = i * N1 + j;
            do {
                if (i + j < N) {
                    ++j;
                    ++index;
                }
                else {
                    ++i;
                    index += (N + 2) - j;
                    j = 1;
                }
            } while (i < N && X[index] == 0);
            if (i >= N)
                break;
            decrease = true;
        }
    }
}
}  // namespace steenrod

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
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_time (t SMALLINT, time REAL);");
    }

    void save_time(const std::string& table_prefix, int t, double time)
    {
        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_time (t, time) VALUES (?1, ?2);");
        stmt.bind_and_step(t, time);
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
    Mod2d load_d2(std::string_view table_prefix) const
    {
        Mod2d result;
        Statement stmt(*this, fmt::format("SELECT id, d2 FROM {} ORDER BY id;", table_prefix));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            int s = LocId(id).s;
            int v = LocId(id).v;
            Mod d2;
            d2.data = stmt.column_blob_tpl<MMod>(1);
            ut::get(ut::get(result, s), v) = std::move(d2);
        }
        return result;
    }
};

Mod ddd(const Mod& dg, const Mod2d& diffs, int s)
{
    Mod result;
    Milnor tmp1, tmp2;
    Mod tmp_m1, tmp_m2;
    std::vector<uint64_t> prod_sec;

    for (auto mdg : dg.data) {
        const auto& dv = diffs[size_t(s - 1)][mdg.v()];
        std::array<uint32_t, XI_MAX> Q = mdg.m().ToXi(), Q1;
        if (Contr(1, 0, Q, Q1)) {
            prod_sec.clear();
            for (auto mdv : dv.data) {
                const auto& dv_mdv = diffs[size_t(s - 2)][mdv.v()];
                std::array<uint32_t, XI_MAX> R = mdv.m().ToXi();
                for (auto mdv_mdv : dv_mdv.data) {
                    std::array<uint32_t, XI_MAX> S = mdv_mdv.m().ToXi();
                    MulMilnorMod4(R, S, prod_sec, mdv_mdv.v_raw() << 2);
                }
            }
            SortMod4(prod_sec);
            for (uint64_t x : prod_sec) {
                if ((x & 3) != 2) {
                    fmt::print("x should be a two torsion\n");
                    std::exit(-2);
                }
                x = (x >> 2) | (0x3ULL << 62);
                tmp1.data.clear();
                MulMilnor(Q1, MMilnor(x).ToXi(), tmp1.data);
                for (auto m : tmp1.data)
                    result.data.push_back(MMod(m.data() | MMod(x).v_raw()));
            }
        }
    }

    for (uint32_t m = 1; m <= XI_MAX_MULT; ++m) {
        for (uint32_t n = 1; n <= m; ++n) {
            tmp_m2.data.clear();
            for (auto mdg : dg.data) {
                const auto& dv = diffs[size_t(s - 1)][mdg.v()];
                std::array<uint32_t, XI_MAX> Q = mdg.m().ToXi(), Q1;
                for (uint32_t n1 = 1; n1 <= n; ++n1) {
                    for (uint32_t m1 = 0; m1 < n1; ++m1) {
                        if (Contr(m - m1, m1, n - n1, n1, Q, Q1)) {
                            tmp_m1.data.clear();
                            for (auto mdv : dv.data) {
                                const auto& dv_mdv = diffs[size_t(s - 2)][mdv.v()];
                                std::array<uint32_t, XI_MAX> R = mdv.m().ToXi(), R1;
                                for (auto& mdv_mdv : dv_mdv.data) {
                                    std::array<uint32_t, XI_MAX> S = mdv_mdv.m().ToXi(), S1;
                                    for (uint32_t k = 0; k <= m1; ++k) {
                                        if (Contr(m1 - k, k, n1 - k, k, R, R1) && Contr(k + 1, 0, S, S1)) {
                                            tmp1.data.clear();
                                            MulMilnor(R1, S1, tmp1.data);
                                            for (auto m : tmp1.data)
                                                tmp_m1.data.push_back(MMod(m.data() | mdv_mdv.v_raw()));
                                        }
                                    }
                                }
                            }
                            SortMod2(tmp_m1.data);
                            for (auto& m_R1S1 : tmp_m1.data) {
                                std::array<uint32_t, XI_MAX> R1S1 = m_R1S1.m().ToXi();
                                tmp1.data.clear();
                                MulMilnor(Q1, R1S1, tmp1.data);
                                for (auto m : tmp1.data)
                                    tmp_m2.data.push_back(MMod(m.data() | m_R1S1.v_raw()));
                            }
                        }
                    }
                }
            }
            SortMod2(tmp_m2.data);
            std::array<uint32_t, XI_MAX> M = {};
            ++M[size_t(m - 1)];
            ++M[size_t(n - 1)];
            for (auto& m_Q1R1S1 : tmp_m2.data) {
                std::array<uint32_t, XI_MAX> Q1R1S1 = m_Q1R1S1.m().ToXi();
                tmp1.data.clear();
                MulMilnor(M, Q1R1S1, tmp1.data);
                for (auto m : tmp1.data)
                    result.data.push_back(MMod(m.data() | m_Q1R1S1.v_raw()));
            }
        }
    }
    SortMod2(result.data);
    return result;
}

int GetCoh(int1d& v_degs, Mod1d& rels, int t_max, const std::string& name);

/*
 * Ext^2(cw) --> H^*(cw)
 *
 *  F_s -----f-----> F_{s-2}
 *   |   \            |
 *   d       ddd      d
 *   |            \   |
 *   V                V
 *  F_{s-1} --f--> F_{s-3}
 */
int compute_d2(const std::string& cw, int t_trunc, int stem_trunc)
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
    Mod2d all_f = dbD2.load_d2(table_d2);

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

    auto it = id_deg.begin();
    while (it != id_deg.end() && it->second.t < 2)
        ++it;
    std::vector<AdamsDegV2> degs;
    Mod2d f, fd;
    int3d fh;
    int1d arr_s, arr_v, arr_i;
    int t = 2;
    std::mutex print_mutex = {};
    /*for (size_t s = 1; s <= 2; ++s) {
        for (size_t i = 0; i < 10; ++i)
            fmt::print("s={} i={} d={}\n", s, i, diffs[s][i]);
    }*/

    for (; t <= t_trunc; ++t) {
        int1d ids;
        degs.clear();
        f.clear();
        fd.clear();
        fh.clear();
        arr_s.clear();
        arr_v.clear();
        arr_i.clear();
        for (; it != id_deg.end() && it->second.t < t + 1; ++it) {
            const auto& [id, deg] = *it;
            ids.push_back(id);
            degs.push_back(deg);
            const size_t diffs_d_size = num_diffs.at(deg);
            ut::get(f, deg.s) = Mod1d(diffs_d_size);
            ut::get(fh, deg.s) = int2d(diffs_d_size);
            ut::get(fd, deg.s) = Mod1d(diffs_d_size);

            if (deg.s >= 2) {
                auto v = LocId(id).v;
                for (size_t i = 0; i < diffs_d_size; ++i) {
                    arr_s.push_back(deg.s);
                    arr_v.push_back((int)(v + i));
                    arr_i.push_back((int)i);
                }
            }
        }

        /*# compute fd+A */
        std::atomic<int> threadsLeft = (int)arr_s.size();
        ut::for_each_par32(arr_s.size(), [&arr_s, &arr_v, &arr_i, &fd, &f, &all_f, &diffs, &print_mutex, &threadsLeft, t](size_t i) {  ////
            int s = arr_s[i];
            int v = arr_v[i];
            int j = arr_i[i];
            if (s > 2)
                fd[s][j] = subs(diffs[s][v], all_f[size_t(s - 1)]) + ddd(diffs[s][v], diffs, s);
            else if (s == 2) {
                // if (v == 0)
                //     f[s][j] = MMilnor::P(1, 2) * MMod(MMilnor::P(0, 1), 0);
            }

            {
                std::scoped_lock lock(print_mutex);
                --threadsLeft;
                fmt::print("t={} s={} threadsLeft={}\n", t, s, threadsLeft.load());
                std::fflush(stdout);
            }
        });

        /*# compute f */
        ut::for_each_par32(degs.size(), [&degs, &fd, &f, &gb](size_t i) {
            int s = degs[i].s;
            if (s >= 3)
                gb.DiffInvBatch(fd[s], f[s], size_t(s - 3));
        });

        /*# compute fh */
        for (auto& deg : degs) {
            int s = deg.s;
            for (size_t i = 0; i < f[s].size(); ++i)
                fh[s][i] = HomToK(f[s][i]);
        }

        dbD2.begin_transaction();
        stmt_time.step_and_reset();
        /*# save products to database */
        for (size_t i_id = 0; i_id < ids.size(); ++i_id) {
            int id = ids[i_id];
            int s = degs[i_id].s;
            for (size_t i = 0; i < f[s].size(); ++i) {
                stmt_map.bind_and_step(id + (int)i, f[s][i].data, myio::Serialize(fh[s][i]));
            }
        }
        double time = timer.Elapsed();
        timer.Reset();
        fmt::print("t={} time={}\n", t, time);
        std::fflush(stdout);
        dbD2.save_time(table_d2, t, time);
        stmt_t_max.bind_and_step(t);
        dbD2.end_transaction();

        /*# update all_f */
        for (size_t s = 0; s < f.size(); ++s)
            for (size_t i = 0; i < f[s].size(); ++i)
                ut::get(all_f, s).push_back(std::move(f[s][i]));
    }

    stmt_t_max.bind_and_step(std::max({t - 1, old_t_max_d2, t_trunc}));
    stmt_time.step_and_reset();
    return 0;
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

/* Prevent double run on linux */
#ifdef __linux__
    if (IsAdamsRunning(fmt::format("./Adams d2 {} ", cw))) {
        fmt::print("Error: ./Adams d2 {} is already running.\n", cw);
        return -1;
    }
#endif

    compute_d2(cw, t_max, DEG_MAX);
    return 0;
}