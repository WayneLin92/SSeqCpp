#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/linalg.h"
#include "algebras/utility.h"
#include "groebner_steenrod_const.h"
#include <cstring>
#include <map>

namespace steenrod {

struct AdamsDegV2
{
    int s, t;
    constexpr AdamsDegV2() : s(0), t(0) {}
    constexpr AdamsDegV2(int s_, int t_) : s(s_), t(t_) {}
    bool operator<(const AdamsDegV2& rhs) const
    {
        if (t < rhs.t)
            return true;
        if (t > rhs.t)
            return false;
        if (s < rhs.s)
            return true;
        return false;
    };
};

/** Return i such that mon divides leads[i].
 *
 * Return -1 if not found.
 * @param indices indices[v_raw] stores the indices of elements of leads with the given v_raw.
 */
int IndexOfDivisibleLeading(const MMod1d& leads, const std::unordered_map<uint64_t, int1d>& indices, MMod mon);

/********************************************************
 *                    class AdamsResConst
 ********************************************************/

AdamsResConst::AdamsResConst(DataMResConst2d data, int2d basis_degrees) : gb_(std::move(data)), basis_degrees_(std::move(basis_degrees))
{
    if (basis_degrees_.empty())
        basis_degrees_.push_back({0});
    if (basis_degrees_[0].empty())
        basis_degrees_[0].push_back(0);

    leads_.resize(gb_.size());
    indices_.resize(gb_.size());

    for (size_t s = 0; s < gb_.size(); ++s) {
        for (int j = 0; j < (int)gb_[s].size(); ++j) {
            leads_[s].push_back(gb_[s][j].x1.GetLead());
            indices_[s][gb_[s][j].x1.GetLead().v_raw()].push_back(j);
        }
    }
}

Mod AdamsResConst::DiffInv(Mod x, size_t s) const
{
    Milnor tmp_a;
    Mod result, tmp_x1, tmp_x2;
    tmp_a.data.reserve(50);
    tmp_x1.data.reserve(100);
    tmp_x2.data.reserve(100);

    size_t index;
    index = 0;
    while (index < x.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], x.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(x.data[index], gb_[s][gb_index].x1.data[0]);
            x.iaddmulP(m, gb_[s][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
            result.iaddmulP(m, gb_[s][gb_index].x2, tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }
    if (x)
        throw MyException(0x4b901937U, "Something is wrong: d_inv(x) not well defined.");
    size_t sp1 = s + 1;
    index = 0;
    while (index < result.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[sp1], indices_[sp1], result.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.data[index], gb_[sp1][gb_index].x1.data[0]);
            result.iaddmulP(m, gb_[sp1][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }

    return result;
}

struct IndexMMod
{
    MMod m;
    unsigned i, index;
    bool operator<(IndexMMod rhs)
    {
        return rhs.m < m;
    }
};
using IndexMMod1d = std::vector<IndexMMod>;

void AdamsResConst::DiffInvBatch(Mod1d xs, Mod1d& result, size_t s) const
{
    Mod tmp_x, prod_x1, prod_x2;
    Milnor tmp_a;
    tmp_x.data.reserve(64);
    prod_x1.data.reserve(64);
    prod_x2.data.reserve(64);
    tmp_a.data.reserve(64);

    IndexMMod1d heap;
    for (size_t i = 0; i < xs.size(); ++i)
        if (xs[i])
            heap.push_back(IndexMMod{xs[i].data[0], (unsigned)i, 0});
    std::make_heap(heap.begin(), heap.end());

    while (!heap.empty()) {
        MMod term = heap.front().m;
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], term);
        if (gb_index != -1) {
            MMilnor m = divLF(term, gb_[s][gb_index].x1.data[0]);
            mulP(m, gb_[s][gb_index].x1, prod_x1, tmp_a);
            mulP(m, gb_[s][gb_index].x2, prod_x2, tmp_a);

            // int b = 0;
            while (!heap.empty() && heap.front().m == term) {
                /*if (b++)
                    bench::Counter(0);
                else
                    bench::Counter(1);*/

                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());

                xs[i].iaddP(prod_x1, tmp_x);
                result[i].iaddP(prod_x2, tmp_x);

                if (index < xs[i].data.size()) {
                    heap.back() = IndexMMod{xs[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
        else {
            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());
                if (++index < xs[i].data.size()) {
                    heap.back() = IndexMMod{xs[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
    }

    for (size_t i = 0; i < xs.size(); ++i)
        if (xs[i])
            throw MyException(0x277dc39aU, "Something is wrong: d_inv(x) not well defined.");

    size_t sp1 = s + 1;
    heap.clear();
    for (size_t i = 0; i < result.size(); ++i)
        if (result[i])
            heap.push_back(IndexMMod{result[i].GetLead(), (unsigned)i, 0});
    std::make_heap(heap.begin(), heap.end());

    while (!heap.empty()) {
        MMod term = heap.front().m;
        int gb_index = IndexOfDivisibleLeading(leads_[sp1], indices_[sp1], term);
        if (gb_index != -1) {
            MMilnor m = divLF(term, gb_[sp1][gb_index].x1.data[0]);
            mulP(m, gb_[sp1][gb_index].x1, prod_x1, tmp_a);

            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());

                result[i].iaddP(prod_x1, tmp_x);

                if (index < result[i].data.size()) {
                    heap.back() = IndexMMod{result[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
        else {
            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());
                if (++index < result[i].data.size()) {
                    heap.back() = IndexMMod{result[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
    }
}

class DbAdamsResProd : public myio::Database
{
    using Statement = myio::Statement;

public:
    DbAdamsResProd() = default;
    explicit DbAdamsResProd(const std::string& filename) : Database(filename) {}

public:
    /**
     * loc2glo[s][i] = id
     * glo2loc[id] = (s, i)
     *
     * From SteenrodMRes
     */
    void load_id_converter(const std::string& table_in, int2d& loc2glo, std::vector<std::pair<int, int>>& glo2loc) const
    {
        Statement stmt(*this, "SELECT id, s FROM " + table_in + "_generators ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0), s = stmt.column_int(1);
            if (loc2glo.size() <= (size_t)s)
                loc2glo.resize(size_t(s + 1));
            loc2glo[s].push_back(id);
            glo2loc.push_back(std::make_pair(s, (int)(loc2glo[s].size() - 1)));
        }
    }

    void create_products(const std::string& table_prefix, const GenMRes2d& gens)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, indecomposable TINYINT, s SMALLINT, t SMALLINT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_products (id INTEGER, id_ind INTEGER, prod BLOB, prod_h BLOB, PRIMARY KEY (id, id_ind));");

        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_generators (id, s, t, indecomposable) VALUES (?1, ?2, ?3, ?4);");

        begin_transaction();
        for (size_t s = 0; s < gens.size(); ++s) {
            for (size_t i = 0; i < gens[s].size(); ++i) {
                stmt.bind_int(1, gens[s][i].id);
                stmt.bind_int(2, (int)s);
                stmt.bind_int(3, gens[s][i].t);
                stmt.bind_int(4, 0);
                stmt.step_and_reset();
            }
        }
        end_transaction();
    }

    void create_products(const std::string& table_prefix)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, indecomposable TINYINT, s SMALLINT, t SMALLINT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_products (id INTEGER, id_ind INTEGER, prod BLOB, prod_h BLOB, PRIMARY KEY (id, id_ind));");
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

    void load_products(const std::string& table_prefix, std::map<int, Mod2d>& map, std::map<int, int3d>& map_h) const
    {
        int1d id2s = get_column_int(table_prefix + "_generators", "s", "ORDER BY id;");
        Statement stmt(*this, "SELECT id, id_ind, prod, prod_h FROM " + table_prefix + "_products ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            int s = id2s[id];
            int id_ind = stmt.column_int(1);
            Mod prod;
            prod.data = stmt.column_blob_tpl<MMod>(2);
            int1d prod_h = stmt.column_blob_tpl<int>(3);
            if (map[id_ind].size() <= s) {
                map[id_ind].resize(size_t(s + 1));
                map_h[id_ind].resize(size_t(s + 1));
            }
            map.at(id_ind)[s].push_back(std::move(prod));
            map_h.at(id_ind)[s].push_back(std::move(prod_h));
        }
    }

    /* result[id_ind][index] = image(v_index) in degree s */
    std::map<int, Mod1d> load_products(const std::string& table_prefix, int s, const std::vector<std::pair<int, int>>& glo2loc) const
    {
        std::map<int, Mod1d> result;
        Statement stmt(*this, "SELECT id, id_ind, prod FROM " + table_prefix + "_products LEFT JOIN " + table_prefix + "_generators USING(id) WHERE s=" + std::to_string(s) + " ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            int id_ind = stmt.column_int(1);
            Mod prod;
            prod.data = stmt.column_blob_tpl<MMod>(2);

            int index = glo2loc.at(id).second;
            if (result[id_ind].size() <= index)
                result.at(id_ind).resize(size_t(index + 1));
            result[id_ind][index] = std::move(prod);
        }
        return result;
    }

    int2d load_basis_degrees(const std::string& table_prefix, int t_trunc) const
    {
        int2d result;
        Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int s = stmt.column_int(0), t = stmt.column_int(1);
            if (result.size() <= (size_t)s)
                result.resize(size_t(s + 1));
            result[s].push_back(t);
        }
        return result;
    }

    GenMRes2d load_generators(const std::string& table_prefix, int t_trunc) const
    {
        GenMRes2d result;

        int id_min = get_int("SELECT MIN(id) FROM " + table_prefix + "_generators;");
        Statement stmt(*this, "SELECT id, s, t, diff FROM " + table_prefix + "_generators WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0) - id_min, s = stmt.column_int(1), t = stmt.column_int(2);
            Mod diff;
            diff.data = stmt.column_blob_tpl<MMod>(3);

            if (result.size() <= (size_t)s)
                result.resize(size_t(s + 1));
            result[s].push_back(GenMRes{id, t, std::move(diff)});
        }
        return result;
    }

    void load_generators(const std::string& table_prefix, std::map<int, AdamsDegV2>& id_st, std::map<AdamsDegV2, int>& vid_max, std::map<AdamsDegV2, Mod1d>& diffs, int t_trunc) const
    {
        Statement stmt(*this, "SELECT id, s, t, diff FROM " + table_prefix + "_generators WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
        AdamsDegV2 d_prev = AdamsDegV2(-1, -1);
        std::map<int, int> vid_max_v2;
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            AdamsDegV2 d = AdamsDegV2(stmt.column_int(1), stmt.column_int(2));
            Mod diff;
            diff.data = stmt.column_blob_tpl<MMod>(3);
            diffs[d].push_back(std::move(diff));

            ++vid_max_v2[d.s];
            if (d.s != d_prev.s || d.t != d_prev.t) {
                id_st[id] = d;
                if (d_prev.s != -1)
                    vid_max[d_prev] = vid_max_v2[d_prev.s];
                d_prev = d;
            }
        }
        if (d_prev.s != -1)
            vid_max[d_prev] = vid_max_v2[d_prev.s];
    }

    DataMResConst2d load_data(const std::string& table_prefix, int t_trunc) const
    {
        DataMResConst2d data;
        Statement stmt(*this, "SELECT x1, x2, s FROM " + table_prefix + "_relations WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            DataMResConst x;
            x.x1.data = stmt.column_blob_tpl<MMod>(0);
            x.x2.data = stmt.column_blob_tpl<MMod>(1);

            size_t s = (size_t)stmt.column_int(2);
            if (data.size() <= s)
                data.resize(s + 1);
            data[s].push_back(std::move(x));
        }
        return data;
    }
};

AdamsResConst AdamsResConst::load(const DbAdamsResProd& db, int t_trunc)
{
    DataMResConst2d data = db.load_data("SteenrodMRes", t_trunc);
    int2d basis_degrees = db.load_basis_degrees("SteenrodMRes", t_trunc);
    return AdamsResConst(std::move(data), std::move(basis_degrees));
}

int1d HomToK(const Mod& x)
{
    int1d result;
    for (MMod m : x.data)
        if (m.deg_m() == 0)
            result.push_back((int)m.v());
    return result;
}

/* Compute the products of gens with gens[[leftFactors]] */
void compute_products_batch(const GenMRes2d& gens, const AdamsResConst& gb, const std::map<int, std::vector<std::pair<int, int>>>& leftFactors, std::map<int, Mod2d>& map, std::map<int, int3d>& map_h, DbAdamsResProd& dbProd,
                            const std::string& table_prefix)
{
    static Mod1d tmp_x1(ut::FUTURE_NUM_THREADS);
    static Mod1d tmp_x2(ut::FUTURE_NUM_THREADS);
    static Milnor1d tmp_a(ut::FUTURE_NUM_THREADS);

    bench::Timer timer;
    timer.SuppressPrint();

    myio::Statement stmt_prod(dbProd, "INSERT OR IGNORE INTO " + table_prefix + "_products (id_ind, id, prod, prod_h) VALUES (?1, ?2, ?3, ?4);");
    const size_t gens_size = gens.size();

    /* multiply with one */
    {
        dbProd.begin_transaction();
        myio::Statement stmt_ind(dbProd, "UPDATE OR IGNORE " + table_prefix + "_generators SET indecomposable=1 WHERE id=?1 and indecomposable=0;");
        for (auto& [s, p] : leftFactors) {
            for (auto& [id, k] : p) {
                map[id].resize(gens_size);
                map_h[id].resize(gens_size);

                map[id][s].resize(gens[s].size());
                map_h[id][s].resize(gens[s].size());

                map[id][s][k] = MMod(MMilnor(), 0);
                map_h[id][s][k] = HomToK(map[id][s][k]);

                /* Save to database */
                stmt_ind.bind_int(1, id);
                stmt_ind.step_and_reset();
                for (size_t i = 0; i < gens[s].size(); ++i) {
                    stmt_prod.bind_int(1, id);
                    stmt_prod.bind_int(2, gens[s][i].id);
                    stmt_prod.bind_blob(3, map[id][s][i].data);
                    stmt_prod.bind_blob(4, map_h[id][s][i]);
                    stmt_prod.step_and_reset();
                }
            }
        }
        dbProd.end_transaction();
    }

    std::vector<unsigned> arr_s1;
    std::vector<unsigned> arr_i1;
    std::vector<unsigned> arr_id;
    std::vector<unsigned> arr_i2;
    for (size_t s2 = 0; s2 < gens_size - 2; ++s2) {
        size_t size_elements = 0;
        for (auto& [s, p] : leftFactors) {
            size_t s1 = s2 + 1 + (size_t)s;
            if (s1 < gens_size)
                size_elements += gens[s1].size() * p.size();
        }

        /* Compute fdv */
        Mod1d image_diff(size_elements);
        size_t start = 0;
        arr_s1.clear();
        arr_i1.clear();
        arr_id.clear();
        arr_i2.clear();
        for (auto& [s_, p_] : leftFactors) {
            int s = s_;
            auto& p = p_;
            size_t s1 = s2 + 1 + (size_t)s;
            if (s1 < gens_size) {
                for (auto& [id, k] : p) {
                    size_t old_map_s1_size = map[id][s1].size();
                    /* for (size_t i = old_map_s1_size; i < gens[s1].size(); ++i)
                         subsP(gens[s1][i].diff, map[id][s1 - 1], image_diff[start + i], tmp_a, tmp_x1, tmp_x2);*/
                    for (size_t i = old_map_s1_size; i < gens[s1].size(); ++i) {
                        arr_s1.push_back(unsigned(s1));
                        arr_i1.push_back(unsigned(i));
                        arr_id.push_back(unsigned(id));
                        arr_i2.push_back(unsigned(start + i));
                    }
                    start += gens[s1].size();
                }
            }
        }
        size_t nThreads = std::min(ut::FUTURE_NUM_THREADS, arr_s1.size());
        ut::for_each_par(nThreads, [&gens, &map, &image_diff, &arr_s1, &arr_i1, &arr_id, &arr_i2](size_t j) {
            for (size_t i = j; i < arr_s1.size(); i += ut::FUTURE_NUM_THREADS)
                subsP(gens[arr_s1[i]][arr_i1[i]].diff, map[arr_id[i]][size_t(arr_s1[i] - 1)], image_diff[arr_i2[i]], tmp_a[j], tmp_x1[j], tmp_x2[j]);
        });

        /* Compute d^{-1}fdv */
        Mod2d image_diff_by_t(DIM_MAX_RES);
        std::vector<std::pair<unsigned, unsigned>> indices;
        for (size_t i = 0; i < size_elements; ++i) {
            unsigned t = 0;
            if (image_diff[i]) {
                MMod mv = image_diff[i].GetLead();
                t = unsigned(mv.deg_m() + gb.basis_degrees(s2)[mv.v()]);
            }
            indices.push_back(std::make_pair(t, (unsigned)image_diff_by_t[t].size()));
            image_diff_by_t[t].push_back(std::move(image_diff[i]));
        }
        Mod2d image(DIM_MAX_RES);
        for (size_t t = 0; t < DIM_MAX_RES; ++t)
            image[t].resize(image_diff_by_t[t].size());
        size_t t_max = (size_t)DIM_MAX_RES;
        while (t_max > 0 && image[t_max - 1].empty())
            --t_max;

        ut::for_each_par(t_max, [&gb, &image_diff_by_t, &image, s2, t_max](size_t t) { gb.DiffInvBatch(image_diff_by_t[t_max - t], image[t_max - t], s2); });

        /* Save to map and map_h */
        start = 0;
        dbProd.begin_transaction();
        for (auto& [s, p] : leftFactors) {
            size_t s1 = s2 + 1 + (size_t)s;
            if (s1 < gens_size) {
                for (auto& [id, k] : p) {
                    size_t old_map_s1_size = map[id][s1].size();
                    map[id][s1].resize(gens[s1].size());
                    for (size_t i = old_map_s1_size; i < gens[s1].size(); ++i) {
                        auto& [t, index] = indices[start + i];
                        map[id][s1][i] = std::move(image[t][index]);
                    }
                    start += gens[s1].size();
                    for (size_t i = old_map_s1_size; i < gens[s1].size(); ++i) {
                        map_h[id][s1].push_back(HomToK(map[id][s1][i]));

                        /* Save to database */
                        stmt_prod.bind_int(1, id);
                        stmt_prod.bind_int(2, gens[s1][i].id);
                        stmt_prod.bind_blob(3, map[id][s1][i].data);
                        stmt_prod.bind_blob(4, map_h[id][s1][i]);
                        stmt_prod.step_and_reset();
                        dbProd.reg_transaction();

                        map[id][s1 - 1].clear();
                    }
                }
            }
        }
        dbProd.end_transaction(1000);

        double time = timer.Elapsed();
        timer.Reset();
        if (time > 0.5)
            std::cout << "    s2=" << s2 << ' ' << time << 's' << std::endl;
    }
}

void compute_products_by_t(int t_trunc, const std::string& db_in, const std::string& table_in, const std::string& db_out)
{
    const std::string table_out = "S0_Adams_res";

    DbAdamsResProd dbRes(db_in);
    {  ////////////////// Convert old version to new
        int id_min = dbRes.get_int("SELECT MIN(id) FROM " + table_in + "_generators;");
        if (id_min != 0) {
            if (id_min != 1)
                throw MyException(0x4473a085U, "Table SteenrodMRes_generators is probably broken.");
            dbRes.execute_cmd("update SteenrodMRes_generators set id=id-1;");
        }
    }

    auto gb = AdamsResConst::load(dbRes, t_trunc);
    std::map<int, AdamsDegV2> id_s_t;
    std::map<AdamsDegV2, int> vid_max;
    std::map<AdamsDegV2, Mod1d> diffs;
    dbRes.load_generators(table_in, id_s_t, vid_max, diffs, t_trunc);
    int2d loc2glo;
    std::vector<std::pair<int, int>> glo2loc;
    dbRes.load_id_converter(table_in, loc2glo, glo2loc);

    DbAdamsResProd dbProd(db_out);
    dbProd.create_products(table_out);
    dbProd.create_time(table_out);
    int id_start = dbProd.get_int("SELECT COALESCE(MAX(id), -1) FROM " + table_out + "_products") + 1;

    myio::Statement stmt_gen(dbProd, "INSERT INTO " + table_out + "_generators (id, indecomposable, s, t) values (?1, ?2, ?3, ?4);");
    myio::Statement stmt_ind(dbProd, "UPDATE OR IGNORE " + table_out + "_generators SET indecomposable=1 WHERE id=?1 and indecomposable=0;");
    myio::Statement stmt_prod(dbProd, "INSERT INTO " + table_out + "_products (id, id_ind, prod, prod_h) VALUES (?1, ?2, ?3, ?4);");

    const Mod one = MMod(MMilnor(), 0);
    const int1d one_h = {0};

    bench::Timer timer;
    timer.SuppressPrint();

    for (auto& [id, deg] : id_s_t) {
        if (id < id_start)
            continue;
        
        const auto& diffs_d = diffs.at(deg);
        const size_t diffs_size = diffs_d.size();

        /* save generators to database */
        for (size_t i = 0; i < diffs_size; ++i) {
            stmt_gen.bind_int(1, id + (int)i);
            stmt_gen.bind_int(2, 0);
            stmt_gen.bind_int(3, deg.s);
            stmt_gen.bind_int(4, deg.t);
            stmt_gen.step_and_reset();
        }

        if (deg.t > 0) {
            dbProd.begin_transaction();

            std::map<int, Mod1d> f_sm1 = dbProd.load_products(table_out, deg.s - 1, glo2loc);

            /* compute fd */
            std::map<int, Mod1d> fd;
            int1d id_inds;
            for (auto& [id_ind, _] : f_sm1)
                id_inds.push_back(id_ind);

            for (auto& [id_ind, f_sm1_id_ind] : f_sm1) {
                int vid_max_sm1 = 0;
                for (int t1 = deg.t; t1-- > 0;) {
                    AdamsDegV2 d1 = AdamsDegV2(deg.s - 1, t1);
                    if (vid_max.find(d1) != vid_max.end()) {
                        vid_max_sm1 = vid_max.at(d1);
                        break;
                    }
                }

                f_sm1_id_ind.resize(size_t(vid_max_sm1));
                fd[id_ind].resize(diffs_size);
            }

            ut::for_each_par(diffs_size * id_inds.size(), [&id_inds, &fd, &diffs_d, &f_sm1, diffs_size](size_t i) {
                int id_ind = id_inds[i / diffs_size];
                size_t j = i % diffs_size;
                fd.at(id_ind)[j] = subs(diffs_d[j], f_sm1.at(id_ind));
            });

            /* compute f */
            std::map<int, Mod1d> f;
            std::vector<size_t> s1;
            for (auto& [id_ind, fd_id_ind] : fd) {
                s1.push_back(size_t(deg.s - 1 - glo2loc.at(id_ind).first));
                f[id_ind].resize(diffs_size);
            }

            ut::for_each_par(id_inds.size(), [&id_inds, &gb, &fd, &f, &s1](size_t i) { gb.DiffInvBatch(fd[id_inds[i]], f[id_inds[i]], s1[i]); });

            /* compute fh */
            std::map<int, int2d> fh;
            for (auto& [id_ind, f_id_ind] : f)
                for (size_t i = 0; i < diffs_size; ++i)
                    fh[id_ind].push_back(HomToK(f_id_ind[i]));

            /* save products to database */
            for (auto& [id_ind, f_id_ind] : f) {
                for (size_t i = 0; i < diffs_size; ++i) {
                    if (f_id_ind[i]) {
                        stmt_prod.bind_int(1, id + (int)i);
                        stmt_prod.bind_int(2, id_ind);
                        stmt_prod.bind_blob(3, f_id_ind[i].data);
                        stmt_prod.bind_blob(4, fh.at(id_ind)[i]);
                        stmt_prod.step_and_reset();
                        dbProd.reg_transaction();
                    }
                }
            }

            /* find indecomposables */
            int2d fx;
            for (const auto& [id_ind, fh_id_ind] : fh) {
                size_t offset = fx.size();
                for (size_t i = 0; i < fh_id_ind.size(); ++i) {
                    for (int k : fh_id_ind[i]) {
                        if (fx.size() <= offset + (size_t)k)
                            fx.resize(offset + (size_t)k + 1);
                        fx[offset + (size_t)k].push_back((int)i);
                    }
                }
            }
            int1d lead_image = lina::GetLeads(lina::GetSpace(fx));
            int1d indices = lina::AddVectors(ut::int_range((int)diffs_size), lead_image);

            /* mark indicomposables in database */
            for (int i : indices) {
                stmt_ind.bind_int(1, id + i);
                stmt_ind.step_and_reset();
            }

            /* save product x*1 to database */
            for (int i : indices) {
                stmt_prod.bind_int(1, id + (int)i);
                stmt_prod.bind_int(2, id + (int)i);
                stmt_prod.bind_blob(3, one.data);
                stmt_prod.bind_blob(4, one_h);
                stmt_prod.step_and_reset();
            }

            dbProd.end_transaction(1000);
        }

        double time = timer.Elapsed();
        timer.Reset();
        std::cout << "t=" << deg.t << " s=" << deg.s << " time=" << time << std::endl;
        dbProd.save_time(table_out, deg.s, deg.t, time);

        diffs.erase(deg);
    }
}

/**
 * Compute products with the loaded indecomposables
 */
void compute_products_ind(int t_trunc, std::string& db_in, std::string& db_out)
{
    std::cout << "Compute with loaded indecomposables" << std::endl;
    DbAdamsResProd dbRes(db_in);

    std::string table_SteenrodMRes = "SteenrodMRes";
    auto gb = AdamsResConst::load(dbRes, DEG_MAX);
    auto gens = dbRes.load_generators(table_SteenrodMRes, t_trunc);

    int num_gens = 0;
    for (size_t i = 0; i < gens.size(); ++i)
        num_gens += (int)gens[i].size();

    std::map<int, Mod2d> map; /* `map` and `map_h` should be destructed after dbProd._future */
    std::map<int, int3d> map_h;

    DbAdamsResProd dbProd(db_out);
    dbProd.create_products(table_SteenrodMRes, gens);
    dbProd.load_products(table_SteenrodMRes, map, map_h);

    int1d id_inds = dbProd.get_column_int(table_SteenrodMRes + "_generators", "id", "WHERE indecomposable=1 ORDER BY id");
    std::map<int, std::vector<std::pair<int, int>>> leftFactors;
    int1d ids;
    for (size_t s = 1; s < gens.size(); ++s) {
        for (size_t i = 0; i < gens[s].size(); ++i) {
            int id = gens[s][i].id;
            if (std::binary_search(id_inds.begin(), id_inds.end(), id)) {
                leftFactors[(int)s].push_back(std::make_pair((int)id, (int)i));
                ids.push_back(id);
            }
        }
    }

    compute_products_batch(gens, gb, leftFactors, map, map_h, dbProd, table_SteenrodMRes);
}

void compute_products(int t_trunc, std::string& db_in, std::string& db_out)
{
    std::cout << "Compute in t<=" << std::to_string(t_trunc) << std::endl;
    DbAdamsResProd dbRes(db_in);

    std::string table_SteenrodMRes = "SteenrodMRes";
    auto gb = AdamsResConst::load(dbRes, DEG_MAX);
    auto gens = dbRes.load_generators(table_SteenrodMRes, t_trunc);

    int num_gens = 0;
    for (size_t i = 0; i < gens.size(); ++i)
        num_gens += (int)gens[i].size();

    std::map<int, Mod2d> map; /* `map` and `map_h` should be destructed after  dbProd._future */
    std::map<int, int3d> map_h;

    DbAdamsResProd dbProd(db_out);
    dbProd.create_products(table_SteenrodMRes, gens);
    dbProd.load_products(table_SteenrodMRes, map, map_h);

    int1d id_inds = dbProd.get_column_int(table_SteenrodMRes + "_generators", "id", "WHERE indecomposable=1 and t<=" + std::to_string(t_trunc) + " ORDER BY id");
    for (size_t s = 1; s < gens.size(); ++s) {
        std::cout << "s1=" << s << std::endl;

        /* Compute the indecomposables in s */
        int2d fx;
        for (const auto& [i, map_h_i] : map_h) {
            if (map_h_i.size() <= (size_t)s)
                continue;
            size_t offset = fx.size();
            for (size_t k = 0; k < map_h_i[s].size(); ++k) {
                for (int l : map_h_i[s][k]) {
                    if (fx.size() <= offset + (size_t)l)
                        fx.resize(offset + (size_t)l + 1);
                    fx[offset + (size_t)l].push_back((int)k);
                }
            }
        }
        int2d image = lina::GetSpace(fx);
        int1d r = ut::int_range((int)gens[s].size());
        int1d lead_image;
        for (const int1d& a : image)
            lead_image.push_back(a[0]);
        std::sort(lead_image.begin(), lead_image.end());
        int1d indices = lina::AddVectors(r, lead_image);
        for (size_t i = 0; i < gens[s].size(); ++i)
            if (std::binary_search(id_inds.begin(), id_inds.end(), gens[s][i].id))
                indices.push_back((int)i);
        std::sort(indices.begin(), indices.end());

        std::map<int, std::vector<std::pair<int, int>>> leftFactors;
        int1d ids;
        for (size_t i = 0; i < indices.size(); ++i) {
            if (indices[i] >= (int)gens[s].size())
                break;
            int id = gens[s][indices[i]].id;
            leftFactors[(int)s].push_back(std::make_pair((int)id, (int)indices[i]));
            ids.push_back(id);
        }

        compute_products_batch(gens, gb, leftFactors, map, map_h, dbProd, table_SteenrodMRes);
    }
}

}  // namespace steenrod

// std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
int main_prod(int argc, char** argv, int index)
{
    int t_max = 392;
    std::string db_in = "AdamsE2.db";
    std::string table_in = "SteenrodMRes";
    std::string db_out = "S0_Adams_res_prod.db";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  Adams prod [t_max] [db_in] [table_in] [db_out]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  t_max = " << t_max << "\n";
        std::cout << "  db_in = " << db_in << "\n";
        std::cout << "  table_in = " << table_in << "\n";
        std::cout << "  db_out = " << db_out << "\n";

        std::cout << "Version:\n  2.0 (2022-08-07)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_in", db_in))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_in", table_in))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_out", db_out))
        return index;
    bench::Timer timer;

    // steenrod::compute_products_ind(t_max, db_in, db_out);
    // steenrod::compute_products(t_max, db_in, db_out);

    steenrod::compute_products_by_t(t_max, db_in, table_in, db_out);

    // bench::Counter::print();
    return 0;
}
