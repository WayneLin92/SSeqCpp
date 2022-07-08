#include "groebner_steenrod_const.h"
#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/linalg.h"
#include "algebras/utility.h"
#include <cstring>
#include <map>

namespace steenrod {

/** Return i such that mon divides leads[i].
 *
 * Return -1 if not found.
 * @param indices indices[v_raw] stores the indices of elements of leads with the given v_raw.
 */
int IndexOfDivisibleLeading(const MMod1d& leads, const std::unordered_map<uint64_t, int1d>& indices, MMod mon)
{
    auto key = mon.v_raw();
    auto p = indices.find(key);
    if (p != indices.end())
        for (int k : p->second)
            if (divisibleLF(leads[k], mon))
                return k;
    return -1;
}

/********************************************************
 *                    class SteenrodMResConst
 ********************************************************/

SteenrodMResConst::SteenrodMResConst(int t_trunc, int s_trunc, DataMResConst2d data, int2d basis_degrees) : t_trunc_(t_trunc), s_trunc_(s_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees))
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

Mod SteenrodMResConst::DiffInv(Mod x, size_t s) const
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
            x.iaddmul(m, gb_[s][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
            result.iaddmul(m, gb_[s][gb_index].x2, tmp_a, tmp_x1, tmp_x2);
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
            result.iaddmul(m, gb_[sp1][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
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

void SteenrodMResConst::DiffInvBatch(Mod1d xs, Mod1d& result, size_t s) const
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

class DbMResProd : public myio::Database
{
    using Statement = myio::Statement;

public:
    DbMResProd() = default;
    explicit DbMResProd(const std::string& filename) : Database(filename) {}

public:
    void create_generators(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, repr BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (id INTEGER PRIMARY KEY, g BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_dual(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, repr BLOB, s SMALLINT, t SMALLINT);");
    }

    void create_products(const std::string& table_prefix, const GenMRes2d& gens)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, s SMALLINT, t SMALLINT, indecomposable TINYINT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators_products (id_ind INTEGER, id INTEGER, prod BLOB, prod_h BLOB, PRIMARY KEY (id_ind, id));");

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

    void create_generators_and_delete(const std::string& table_prefix) const
    {
        create_generators(table_prefix);
        delete_from(table_prefix + "_generators");
    }
    void create_relations_and_delete(const std::string& table_prefix) const
    {
        create_relations(table_prefix);
        delete_from(table_prefix + "_relations");
    }

    void save_generators(const std::string& table_prefix, const DataMResConst1d& kernel, int s, int t) const  ////
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (diff, s, t) VALUES (?1, ?2, ?3);");

        for (auto& k : kernel) {
            stmt.bind_blob(1, k.x1.data.data(), int(k.x1.data.size() * sizeof(MMod)));
            stmt.bind_int(2, s);
            stmt.bind_int(3, t);
            stmt.step_and_reset();
        }
    }

    void load_products(const std::string& table_prefix, std::map<int, Mod2d>& map, std::map<int, int3d>& map_h) const
    {
        int1d id2s = get_column_int(table_prefix + "_generators", "s", "ORDER BY id;");
        Statement stmt(*this, "SELECT id, id_ind, prod, prod_h FROM " + table_prefix + "_generators_products ORDER BY id;");
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

    int2d load_basis_degrees(const std::string& table_prefix) const
    {
        int2d result;
        Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators ORDER BY id;");
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
        Statement stmt(*this, "SELECT id, s, t, diff FROM " + table_prefix + "_generators" + (t_trunc == DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_trunc)) + " ORDER BY id;");
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

    DataMResConst2d load_data(const std::string& table_prefix) const
    {
        DataMResConst2d data;
        Statement stmt(*this, "SELECT x1, x2, s FROM " + table_prefix + "_relations ORDER BY id;");
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

    void latest_st(const std::string& table_prefix, int& s, int& t) const
    {
        t = get_int("SELECT coalesce(max(t), 0) FROM " + table_prefix + "_relations;");
        s = get_int("SELECT coalesce(min(s), 1000) FROM " + table_prefix + "_relations WHERE t=" + std::to_string(t) + ";");
    }
};

SteenrodMResConst SteenrodMResConst::load(const std::string& filename)
{
    DbMResProd db(filename);
    DataMResConst2d data = db.load_data("SteenrodMRes");
    int2d basis_degrees = db.load_basis_degrees("SteenrodMRes");
    int latest_s = 0, latest_t = 0;
    db.latest_st("SteenrodMRes", latest_s, latest_t);
    return SteenrodMResConst(latest_t, latest_s, std::move(data), std::move(basis_degrees));
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
void compute_products_batch(const GenMRes2d& gens, const SteenrodMResConst& gb, const std::map<int, std::vector<std::pair<int, int>>>& leftFactors, std::map<int, Mod2d>& map, std::map<int, int3d>& map_h, DbMResProd& dbProd,
                            const std::string& table_prefix)
{
    static Mod1d tmp_x1(ut::FUTURE_NUM_THREADS);
    static Mod1d tmp_x2(ut::FUTURE_NUM_THREADS);
    static Milnor1d tmp_a(ut::FUTURE_NUM_THREADS);

    bench::Timer timer;
    timer.SuppressPrint();

    myio::Statement stmt_prod(dbProd, "INSERT OR IGNORE INTO " + table_prefix + "_generators_products (id_ind, id, prod, prod_h) VALUES (?1, ?2, ?3, ?4);");
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

/**
 * Compute products with the loaded indecomposables
 */
void compute_products_ind(int t_trunc, std::string& db_in, std::string& db_out)
{
    std::cout << "Compute with loaded indecomposables" << std::endl;

    std::string table_SteenrodMRes = "SteenrodMRes";
    auto gb = SteenrodMResConst::load(db_in);
    DbMResProd db(db_in);
    auto gens = db.load_generators(table_SteenrodMRes, t_trunc);

    int num_gens = 0;
    for (size_t i = 0; i < gens.size(); ++i)
        num_gens += (int)gens[i].size();

    std::map<int, Mod2d> map; /* `map` and `map_h` should be destructed after dbProd._future */
    std::map<int, int3d> map_h;

    DbMResProd dbProd(db_out);
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

    std::string table_SteenrodMRes = "SteenrodMRes";
    auto gb = SteenrodMResConst::load(db_in);
    DbMResProd db(db_in);
    auto gens = db.load_generators(table_SteenrodMRes, t_trunc);

    int num_gens = 0;
    for (size_t i = 0; i < gens.size(); ++i)
        num_gens += (int)gens[i].size();

    std::map<int, Mod2d> map; /* `map` and `map_h` should be destructed after  dbProd._future */
    std::map<int, int3d> map_h;

    DbMResProd dbProd(db_out);
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
