#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/linalg.h"
#include "algebras/utility.h"
#include "groebner_steenrod_const.h"
#include "main.h"
#include <cstring>

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

void DbAdamsResLoader::load_id_converter(const std::string& table_in, int2d& loc2glo, std::vector<std::pair<int, int>>& glo2loc) const
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

int2d DbAdamsResLoader::load_basis_degrees(const std::string& table_prefix, int t_trunc) const
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

void assign_vid_num(int2d& vid_num, AdamsDegV2 deg, int id_max)
{
    if ((int)vid_num.size() <= deg.s)
        vid_num.resize(size_t(deg.s + 1));
    size_t old_size = vid_num[deg.s].size();
    int filler = old_size ? vid_num[deg.s].back() : 0;
    vid_num[deg.s].resize(size_t(deg.stem() + 1));
    for (size_t i = old_size; i < deg.stem(); ++i)
        vid_num[deg.s][i] = filler;
    vid_num[deg.s][deg.stem()] = id_max;
}

void DbAdamsResLoader::load_generators(const std::string& table_prefix, std::map<int, AdamsDegV2>& id_st, int2d& vid_num, std::map<AdamsDegV2, Mod1d>& diffs, int t_trunc) const
{
    Statement stmt(*this, "SELECT id, s, t, diff FROM " + table_prefix + "_generators WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
    AdamsDegV2 d_prev = AdamsDegV2(-1, -1);
    std::map<int, int> vid_num_v2;
    while (stmt.step() == MYSQLITE_ROW) {
        int id = stmt.column_int(0);
        AdamsDegV2 d = AdamsDegV2(stmt.column_int(1), stmt.column_int(2));
        Mod diff;
        diff.data = stmt.column_blob_tpl<MMod>(3);
        diffs[d].push_back(std::move(diff));

        ++vid_num_v2[d.s];
        if (d.s != d_prev.s || d.t != d_prev.t) {
            id_st[id] = d;
            if (d_prev.s != -1)
                assign_vid_num(vid_num, d_prev, vid_num_v2[d_prev.s]);
            d_prev = d;
        }
    }
    if (d_prev.s != -1)
        assign_vid_num(vid_num, d_prev, vid_num_v2[d_prev.s]);
    for (int s = 0; s <= t_trunc; ++s) {
        int num = (vid_num.size() > s && vid_num[s].size()) ? vid_num[s].back() : 0;
        assign_vid_num(vid_num, AdamsDegV2(s, t_trunc), num);
    }
}

DataMResConst2d DbAdamsResLoader::load_data(const std::string& table_prefix, int t_trunc) const
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

class DbAdamsResProd : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsResProd(const std::string& filename) : Database(filename) {}

public:
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
};

AdamsResConst AdamsResConst::load(const DbAdamsResLoader& db, const std::string& table, int t_trunc)
{
    DataMResConst2d data = db.load_data(table, t_trunc);
    int2d basis_degrees = db.load_basis_degrees(table, t_trunc);
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

void compute_products_by_t(int t_trunc, const std::string& db_in, const std::string& table_in, const std::string& db_out)
{
    const std::string table_out = "S0_Adams_res";

    DbAdamsResLoader dbRes(db_in);

    auto gb = AdamsResConst::load(dbRes, table_in, t_trunc);
    std::map<int, AdamsDegV2> id_deg;
    int2d vid_num;
    std::map<AdamsDegV2, Mod1d> diffs;
    dbRes.load_generators(table_in, id_deg, vid_num, diffs, t_trunc);
    int2d loc2glo;
    std::vector<std::pair<int, int>> glo2loc;
    dbRes.load_id_converter(table_in, loc2glo, glo2loc);

    DbAdamsResProd dbProd(db_out);
    dbProd.create_products(table_out);
    dbProd.create_time(table_out);
    int id_start = dbProd.get_int("SELECT COALESCE(MAX(id), -1) FROM " + table_out + "_products") + 1;

    myio::Statement stmt_gen(dbProd, "INSERT INTO " + table_out + "_generators (id, indecomposable, s, t) values (?1, ?2, ?3, ?4);");
    myio::Statement stmt_ind(dbProd, "UPDATE " + table_out + "_generators SET indecomposable=1 WHERE id=?1 and indecomposable=0;");
    myio::Statement stmt_prod(dbProd, "INSERT INTO " + table_out + "_products (id, id_ind, prod, prod_h) VALUES (?1, ?2, ?3, ?4);");

    const Mod one = MMod(MMilnor(), 0);
    const int1d one_h = {0};

    bench::Timer timer;
    timer.SuppressPrint();

    for (const auto& [id, deg] : id_deg) {
        if (id < id_start) {
            diffs.erase(deg);
            continue;
        }

        const auto& diffs_d = diffs.at(deg);
        const size_t diffs_d_size = diffs_d.size();

        dbProd.begin_transaction();

        /* save generators to database */
        for (size_t i = 0; i < diffs_d_size; ++i) {
            stmt_gen.bind_int(1, id + (int)i);
            stmt_gen.bind_int(2, 0);
            stmt_gen.bind_int(3, deg.s);
            stmt_gen.bind_int(4, deg.t);
            stmt_gen.step_and_reset();
        }

        if (deg.t > 0) {
            std::map<int, Mod1d> f_sm1 = dbProd.load_products(table_out, deg.s - 1, glo2loc);

            /* compute fd */
            std::map<int, Mod1d> fd;
            int1d id_inds;
            for (auto& [id_ind, _] : f_sm1)
                id_inds.push_back(id_ind);

            int vid_num_sm1 = vid_num[size_t(deg.s - 1)][deg.stem()];
            for (auto& [id_ind, f_sm1_id_ind] : f_sm1) {
                f_sm1_id_ind.resize(size_t(vid_num_sm1));
                fd[id_ind].resize(diffs_d_size);
            }

            ut::for_each_par(diffs_d_size * id_inds.size(), [&id_inds, &fd, &diffs_d, &f_sm1, diffs_d_size](size_t i) {
                int id_ind = id_inds[i / diffs_d_size];
                size_t j = i % diffs_d_size;
                fd.at(id_ind)[j] = subs(diffs_d[j], f_sm1.at(id_ind));
            });

            /* compute f */
            std::map<int, Mod1d> f;
            int1d s1;
            for (auto& [id_ind, _] : fd) {
                s1.push_back(deg.s - 1 - glo2loc.at(id_ind).first);
                f[id_ind].resize(diffs_d_size);
            }

            ut::for_each_par(id_inds.size(), [&id_inds, &gb, &fd, &f, &s1](size_t i) { gb.DiffInvBatch(fd[id_inds[i]], f[id_inds[i]], s1[i]); });

            /* compute fh */
            std::map<int, int2d> fh;
            for (auto& [id_ind, f_id_ind] : f)
                for (size_t i = 0; i < diffs_d_size; ++i)
                    fh[id_ind].push_back(HomToK(f_id_ind[i]));

            /* save products to database */
            for (auto& [id_ind, f_id_ind] : f) {
                for (size_t i = 0; i < diffs_d_size; ++i) {
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
            int1d indices = lina::AddVectors(ut::int_range((int)diffs_d_size), lead_image);

            /* mark indicomposables in database */
            for (int i : indices) {
                stmt_ind.bind_int(1, id + i);
                stmt_ind.step_and_reset();
            }

            /* indicomposable comultiply with itself */
            for (int i : indices) {
                stmt_prod.bind_int(1, id + (int)i);
                stmt_prod.bind_int(2, id + (int)i);
                stmt_prod.bind_blob(3, one.data);
                stmt_prod.bind_blob(4, one_h);
                stmt_prod.step_and_reset();
            }
        }

        double time = timer.Elapsed();
        timer.Reset();
        std::cout << "t=" << deg.t << " s=" << deg.s << " time=" << time << std::endl;
        dbProd.save_time(table_out, deg.s, deg.t, time);

        dbProd.end_transaction(1000);
        diffs.erase(deg);
    }
}

// std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
int main_prod(int argc, char** argv, int index)
{
    int t_max = 100;
#ifdef MYDEPLOY
    std::string db_in = "AdamsE2.db";
    std::string table_in = "SteenrodMRes";
#else
    std::string db_in = "S0_Adams_res.db";
    std::string table_in = "S0_Adams_res";
#endif
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

    compute_products_by_t(t_max, db_in, table_in, db_out);

    // bench::Counter::print();
    return 0;
}
