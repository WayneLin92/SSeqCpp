#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/linalg.h"
#include "algebras/utility.h"
#include "algebras/myhash.h"////
#include "groebner_steenrod_const.h"
#include "main.h"
#include <cstring>
#include <set>

class DbAdamsResProdLoader : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsResProdLoader(const std::string& filename) : Database(filename) {}

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

class DbAdamsResCell : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsResCell(const std::string& filename) : Database(filename) {}

public:
    void create_products(const std::string& table_prefix)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS S0_Adams_hi_products (id INTEGER PRIMARY KEY, prod_h TEXT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_E2 (id INTEGER PRIMARY KEY, cell SMALLINT, indecomposable TINYINT, h_repr TEXT, coh_repr TEXT, s SMALLINT, t SMALLINT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_products_whole (S0_id INTEGER, cell SMALLINT, id_ind INTEGER, prod BLOB, prod_h TEXT, s SMALLINT, PRIMARY KEY (S0_id, cell, id_ind));");
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_products (id INTEGER, id_ind INTEGER, prod_h TEXT, PRIMARY KEY (id, id_ind));");
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
    std::map<int, Mod1d> load_products(const std::string& table_prefix, int cell, int s, const std::vector<std::pair<int, int>>& glo2loc) const
    {
        std::map<int, Mod1d> result;
        Statement stmt(*this, "SELECT S0_id, id_ind, prod FROM " + table_prefix + "_products_whole WHERE cell=" + std::to_string(cell) + " AND s=" + std::to_string(s) + " ORDER BY S0_id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            if (id >= (int)glo2loc.size())
                break;
            int id_ind = stmt.column_int(1);
            Mod prod;
            prod.data = stmt.column_blob_tpl<MMod>(2);

            int index = glo2loc[id].second;
            if (result[id_ind].size() <= index)
                result.at(id_ind).resize(size_t(index + 1));
            result[id_ind][index] = std::move(prod);
        }
        return result;
    }

    /* result[id_ind][index] = image(v_index) in degree s */
    std::map<int, int> load_id_ind_to_s(const std::string& table_prefix) const
    {
        std::map<int, int> result;
        Statement stmt(*this, "SELECT id, s FROM " + table_prefix + "_E2 WHERE indecomposable=1 ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            int s = stmt.column_int(1);
            result[id] = s;
        }
        return result;
    }
};

int1d HomToKCell(const Mod& x, int t_cell)
{
    int1d result;
    for (MMod m : x.data)
        if (m.deg_m() == t_cell && m.w() == 1)
            result.push_back((int)m.v());
    return result;
}

Mod Comult(const Mod& x, int t_cell)
{
    Mod result;
    for (MMod m : x.data) {
        auto xi = m.m().ToXi();
        if ((int)xi[0] >= t_cell) {
            xi[0] -= t_cell;
            MMilnor m1 = MMilnor::Xi(xi.data());
            result += MMod(m1, m.v());
        }
    }
    return result;
}

/* x+=y */
void add_prod_h(int1d& x, const int1d& y, int1d& tmp)
{
    tmp.clear();
    std::set_symmetric_difference(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(tmp), std::greater<int>());
    std::swap(x, tmp);
};

namespace ut {
template <>
inline uint64_t hash<MMod>(const MMod& x)
{
    return uint64_t(x.data());
}
template <>
inline uint64_t hash<Mod>(const Mod& x)
{
    return ut::hash(x.data);
}
template <>
inline uint64_t hash<AdamsDegV2>(const AdamsDegV2& x)
{
    return ut::hash(uint64_t(x.s) + (uint64_t(x.t) << 32));
}
template <>
inline uint64_t hash<std::pair<int, int>>(const std::pair<int, int>& x)
{
    uint64_t seed = (uint64_t)x.first;
    ut::hash_combine(seed, x.second);
    return seed;
}
}

void compute_2cell_products_by_t(int t_trunc, const std::string& complex_name, const std::string& db_in, const std::string& table_in)
{
    int t_cell = 0;
    if (complex_name == "C2")
        t_cell = 1;
    else if (complex_name == "Ceta")
        t_cell = 2;
    else if (complex_name == "Cnu")
        t_cell = 4;
    else if (complex_name == "Csigma")
        t_cell = 8;
    else
        throw MyException(0xd2e1e176U, complex_name + " is not supported");

    //const std::string db_in_prod = "S0_Adams_res_prod.db";
    //const std::string table_in_prod = "S0_Adams_res";
    const std::string db_out = complex_name + "_Adams_chain.db";
    const std::string table_out = complex_name + "_Adams";

    DbAdamsResLoader dbRes(db_in);
    auto gb = AdamsResConst::load(dbRes, t_trunc);
    std::map<int, AdamsDegV2> id_deg;
    int2d vid_num;
    std::map<AdamsDegV2, Mod1d> diffs;
    dbRes.load_generators(table_in, id_deg, vid_num, diffs, t_trunc);
    int2d loc2glo;
    std::vector<std::pair<int, int>> glo2loc;
    dbRes.load_id_converter(table_in, loc2glo, glo2loc);

    //DbAdamsResProdLoader dbResProd(db_in_prod);

    DbAdamsResCell dbProd(db_out);
    dbProd.create_products(table_out);
    dbProd.create_time(table_out);

    int latest_s = dbProd.get_int("SELECT s FROM " + table_out + "_E2 WHERE cell=0 ORDER BY id DESC LIMIT 1", -100);
    int latest_t = dbProd.get_int("SELECT t FROM " + table_out + "_E2 WHERE cell=0 ORDER BY id DESC LIMIT 1", -100);
    int latest_s1 = dbProd.get_int("SELECT s FROM " + table_out + "_E2 WHERE cell=1 ORDER BY id DESC LIMIT 1", -100);
    int latest_t1 = dbProd.get_int("SELECT t FROM " + table_out + "_E2 WHERE cell=1 ORDER BY id DESC LIMIT 1", -100);
    AdamsDegV2 latest_deg = AdamsDegV2(latest_s, latest_t);
    AdamsDegV2 latest_deg1 = AdamsDegV2(latest_s1 + 1, latest_t1);
    if (latest_deg < latest_deg1)
        latest_deg = latest_deg1;
    std::map<int, int> id_ind_to_s = dbProd.load_id_ind_to_s(table_out);
    int gen_id = dbProd.get_int("SELECT COALESCE(MAX(id) + 1, 0) FROM " + table_out + "_E2");

    myio::Statement stmt_hi_prod(dbProd, "INSERT INTO S0_Adams_hi_products (id, prod_h) values (?1, ?2);");
    myio::Statement stmt_gen(dbProd, "INSERT INTO " + table_out + "_E2 (id, cell, indecomposable, h_repr, coh_repr, s, t) values (?1, ?2, ?3, ?4, ?5, ?6, ?7);");
    myio::Statement stmt_ind(dbProd, "UPDATE " + table_out + "_E2 SET indecomposable=1 WHERE id=?1 and indecomposable=0;");
    myio::Statement stmt_prod_whole(dbProd, "INSERT INTO " + table_out + "_products_whole (S0_id, cell, id_ind, prod, prod_h, s) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");
    myio::Statement stmt_prod(dbProd, "INSERT INTO " + table_out + "_products (id, id_ind, prod_h) VALUES (?1, ?2, ?3);");

    const Mod one = MMod(MMilnor(), 0);
    const int1d one_h = {0};
    Mod1d dummy_diffs;

    bench::Timer timer;
    timer.SuppressPrint();

    std::map<AdamsDegV2, int> deg_id;
    for (const auto& [id, deg] : id_deg) {
        deg_id[deg] = id;
        if (deg.t + t_cell <= t_trunc)
            deg_id[AdamsDegV2(deg.s + 1, deg.t + t_cell)] = -1;
    }

    for (const auto& [deg, id] : deg_id) {
        AdamsDegV2 deg1 = AdamsDegV2(deg.s - 1, deg.t - t_cell);

        if (!(latest_deg < deg)) {
            diffs.erase(deg1);
            continue;
        }

        dbProd.begin_transaction();

        const auto& diffs_d = diffs.find(deg) != diffs.end() ? diffs.at(deg) : dummy_diffs;
        const auto& diffs_d1 = diffs.find(deg1) != diffs.end() ? diffs.at(deg1) : dummy_diffs;
        const size_t diffs_d_size = diffs_d.size();
        const size_t diffs_d1_size = diffs_d1.size();

        Mod1d diffs_d_cell1;
        for (size_t i = 0; i < diffs_d_size; ++i)
            diffs_d_cell1.push_back(Comult(diffs_d[i], t_cell));

        if (deg == AdamsDegV2(3, 6))
            std::cout << '\n';

        /* Save the hi (co)products */
        int2d prod_hi;
        std::map<int, int1d> prod_hi_dual;
        if (deg.s - 1 >= 0) {
            int vid_start = deg.stem() >= t_cell ? vid_num[size_t(deg.s - 1)][size_t(deg.stem() - t_cell)] : 0;
            int vid_end = deg.stem() + 1 >= t_cell ? vid_num[size_t(deg.s - 1)][size_t(deg.stem() + 1 - t_cell)] : 0;
            for (int vid = vid_start; vid < vid_end; ++vid)
                prod_hi_dual[vid] = int1d{};
        }
        for (size_t i = 0; i < diffs_d_size; ++i) {
            prod_hi.push_back(HomToKCell(diffs_d[i], t_cell));
            for (int j : prod_hi.back())
                prod_hi_dual[j].push_back((int)i);

            stmt_hi_prod.bind_int(1, id + (int)i);
            stmt_hi_prod.bind_str(2, myio::Serialize(prod_hi.back()));
            stmt_hi_prod.step_and_reset();
            dbProd.reg_transaction();
        }

        /* Compute the kernel of the map (cell 0) and the kernel of the dual map (cell 1) */
        int2d kernel_ht, kernel_ht_dual, _x1, _x2;
        lina::SetLinearMap(prod_hi, _x1, kernel_ht, _x2);
        lina::SimplifySpace(kernel_ht);
        int1d x_prod_hi_dual;
        int2d fx_prod_hi_dual;
        for (auto& [xi, f] : prod_hi_dual) {
            x_prod_hi_dual.push_back(loc2glo[size_t(deg.s - 1)][xi]);
            fx_prod_hi_dual.push_back(std::move(f));
        }
        _x1.clear();
        _x2.clear();
        lina::SetLinearMapV2(x_prod_hi_dual, fx_prod_hi_dual, _x1, kernel_ht_dual, _x2);
        lina::SimplifySpace(kernel_ht_dual);

        /* save generators to database */
        int gen_id_cell0_start = gen_id;
        for (size_t i = 0; i < kernel_ht.size(); ++i) {
            int1d kernel_offset;
            for (size_t j = 0; j < kernel_ht[i].size(); ++j)
                kernel_offset.push_back(id + kernel_ht[i][j]);

            stmt_gen.bind_int(1, gen_id++);
            stmt_gen.bind_int(2, 0);
            stmt_gen.bind_int(3, 0);
            stmt_gen.bind_str(4, myio::Serialize(kernel_offset));
            stmt_gen.bind_str(5, myio::Serialize(int1d{kernel_offset.front()}));
            stmt_gen.bind_int(6, deg.s);
            stmt_gen.bind_int(7, deg.t);
            stmt_gen.step_and_reset();
            dbProd.reg_transaction();
        }
        int gen_id_cell1_start = gen_id;
        for (size_t i = 0; i < kernel_ht_dual.size(); ++i) {
            stmt_gen.bind_int(1, gen_id++);
            stmt_gen.bind_int(2, 1);
            stmt_gen.bind_int(3, 0);
            stmt_gen.bind_str(4, myio::Serialize(int1d{kernel_ht_dual[i].front()}));
            stmt_gen.bind_str(5, myio::Serialize(kernel_ht_dual[i]));
            stmt_gen.bind_int(6, deg.s - 1);
            stmt_gen.bind_int(7, deg.t);
            stmt_gen.step_and_reset();
            dbProd.reg_transaction();
        }

        if (deg.t > 0) {
            if (deg1.t > 0 && diffs_d1_size) {
                std::map<int, Mod1d> f_cell1_sm2 = dbProd.load_products(table_out, 1, deg.s - 2, glo2loc);

                /* compute fd */
                std::map<int, Mod1d> fd;
                int1d id_inds;
                for (auto& [id_ind, _] : f_cell1_sm2)
                    id_inds.push_back(id_ind);

                int vid_num_sm2 = vid_num[size_t(deg1.s - 1)][deg1.stem()];
                for (int id_ind : id_inds) {
                    f_cell1_sm2[id_ind].resize(size_t(vid_num_sm2));
                    fd[id_ind].resize(diffs_d1_size);
                }

                ut::for_each_par(diffs_d1_size * id_inds.size(), [&id_inds, &fd, &diffs_d1, &f_cell1_sm2, diffs_d1_size](size_t i) {
                    int id_ind = id_inds[i / diffs_d1_size];
                    size_t j = i % diffs_d1_size;
                    fd.at(id_ind)[j] = subs(diffs_d1[j], f_cell1_sm2.at(id_ind));
                });

                //std::cout << "ut::hash(id_ind_to_s)=" << ut::hash(id_ind_to_s) << '\n';  //////////////////

                /* compute f */
                std::map<int, Mod1d> f;
                int1d s1;
                for (auto& [id_ind, _] : fd) {
                    s1.push_back(deg1.s - 1 - id_ind_to_s.at(id_ind));
                    f[id_ind].resize(diffs_d1_size);
                }

                ut::for_each_par(id_inds.size(), [&id_inds, &gb, &fd, &f, &s1](size_t i) { gb.DiffInvBatch(fd[id_inds[i]], f[id_inds[i]], s1[i]); });

                /* compute fh */
                std::map<int, int2d> fh;
                for (auto& [id_ind, f_id_ind] : f)
                    for (size_t i = 0; i < diffs_d1_size; ++i)
                        fh[id_ind].push_back(HomToK(f_id_ind[i]));

                /* save products to database */
                int vid_d1 = deg1.stem() > 0 ? vid_num[deg1.s][size_t(deg1.stem() - 1)] : 0;
                int id_d1 = loc2glo[deg1.s][vid_d1];
                for (auto& [id_ind, f_id_ind] : f) {
                    for (size_t i = 0; i < diffs_d1_size; ++i) {
                        if (f_id_ind[i]) {
                            stmt_prod_whole.bind_int(1, id_d1 + (int)i);
                            stmt_prod_whole.bind_int(2, 1);
                            stmt_prod_whole.bind_int(3, id_ind);
                            stmt_prod_whole.bind_blob(4, f_id_ind[i].data);
                            stmt_prod_whole.bind_str(5, myio::Serialize(fh.at(id_ind)[i]));
                            stmt_prod_whole.bind_int(6, deg1.s); 
                            stmt_prod_whole.step_and_reset();
                            dbProd.reg_transaction();
                        }
                    }
                }

                int2d fx;
                for (auto& [id_ind, fh_id_ind] : fh) {
                    size_t offset = fx.size();
                    for (size_t i = 0; i < kernel_ht_dual.size(); ++i) {
                        int index = kernel_ht_dual[i].front() - id_d1;
                        if (!fh_id_ind[index].empty()) {
                            stmt_prod.bind_int(1, gen_id_cell1_start + (int)i);
                            stmt_prod.bind_int(2, id_ind);

                            stmt_prod.bind_str(3, myio::Serialize(fh_id_ind[index]));
                            stmt_prod.step_and_reset();
                            dbProd.reg_transaction();

                            for (int k : fh_id_ind[index]) {
                                if (fx.size() <= offset + (size_t)k)
                                    fx.resize(offset + (size_t)k + 1);
                                fx[offset + (size_t)k].push_back((int)i);
                            }
                        }
                    }
                }

                /* find indecomposables */
                int1d lead_image = lina::GetLeads(lina::GetSpace(fx));
                int1d indices = lina::AddVectors(ut::int_range((int)kernel_ht_dual.size()), lead_image);

                /* mark indicomposables in database */
                for (int i : indices) {
                    stmt_ind.bind_int(1, gen_id_cell1_start + i);
                    stmt_ind.step_and_reset();

                    id_ind_to_s[gen_id_cell1_start + i] = deg1.s;
                }

                /* cell1 id_ind comultiplies with itself */
                for (int i : indices) {
                    for (int j : kernel_ht_dual[i]) {
                        stmt_prod_whole.bind_int(1, j);
                        stmt_prod_whole.bind_int(2, 1);
                        stmt_prod_whole.bind_int(3, gen_id_cell1_start + (int)i);
                        stmt_prod_whole.bind_blob(4, one.data);
                        stmt_prod_whole.bind_str(5, myio::Serialize(one_h));
                        stmt_prod_whole.bind_int(6, deg1.s);
                        stmt_prod_whole.step_and_reset();
                    }

                    stmt_prod.bind_int(1, gen_id_cell1_start + (int)i);
                    stmt_prod.bind_int(2, gen_id_cell1_start + (int)i);
                    stmt_prod.bind_str(3, myio::Serialize(one_h));
                    stmt_prod.step_and_reset();
                    dbProd.reg_transaction();
                }
            }

            if (diffs_d_size) {
                std::map<int, Mod1d> f_sm1 = dbProd.load_products(table_out, 0, deg.s - 1, glo2loc);
                std::map<int, Mod1d> f_cell1_sm1 = dbProd.load_products(table_out, 1, deg.s - 1, glo2loc);

                std::set<int> set_id_inds;
                for (auto& [id_ind, _] : f_sm1)
                    set_id_inds.insert(id_ind);
                for (auto& [id_ind, _] : f_cell1_sm1)
                    set_id_inds.insert(id_ind);
                int1d id_inds;
                for (int id_ind : set_id_inds)
                    id_inds.push_back(id_ind);

                /* compute fd */
                std::map<int, Mod1d> fd;
                int vid_num_sm1 = vid_num[size_t(deg.s - 1)][deg.stem()];
                for (int id_ind : id_inds) {
                    f_sm1[id_ind].resize(size_t(vid_num_sm1));
                    f_cell1_sm1[id_ind].resize(size_t(vid_num_sm1));
                    fd[id_ind].resize(diffs_d_size);
                }

                ut::for_each_par(diffs_d_size * id_inds.size(), [&id_inds, &fd, &diffs_d, &diffs_d_cell1, &f_sm1, &f_cell1_sm1, diffs_d_size](size_t i) {
                    int id_ind = id_inds[i / diffs_d_size];
                    size_t j = i % diffs_d_size;
                    fd.at(id_ind)[j] = subs(diffs_d[j], f_sm1.at(id_ind)) + subs(diffs_d_cell1[j], f_cell1_sm1.at(id_ind));
                });

                /* compute f */
                std::map<int, Mod1d> f;
                std::vector<size_t> s1;
                for (auto& [id_ind, fd_id_ind] : fd) {
                    s1.push_back(size_t(deg.s - 1 - id_ind_to_s.at(id_ind)));
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
                            stmt_prod_whole.bind_int(1, id + (int)i);
                            stmt_prod_whole.bind_int(2, 0);
                            stmt_prod_whole.bind_int(3, id_ind);
                            stmt_prod_whole.bind_blob(4, f_id_ind[i].data);
                            stmt_prod_whole.bind_str(5, myio::Serialize(fh.at(id_ind)[i]));
                            stmt_prod_whole.bind_int(6, deg.s);
                            stmt_prod_whole.step_and_reset();
                            dbProd.reg_transaction();
                        }
                    }
                }
                int1d tmp_prod_h;
                for (auto& [id_ind, fh_id_ind] : fh) { /* multiply with id_ind */
                    for (size_t i = 0; i < kernel_ht.size(); ++i) {
                        int1d prod_h;
                        for (int j : kernel_ht[i])
                            add_prod_h(prod_h, fh.at(id_ind)[j], tmp_prod_h);
                        if (!prod_h.empty()) {
                            stmt_prod.bind_int(1, gen_id_cell0_start + (int)i);
                            stmt_prod.bind_int(2, id_ind);
                            stmt_prod.bind_str(3, myio::Serialize(prod_h));
                            stmt_prod.step_and_reset();
                            dbProd.reg_transaction();
                        }
                    }
                }
                for (size_t i = 0; i < kernel_ht.size(); ++i) { /* multiply with the image of one */
                    int1d kernel_offset;
                    for (size_t j = 0; j < kernel_ht[i].size(); ++j)
                        kernel_offset.push_back(id + kernel_ht[i][j]);
                    stmt_prod.bind_int(1, gen_id_cell0_start + (int)i);
                    stmt_prod.bind_int(2, 0);
                    stmt_prod.bind_str(3, myio::Serialize(kernel_offset));
                    stmt_prod.step_and_reset();
                    dbProd.reg_transaction();
                }
            }
        }

        double time = timer.Elapsed();
        timer.Reset();
        std::cout << "t=" << deg.t << " s=" << deg.s << " time=" << time << std::endl;
        dbProd.save_time(table_out, deg.s, deg.t, time);

        dbProd.end_transaction(); ////
        diffs.erase(deg1);
    }
}

// std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
int main_2cell_prod(int argc, char** argv, int index)
{
    int t_max = 100;
    std::string complex_name = "C2";
    std::string db_in = "AdamsE2.db";
    std::string table_in = "SteenrodMRes";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  Adams complex [complex_name] [t_max] [db_in] [table_in]\n\n";

        std::cout << "complex_name=C2/Ceta/Cnu/Csigma\n\n";

        std::cout << "Default values:\n";
        std::cout << "  complex_name = " << complex_name << "\n";
        std::cout << "  t_max = " << t_max << "\n";
        std::cout << "  db_in = " << db_in << "\n";
        std::cout << "  table_in = " << table_in << "\n\n";

        std::cout << "Version:\n  2.0 (2022-08-17)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "complex_name", complex_name))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_in", db_in))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_in", table_in))
        return index;
    bench::Timer timer;

    compute_2cell_products_by_t(t_max, complex_name, db_in, table_in);

    return 0;
}

int main_2cell(int argc, char** argv, int index)
{
    std::string cmd;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Usage:\n  Adams 2cell <cmd> [-h] ...\n\n";

        std::cout << "<cmd> can be one of the following:\n";
        std::cout << "  prod: Compute the multiplications\n\n";

        std::cout << "Version:\n  2.0 (2022-08-07)" << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cmd", cmd))
        return index;

    if (cmd == "prod")
        return main_2cell_prod(argc, argv, index);
    else
        std::cerr << "Invalid cmd: " << argv[1] << '\n';

    return 0;
}
