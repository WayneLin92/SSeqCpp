#include "groebner_steenrod_const.h"
#include "algebras/database.h"
#include "algebras/linalg.h"
#include <cstring>

namespace steenrod {

/** Return i such that mon divides leads[i].
 *
 * Return -1 if not found.
 * @param indices indices[v_raw] stores the indices of elements of leads with the given v_raw.
 */
int IndexOfDivisibleLeading(const MMod1d& leads, const std::unordered_map<uint64_t, array>& indices, MMod mon)
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
 *                    class GroebnerMResConst
 ********************************************************/

GroebnerMResConst::GroebnerMResConst(int t_trunc, int s_trunc, DataMResConst2d data, array2d basis_degrees) : t_trunc_(t_trunc), s_trunc_(s_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees))
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

Mod GroebnerMResConst::DiffInv(Mod x, size_t s) const
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

class DbSteenrod : public myio::Database
{
    using Statement = myio::Statement;

private:
    std::future<void> f_;

public:
    DbSteenrod() = default;
    explicit DbSteenrod(const std::string& filename) : Database(filename) {}

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

    void create_products(const std::string& table_prefix, const GenMRes2d& gens) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators_products (id INTEGER PRIMARY KEY, s SMALLINT, t SMALLINT, indecomposable TINYINT);");

        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_generators_products (id, s, t, indecomposable) VALUES (?1, ?2, ?3, ?4);");

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

    void add_product_column(const std::string& table_prefix, int id) const
    {
        std::string id_str = std::to_string(id);
        int added = get_int("SELECT indecomposable FROM " + table_prefix + "_generators_products WHERE id=" + id_str + ";");
        if (!added) {
            execute_cmd("ALTER TABLE " + table_prefix + "_generators_products ADD COLUMN m" + id_str + " BLOB;");
            execute_cmd("ALTER TABLE " + table_prefix + "_generators_products ADD COLUMN mh" + id_str + " BLOB;");
            execute_cmd("UPDATE " + table_prefix + "_generators_products SET indecomposable=1 WHERE id=" + id_str + ";");
        }
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

    void save_products(const std::string& table_prefix, const array& ids_ind, const GenMRes2d& gens, const std::map<int, Mod2d>& map, const std::map<int, array3d>& map_h)
    {
        if (f_.valid())
            f_.wait();
        f_ = std::async(std::launch::async, [this, ids_ind, &gens, &map, &map_h, &table_prefix]() {
            begin_transaction();
            for (int id : ids_ind)
                add_product_column("SteenrodMRes", id);
            for (int id : ids_ind) {
                auto& map_id = map.at(id);
                auto& map_h_id = map_h.at(id);
                std::string id_str = std::to_string(id);
                Statement stmt(*this, "UPDATE OR IGNORE " + table_prefix + "_generators_products SET m" + id_str + "=?1, mh" + id_str + "=?2 WHERE id=?3 and m" + id_str + " IS NULL;");
                for (size_t s = 0; s < map_id.size(); ++s) {
                    for (size_t i = 0; i < map_id[s].size(); ++i) {
                        stmt.bind_blob(1, map_id[s][i].data);
                        stmt.bind_blob(2, map_h_id[s][i]);
                        stmt.bind_int(3, gens[s][i].id);
                        stmt.step_and_reset();
                    }
                }
            }
            end_transaction();
        });
    }

    void load_products(const std::string& table_prefix, std::map<int, Mod2d>& map, std::map<int, array3d>& map_h) const
    {
        array id_inds = get_column_int(table_prefix + "_generators_products", "id", "WHERE indecomposable=1 ORDER BY id");
        for (int id_ind : id_inds) {
            std::string id_ind_str = std::to_string(id_ind);
            Statement stmt(*this, "SELECT id, s, m" + id_ind_str + ", mh" + id_ind_str + " FROM " + table_prefix + "_generators_products WHERE m" + id_ind_str + " is not NULL ORDER BY id;");
            while (stmt.step() == MYSQLITE_ROW) {
                int id = stmt.column_int(0);
                int s = stmt.column_int(1);
                Mod prod;
                prod.data = stmt.column_blob_tpl<MMod>(2);
                array prod_h = stmt.column_blob_tpl<int>(3);
                if (map[id_ind].size() <= s) {
                    map[id_ind].resize(size_t(s + 1));
                    map_h[id_ind].resize(size_t(s + 1));
                }
                map.at(id_ind)[s].push_back(std::move(prod));
                map_h.at(id_ind)[s].push_back(std::move(prod_h));
            }
        }
    }

    array2d load_basis_degrees(const std::string& table_prefix) const
    {
        array2d result;
        Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int s = stmt.column_int(0), t = stmt.column_int(1);
            if (result.size() <= (size_t)s)
                result.resize(size_t(s + 1));
            result[s].push_back(t);
        }
        return result;
    }

    GenMRes2d load_generators(const std::string& table_prefix) const
    {
        GenMRes2d result;
        Statement stmt(*this, "SELECT id - 1, s, t, diff FROM " + table_prefix + "_generators ORDER BY id;"); //TODO: We use id-1 here because the original id starts with 1
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0), s = stmt.column_int(1), t = stmt.column_int(2);
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

GroebnerMResConst GroebnerMResConst::load(const std::string& filename)
{
    DbSteenrod db(filename);
    DataMResConst2d data = db.load_data("SteenrodMRes");
    array2d basis_degrees = db.load_basis_degrees("SteenrodMRes");
    int latest_s = 0, latest_t = 0;
    db.latest_st("SteenrodMRes", latest_s, latest_t);
    return GroebnerMResConst(latest_t, latest_s, std::move(data), std::move(basis_degrees));
}

array HomToK(const Mod& x)
{
    array result;
    for (MMod m : x.data)
        if (m.deg_m() == 0)
            result.push_back((int)m.v());
    return result;
}

/* Compute the products of gens with gens[s][k] */
void compute_products_sk(const GenMRes2d& gens, const GroebnerMResConst& gb, size_t s, size_t k, Mod2d& map, array3d& map_h)
{
    const size_t gens_size = gens.size();
    map.resize(gens_size);
    map_h.resize(gens_size);

    size_t old_map_s_size = map[s].size();

    map[s].resize(gens[s].size());
    map_h[s].resize(gens[s].size()); 

    map[s][k] = MMod(MMilnor(), 0);
    map_h[s][k] = HomToK(map[s][k]);

    Mod map_diff, tmp;
    for (size_t s1 = s + 1; s1 < gens_size; ++s1) {
        std::cout << "    s2=" << s1 << "          \r";
        size_t old_map_s1_size = map[s1].size();
        map[s1].resize(gens[s1].size());
        map_h[s1].resize(gens[s1].size());
        ut::for_each_par(gens[s1].size() - old_map_s1_size, [&](size_t i) {
            map[s1][old_map_s1_size + i] = gb.DiffInv(subs(gens[s1][old_map_s1_size + i].diff, map[s1 - 1]), s1 - s - 1);
            map_h[s1][old_map_s1_size + i] = HomToK(map[s1][old_map_s1_size + i]);
        });
    }
}

void compute_products()
{
    std::string filename = "C:\\Users\\lwnpk\\Documents\\Projects\\algtop_cpp_to_GZ_3\\AdamsE2_t100.db";
    std::string table_SteenrodMRes = "SteenrodMRes";
    auto gb = GroebnerMResConst::load(filename);
    DbSteenrod db(filename);
    auto gens = db.load_generators(table_SteenrodMRes);

    int num_gens = 0;
    for (size_t i = 0; i < gens.size(); ++i)
        num_gens += (int)gens[i].size();
    DbSteenrod dbProd("AdamsE2Prod.db");
    dbProd.create_products(table_SteenrodMRes, gens);

    std::map<int, Mod2d> map;
    std::map<int, array3d> map_h;
    dbProd.load_products(table_SteenrodMRes, map, map_h);

    array id_inds = dbProd.get_column_int(table_SteenrodMRes + "_generators_products", "id", "WHERE indecomposable=1 ORDER BY id");
    for (size_t s = 1; s < gens.size(); ++s) {
        std::cout << "s1=" << s << "          \n";
        array2d fx;
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
        array2d image = lina::GetSpace(fx);
        array r = ut::int_range((int)gens[s].size());
        array lead_image;
        for (const array& a : image)
            lead_image.push_back(a[0]);
        std::sort(lead_image.begin(), lead_image.end());
        array indices = lina::AddVectors(r, lead_image);
        for (size_t i = 0; i < gens[s].size(); ++i)
            if (std::binary_search(id_inds.begin(), id_inds.end(), gens[s][i].id))
                indices.push_back((int)i);
        std::sort(indices.begin(), indices.end());

        array ids_ind;
        for (size_t i = 0; i < indices.size(); ++i) {
            int id = gens[s][indices[i]].id;
            ids_ind.push_back(id);
            std::cout << "  i1=" << i << '/' << indices.size() << "          \n";
            compute_products_sk(gens, gb, s, indices[i], map[id], map_h[id]);
        }

        /* Save to database */
        dbProd.save_products(table_SteenrodMRes, ids_ind, gens, map, map_h);
    }
}

}  // namespace steenrod
