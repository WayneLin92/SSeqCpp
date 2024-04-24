#include "algebras/benchmark.h"
#include "algebras/dbAdamsSS.h"
#include "algebras/groebner.h"
#include "algebras/linalg.h"
#include "main.h"
#include <cstring>
#include <map>
#include <regex>

class MyDB : public myio::DbAdamsSS
{
    using Statement = myio::Statement;

public:
    MyDB() = default;
    explicit MyDB(const std::string& filename) : DbAdamsSS(filename)
    {
        // if (newFile_)
        SetVersion();
    }

    void SetVersion()
    {
        create_db_version(*this);
        Statement stmt(*this, "INSERT INTO version (id, name, value) VALUES (?1, ?2, ?3) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
        stmt.bind_and_step(0, std::string("version"), DB_ADAMS_VERSION);
        stmt.bind_and_step(1, std::string("change notes"), std::string(DB_VERSION_NOTES));
    }

public:
    void create_generators_cell(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, to_S0 TEXT, s SMALLINT, t SMALLINT);");
    }

    void drop_and_create_generators_cell(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_generators");
        create_generators_cell(table_prefix);
    }

    void save_generators_v2(const std::string& table_prefix, const alg2::AdamsDeg1d& gen_degs) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (id, s, t) VALUES (?1, ?2, ?3);");
        for (size_t i = 0; i < gen_degs.size(); ++i) {
            stmt.bind_and_step((int)i, gen_degs[i].s, gen_degs[i].t);
        }
    }

    void save_basis_mod_v2(const std::string& table_prefix, const std::map<alg2::AdamsDeg, alg2::MMod1d>& basis) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (id, mon, s, t) VALUES (?1, ?2, ?3, ?4);");

        int count = 0;
        auto degs = ut::get_keys(basis);
        std::sort(degs.begin(), degs.end(), [](const alg::AdamsDeg& d1, const alg::AdamsDeg& d2) { return d1.t < d2.t || (d1.t == d2.t && d1.s > d2.s); });
        for (auto& deg : degs) {
            auto& basis_d = basis.at(deg);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                stmt.bind_and_step(count, myio::Serialize(basis_d[i]), deg.s, deg.t);
                ++count;
            }
        }
    }

    void save_generators_cell(const std::string& table_prefix, const alg2::AdamsDeg1d& gen_degs, const alg2::Poly1d& toS0) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (id, to_S0, s, t) VALUES (?1, ?2, ?3, ?4);");

        for (size_t i = 0; i < gen_degs.size(); ++i) {
            stmt.bind_and_step((int)i, myio::Serialize(toS0[i]), gen_degs[i].s, gen_degs[i].t);
        }
    }

    void load_indecomposables(const std::string& table, alg2::int1d& ids, alg2::AdamsDeg1d& gen_degs, int t_trunc, int stem_trunc) const
    {
        Statement stmt(*this, fmt::format("SELECT id, s, t FROM {} WHERE indecomposable=1 AND t<={} AND t-s<={} ORDER BY t,-s;", table, t_trunc, stem_trunc));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0), s = stmt.column_int(1), t = stmt.column_int(2);
            ids.push_back(id);
            gen_degs.push_back(alg2::AdamsDeg{s, t});
        }
    }

    void load_indecomposables_cell(const std::string& table, alg2::int1d& ids, alg2::int2d& toS0res, alg2::AdamsDeg1d& gen_degs, int t_trunc) const
    {
        Statement stmt(*this, "SELECT id, coh_repr, s, t FROM " + table + " WHERE indecomposable=1 AND t<=" + std::to_string(t_trunc) + " ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0), s = stmt.column_int(2), t = stmt.column_int(3);
            if (id == 0)
                toS0res.push_back(alg2::int1d{});
            else
                toS0res.push_back(myio::Deserialize<alg2::int1d>(stmt.column_str(1)));
            ids.push_back(id);
            gen_degs.push_back(alg2::AdamsDeg{s, t});
        }
    }

    /* id_ind * id = prod */
    std::map<std::pair<int, int>, alg2::int1d> load_products_h(const std::string& table_prefix, int t_trunc, int stem_trunc) const
    {
        std::map<std::pair<int, int>, alg2::int1d> result;
        Statement stmt(*this, fmt::format("SELECT id, id_ind, prod_h FROM {}_products LEFT JOIN {}_generators USING(id) WHERE length(prod_h)>0 AND t<={} AND t-s<={} ORDER BY id;", table_prefix, table_prefix, t_trunc, stem_trunc));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            int g = stmt.column_int(1);
            alg2::int1d prod_h = myio::Deserialize<alg2::int1d>(stmt.column_str(2));
            result[std::make_pair(g, id)] = std::move(prod_h);
        }
        return result;
    }

    std::map<int, alg2::int1d> load_d2(const std::string& table) const
    {
        std::map<int, alg2::int1d> result;
        Statement stmt(*this, fmt::format("SELECT id, d2_h FROM {}", table));
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            auto d2_dual = myio::Deserialize<alg2::int1d>(stmt.column_str(1));
            int s = LocId(id).s - 2;
            for (auto& id1 : d2_dual)
                result[LocId(s, id1).id()].push_back(id);
        }
        return result;
    }

    /** id_ind * id = prod */
    std::map<std::pair<int, int>, alg2::int1d> load_cell_products_h(const std::string& table_prefix, int t_trunc) const
    {
        std::map<std::pair<int, int>, alg2::int1d> result;
        Statement stmt(*this, "SELECT id, id_ind, prod_h FROM " + table_prefix + "_products LEFT JOIN " + table_prefix + "_E2 USING(id) WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0);
            int id_ind = stmt.column_int(1);
            alg2::int1d prod_h = myio::Deserialize<alg2::int1d>(stmt.column_str(2));
            result[std::make_pair(id_ind, id)] = std::move(prod_h);
        }
        return result;
    }

    std::map<alg2::AdamsDeg, alg2::int2d> load_basis_repr(const std::string& table_prefix) const
    {
        std::map<alg2::AdamsDeg, alg2::int2d> result;
        Statement stmt(*this, "SELECT s, t, repr FROM " + table_prefix + "_basis ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            ++count;
            alg2::AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
            result[deg].push_back(myio::Deserialize<alg2::int1d>(stmt.column_str(2)));
        }
        // myio::Logger::out() << "repr loaded from " << table_prefix << "_basis, size=" << count << '\n';
        return result;
    }
};

class DbMapAdamsE2 : public myio::DbAdamsSS
{
    using Statement = myio::Statement;

public:
    DbMapAdamsE2() = default;
    explicit DbMapAdamsE2(const std::string& filename) : DbAdamsSS(filename)
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
    }

public:
    void recreate_tables(const std::string& table_prefix)
    {
        drop_table(table_prefix);
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + " (id INTEGER PRIMARY KEY, map TEXT);");
    }
    void save_map(const std::string& table_prefix, const alg2::Poly1d& map)
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + " (id, map) VALUES (?1, ?2);");
        for (size_t i = 0; i < map.size(); ++i)
            stmt.bind_and_step((int)i, myio::Serialize(map[i]));
    }
    void save_map(const std::string& table_prefix, const alg2::Mod1d& map)
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + " (id, map) VALUES (?1, ?2);");
        for (size_t i = 0; i < map.size(); ++i)
            stmt.bind_and_step((int)i, myio::Serialize(map[i]));
    }
};

alg2::int1d mul(const std::map<std::pair<int, int>, alg2::int1d>& map_h, int id_ind, const alg2::int1d& repr)
{
    alg2::int1d result;
    for (int id : repr) {
        auto it = map_h.find(std::make_pair(id_ind, id));
        if (it != map_h.end())
            result = lina::add(result, it->second);
    }
    return result;
}

void ExportRingAdamsE2(std::string_view ring, int t_trunc, int stem_trunc)
{
    using namespace alg2;

    std::string db_in = fmt::format("{}_Adams_res_prod.db", ring);
    std::string table_in = fmt::format("{}_Adams_res", ring);
    std::string db_out = fmt::format("{}_AdamsSS.db", ring);
    std::string table_out = fmt::format("{}_AdamsE2", ring);

    myio::AssertFileExists(db_in);
    MyDB dbProd(db_in);
    int t_max_prod = get_db_t_max(dbProd);
    if (t_trunc > t_max_prod) {
        t_trunc = t_max_prod;
        fmt::print("t_max is truncated to {}\n", t_max_prod);
    }

    AdamsDeg1d gen_degs;
    int1d gen_reprs;
    dbProd.load_indecomposables(table_in + "_generators", gen_reprs, gen_degs, t_trunc, stem_trunc);
    auto map_h_dual = dbProd.load_products_h(table_in, t_trunc, stem_trunc); /* (g, gx) -> x */
    std::map<std::pair<int, int>, int1d> map_h;
    for (auto& [p, arr] : map_h_dual) {
        int s_i = LocId(p.second).s - LocId(p.first).s;
        for (int i : arr)
            map_h[std::make_pair(p.first, LocId(s_i, i).id())].push_back(p.second);
    }

    std::map<AdamsDeg, Poly1d> gb;
    Mon2d leads;
    std::map<AdamsDeg, Mon1d> basis;
    std::map<AdamsDeg, int2d> repr;

    /* Add new basis */
    basis[AdamsDeg{0, 0}].push_back(Mon());
    repr[AdamsDeg{0, 0}].push_back({0});

    /* Add new basis */
    for (int t = 1; t <= t_trunc; t++) {
        std::map<AdamsDeg, Mon1d> basis_new;
        std::map<AdamsDeg, int2d> repr_new;

        /* Consider all possible basis in degree t */
        for (size_t gen_id = gen_degs.size(); gen_id-- > 0;) {
            auto repr_g = gen_reprs[gen_id];
            int t1 = t - gen_degs[gen_id].t;
            if (t1 >= 0) {
                auto p1 = basis.lower_bound(AdamsDeg{0, t1});
                auto p2 = basis.lower_bound(AdamsDeg{0, t1 + 1});
                for (auto p = p1; p != p2; ++p) {
                    auto deg_mon = p->first + gen_degs[gen_id];
                    if (deg_mon.stem() > stem_trunc)
                        continue;
                    for (size_t i = 0; i < p->second.size(); ++i) {
                        Mon& m = p->second[i];
                        auto& repr_m = repr.at(p->first)[i];
                        if (!m || (int)gen_id <= m[0].g()) {
                            Mon mon = m * Mon::Gen((uint32_t)gen_id);
                            auto repr_mon = mul(map_h, repr_g, repr_m);
                            if (gen_id >= leads.size() || std::none_of(leads[gen_id].begin(), leads[gen_id].end(), [&mon](const Mon& _m) { return divisible(_m, mon); })) {
                                basis_new[deg_mon].push_back(std::move(mon));
                                repr_new[deg_mon].push_back(std::move(repr_mon));
                            }
                        }
                    }
                }
            }
        }

        /* Compute groebner and basis in degree t */
        for (auto it = basis_new.begin(); it != basis_new.end(); ++it) {
            auto indices = ut::size_t_range(it->second.size());
            std::sort(indices.begin(), indices.end(), [&it](size_t a, size_t b) { return it->second[a] < it->second[b]; });
            Mon1d basis_new_d;
            int2d repr_new_d;
            for (size_t i = 0; i < indices.size(); ++i) {
                basis_new_d.push_back(it->second[indices[i]]);
                repr_new_d.push_back(repr_new.at(it->first)[indices[i]]);
            }

            int2d image, kernel, g;
            lina::SetLinearMap(repr_new_d, image, kernel, g);
            int1d lead_kernel;

            /* Add to groebner */
            for (const int1d& k : kernel) {
                lead_kernel.push_back(k[0]);
                gb[it->first].push_back(Indices2Poly(k, basis_new_d));
                auto& poly = gb.at(it->first).back();
                int index = poly.GetLead().front().g();
                if (leads.size() <= (size_t)index)
                    leads.resize(size_t(index + 1));
                leads[index].push_back(poly.GetLead());
            }

            /* Add to basis_H */
            std::sort(lead_kernel.begin(), lead_kernel.end());
            int1d index_basis = lina::add(ut::int_range(int(basis_new_d.size())), lead_kernel);
            for (int i : index_basis) {
                basis[it->first].push_back(std::move(basis_new_d[i]));
                repr[it->first].push_back(std::move(repr_new_d[i]));
            }
        }
    }

    /* Check basis size */
    size_t size_basis = 0;
    for (auto it = basis.begin(); it != basis.end(); ++it)
        size_basis += it->second.size();
    int size_E2_res = dbProd.get_int(fmt::format("SELECT count(*) from {}_generators WHERE t<={} and t-s<={}", table_in, t_trunc, stem_trunc));
    if (size_basis != size_E2_res) {
        fmt::print("size_basis={}\n", size_basis);
        fmt::print("size_E2_res={}\n", size_E2_res);
        fmt::print("Error: sizes of bases do not match.");
        throw MyException(0x8334a2c1U, "Error: sizes of bases do not match.");
    }

    /* Save the results */
    MyDB dbE2(db_out);
    dbE2.drop_and_create_generators(table_out);
    dbE2.drop_and_create_relations(table_out);
    dbE2.drop_and_create_basis(table_out);

    dbE2.begin_transaction();

    dbE2.save_generators(table_out, gen_degs, gen_reprs);
    dbE2.save_gb(table_out, gb);
    dbE2.save_basis(table_out, basis, repr);
    set_db_t_max(dbE2, t_trunc);
    set_db_time(dbE2);

    dbE2.end_transaction();
}

void ExportModAdamsE2(std::string_view mod, std::string_view ring, int t_trunc, int stem_trunc)
{
    using namespace alg2;

    std::string db_mod = fmt::format("{}_Adams_res_prod.db", mod);
    std::string table_mod = fmt::format("{}_over_{}_Adams_res", mod, ring);
    const std::string db_ring = fmt::format("{}_AdamsSS.db", ring);
    const std::string table_ring = fmt::format("{}_AdamsE2", ring);
    myio::AssertFileExists(db_mod);
    myio::AssertFileExists(db_ring);

    std::string db_out = fmt::format("{}_AdamsSS.db", mod);
    std::string table_out = fmt::format("{}_AdamsE2", mod);

    MyDB dbProd(db_mod);
    int t_max_prod = get_db_t_max(dbProd);
    MyDB dbRing(db_ring);
    int t_max_ring = get_db_t_max(dbRing);
    int t_max_out = std::min(t_max_prod, t_max_ring);
    if (t_trunc > t_max_out) {
        t_trunc = t_max_out;
        fmt::print("t_max is truncated to {}\n", t_max_out);
    }
    if (!dbProd.has_table(table_mod + "_generators")) {
        constexpr std::array<std::string_view, 3> postfixes = {"generators", "products", "products_time"};
        for (auto& postfix : postfixes)
            dbProd.rename_table(fmt::format("{}_Adams_res_{}", mod, postfix), fmt::format("{}_over_{}_Adams_res_{}", mod, ring, postfix));
        fmt::print("Renamming success!\n");
    }
    AdamsDeg1d v_degs;
    int1d gen_reprs;
    dbProd.load_indecomposables(table_mod + "_generators", gen_reprs, v_degs, t_trunc, stem_trunc);
    auto map_h_dual = dbProd.load_products_h(table_mod, t_trunc, stem_trunc);
    std::map<std::pair<int, int>, int1d> map_h;
    for (auto& [p, arr] : map_h_dual) {
        int s_i = LocId(p.second).s - LocId(p.first).s;
        for (int i : arr)
            map_h[std::make_pair(p.first, LocId(s_i, i).id())].push_back(p.second);
    }

    auto ring_basis = dbRing.load_basis(table_ring);
    auto ring_basis_repr = dbRing.load_basis_repr(table_ring);

    std::map<AdamsDeg, Mod1d> gbm;
    Mon2d leads;
    std::map<AdamsDeg, MMod1d> basis;
    std::map<AdamsDeg, int2d> repr;

    /* Add new basis */
    for (int t = 0; t <= t_trunc; ++t) {
        std::map<AdamsDeg, MMod1d> basis_new;
        std::map<AdamsDeg, int2d> repr_new;

        /* Consider all possible basis in degree t */
        for (size_t gen_id = v_degs.size(); gen_id-- > 0;) {
            auto repr_g = gen_reprs[gen_id];
            int t1 = t - v_degs[gen_id].t;
            if (t1 >= 0) {
                auto p1 = ring_basis.lower_bound(AdamsDeg{0, t1});
                auto p2 = ring_basis.lower_bound(AdamsDeg{0, t1 + 1});
                for (auto p = p1; p != p2; ++p) {
                    auto deg_mon = p->first + v_degs[gen_id];
                    if (deg_mon.stem() > stem_trunc)
                        continue;
                    for (size_t i = 0; i < p->second.size(); ++i) {
                        MMod m = MMod(p->second[i], (uint32_t)gen_id);
                        if (size_t(gen_id) >= leads.size() || std::none_of(leads[gen_id].begin(), leads[gen_id].end(), [&m](const Mon& _m) { return divisible(_m, m.m); })) {
                            basis_new[deg_mon].push_back(m);
                            auto repr_mon = mul(map_h, repr_g, ring_basis_repr.at(p->first)[i]);
                            repr_new[deg_mon].push_back(std::move(repr_mon));
                        }
                    }
                }
            }
        }

        /* Compute groebner and basis in degree t */
        for (auto it = basis_new.begin(); it != basis_new.end(); ++it) {
            auto indices = ut::size_t_range(it->second.size());
            std::sort(indices.begin(), indices.end(), [&it](size_t a, size_t b) { return it->second[a] < it->second[b]; });
            MMod1d basis_new_d;
            int2d repr_new_d;
            for (size_t i = 0; i < indices.size(); ++i) {
                basis_new_d.push_back(it->second[indices[i]]);
                repr_new_d.push_back(repr_new.at(it->first)[indices[i]]);
            }

            int2d image, kernel, g;
            lina::SetLinearMap(repr_new_d, image, kernel, g);
            int1d lead_kernel;

            /* Add to groebner */
            for (const int1d& k : kernel) {
                lead_kernel.push_back(k[0]);
                gbm[it->first].push_back(Indices2Mod(k, basis_new_d));
                auto& x = gbm.at(it->first).back();
                size_t index = (size_t)x.GetLead().v;
                if (leads.size() <= index)
                    leads.resize(index + 1);
                leads[index].push_back(x.GetLead().m);
            }

            /* Add to basis_H */
            std::sort(lead_kernel.begin(), lead_kernel.end());
            int1d index_basis = lina::add(ut::int_range(int(basis_new_d.size())), lead_kernel);
            for (int i : index_basis) {
                basis[it->first].push_back(std::move(basis_new_d[i]));
                repr[it->first].push_back(std::move(repr_new_d[i]));
            }
        }
    }

    size_t size_basis = 0;
    for (auto it = basis.begin(); it != basis.end(); ++it) {
        size_basis += it->second.size();
    }
    int size_E2 = dbProd.get_int(fmt::format("SELECT count(*) FROM {}_generators WHERE t<={} AND t-s<={}", table_mod, t_trunc, stem_trunc));
    if (size_basis != size_t(size_E2)) {
        std::cout << "size_basis=" << size_basis << '\n';
        std::cout << "size_E2=" << size_E2 << '\n';
        throw MyException(0x571c8119U, "Error: sizes of bases do not match.");
    }

    /* Save the results */
    MyDB dbE2(db_out);
    dbE2.begin_transaction();

    dbE2.drop_and_create_generators(table_out);
    dbE2.drop_and_create_relations(table_out);
    dbE2.drop_and_create_basis(table_out);

    dbE2.save_generators(table_out, v_degs, gen_reprs);
    dbE2.save_gb_mod(table_out, gbm);
    dbE2.save_basis_mod(table_out, basis, repr);
    set_db_t_max(dbE2, t_trunc);
    set_db_time(dbE2);
    set_db_over(dbE2, std::string(ring));

    dbE2.end_transaction();
}

void ExportAdamsD2(std::string_view cw)
{
    using namespace alg2;

    std::string db_E2 = fmt::format("{}_AdamsSS.db", cw);
    std::string table_E2 = fmt::format("{}_AdamsE2", cw);
    std::string db_d2 = fmt::format("{}_Adams_d2.db", cw);
    std::string table_d2 = fmt::format("{}_Adams_d2", cw);
    myio::AssertFileExists(db_E2);
    myio::AssertFileExists(db_d2);

    MyDB dbE2(db_E2);
    int t_max_cw = get_db_t_max(dbE2);
    MyDB dbD2(db_d2);
    int t_max_d2 = get_db_t_max(dbD2);
    int t_max_out = std::min(t_max_cw, t_max_d2);

    auto d2E2 = dbD2.load_d2(table_d2);

    std::map<AdamsDeg, int2d> basis_to_res_id;
    {
        myio::Statement stmt(dbE2, fmt::format("select s, t, repr from {}_basis order by id;", table_E2));
        while (stmt.step() == MYSQLITE_ROW) {
            int s = stmt.column_int(0), t = stmt.column_int(1);
            AdamsDeg deg = {s, t};
            int1d repr = myio::Deserialize<int1d>(stmt.column_str(2));
            basis_to_res_id[deg].push_back(std::move(repr));
        }
    }
    std::map<AdamsDeg, int> min_res_id;
    for (auto& [deg, map] : basis_to_res_id) {
        int1d min_res_id1d;
        for (auto& map_i : map)
            if (!map_i.empty())
                min_res_id1d.push_back(*std::min_element(map_i.begin(), map_i.end()));
        min_res_id[deg] = *std::min_element(min_res_id1d.begin(), min_res_id1d.end());
    }
    std::map<int, int1d> res_id_to_basis;
    for (auto& [deg, map] : basis_to_res_id) {
        for (int i = 0; i < map.size(); ++i) {
            int2d image, kernel, g;
            lina::SetLinearMap(map, image, kernel, g);
            int res_id = min_res_id[deg] + i;
            res_id_to_basis[res_id] = lina::GetImage(image, g, {res_id});
        }
    }

    dbE2.begin_transaction();
    if (!dbE2.has_column(table_E2 + "_basis", "d2"))
        dbE2.add_column(table_E2 + "_basis", "d2 TEXT");

    std::vector<std::pair<int, int1d>> monid_repr;
    {
        myio::Statement stmt(dbE2, fmt::format("select id, repr from {}_AdamsE2_basis WHERE t<={} order by id", cw, t_max_out - 1));
        while (stmt.step() == MYSQLITE_ROW) {
            monid_repr.push_back(std::make_pair(stmt.column_int(0), myio::Deserialize<int1d>(stmt.column_str(1))));
        }
    }
    {
        myio::Statement stmt(dbE2, fmt::format("UPDATE {}_AdamsE2_basis SET d2=?1 WHERE id=?2;", cw));
        for (auto& [id, repr] : monid_repr) {
            int1d d2;
            for (int i : repr)
                if (ut::has(d2E2, i))
                    for (int id_d2 : d2E2.at(i))
                        d2 = lina::add(d2, res_id_to_basis.at(id_d2));
            stmt.bind_and_step(myio::Serialize(d2), id);
        }
    }

    set_db_d2_t_max(dbE2, t_max_out);
    set_db_time(dbE2);

    dbE2.end_transaction();
}

void ExportFreeModAdamsE2(std::string_view mod, std::string_view ring, const alg2::int1d& cells, int t_trunc, int stem_trunc)
{
    using namespace alg2;

    const std::string db_ring = fmt::format("{}_AdamsSS.db", ring);
    const std::string table_ring = fmt::format("{}_AdamsE2", ring);
    myio::AssertFileExists(db_ring);

    std::string db_out = fmt::format("{}_AdamsSS.db", mod);
    std::string table_out = fmt::format("{}_AdamsE2", mod);

    MyDB dbRing(db_ring);
    int t_max_ring = get_db_t_max(dbRing);
    if (t_trunc > t_max_ring) {
        t_trunc = t_max_ring;
        fmt::print("t_max is truncated to {}\n", t_max_ring);
    }

    AdamsDeg1d v_degs;
    for (int i : cells)
        v_degs.push_back(AdamsDeg(0, i));

    auto ring_basis = dbRing.load_basis(table_ring);
    auto ring_basis_repr = dbRing.load_basis_repr(table_ring);

    std::map<AdamsDeg, MMod1d> basis;

    /* Add new basis */
    for (int t = 0; t <= t_trunc; ++t) {
        /* Consider all possible basis in degree t */
        for (size_t gen_id = v_degs.size(); gen_id-- > 0;) {
            int t1 = t - v_degs[gen_id].t;
            if (t1 >= 0) {
                auto p1 = ring_basis.lower_bound(AdamsDeg{0, t1});
                auto p2 = ring_basis.lower_bound(AdamsDeg{0, t1 + 1});
                for (auto p = p1; p != p2; ++p) {
                    auto deg_mon = p->first + v_degs[gen_id];
                    if (deg_mon.stem() > stem_trunc)
                        continue;
                    for (size_t i = 0; i < p->second.size(); ++i) {
                        MMod m = MMod(p->second[i], (uint32_t)gen_id);
                        basis[deg_mon].push_back(m);
                    }
                }
            }
        }
    }

    /* Save the results */
    MyDB dbE2(db_out);
    dbE2.begin_transaction();

    dbE2.drop_and_create_generators(table_out);
    dbE2.drop_and_create_relations(table_out);
    dbE2.drop_and_create_basis(table_out);

    dbE2.save_generators_v2(table_out, v_degs);
    dbE2.save_gb_mod(table_out, {});
    dbE2.save_basis_mod_v2(table_out, basis);
    set_db_t_max(dbE2, t_trunc);
    set_db_time(dbE2);

    dbE2.end_transaction();
}

void load_map(const myio::Database& db, std::string_view table_prefix, std::map<int, alg2::int1d>& map_h, int fil)
{
    using namespace alg2;
    std::map<int, alg2::int1d> map_h_dual;
    myio::Statement stmt(db, fmt::format("SELECT id, map_h FROM {} ORDER BY id", table_prefix));
    while (stmt.step() == MYSQLITE_ROW) {
        int id = stmt.column_int(0);
        int1d map_ = myio::Deserialize<int1d>(stmt.column_str(1));
        map_h_dual[id] = std::move(map_);
    }
    for (auto& [i, arr] : map_h_dual) {
        int s = LocId(i).s - fil;
        for (int v : arr)
            map_h[LocId(s, v).id()].push_back(i);
    }
}

void ExportMapAdamsE2(std::string_view cw1, std::string_view cw2, int t_trunc)
{
    using namespace alg2;
    /* A sorted list of ring spectra */
    std::vector<std::string> RING_SPECTRA = {"S0", "tmf", "ko", "X2", "A_A3", "A_A4", "A_A5"};
    std::sort(RING_SPECTRA.begin(), RING_SPECTRA.end());

    std::string db_map = fmt::format("map_Adams_res_{}__{}.db", cw1, cw2);
    std::string table_map = fmt::format("map_Adams_res_{}__{}", cw1, cw2);
    myio::AssertFileExists(db_map);
    myio::Database dbResMap(db_map);
    int t_max_map = get_db_t_max(dbResMap);

    auto from = dbResMap.get_str("select value from version where id=446174262");
    auto to = dbResMap.get_str("select value from version where id=1713085477");
    int fil = 0;
    try { /* For compatibility */
        fil = dbResMap.get_int("select value from version where id=651971502");
    }
    catch (MyException&) {
    }
    int sus = dbResMap.get_int("select value from version where id=1585932889"); /* cw1->Sigma^sus cw2 */

    const std::string db_cw1 = fmt::format("{}_AdamsSS.db", from);
    const std::string table_cw1 = fmt::format("{}_AdamsE2", from);
    const std::string db_cw2 = fmt::format("{}_AdamsSS.db", to);
    const std::string table_cw2 = fmt::format("{}_AdamsE2", to);
    myio::AssertFileExists(db_cw1);
    myio::AssertFileExists(db_cw2);

    MyDB dbCw1(db_cw1);
    int t_max_cw1 = get_db_t_max(dbCw1);
    MyDB dbCw2(db_cw2);
    int t_max_cw2 = get_db_t_max(dbCw2);
    int t_max_out = std::min({t_max_map, t_max_cw1, t_max_cw2 + sus - fil});
    if (t_trunc > t_max_out) {
        t_trunc = t_max_out;
        fmt::print("t_max is truncated to {}\n", t_max_out);
    }

    auto gen_degs = dbCw1.load_gen_adamsdegs(table_cw1);
    auto gen_reprs1 = dbCw1.get_column_int(table_cw1 + "_generators", "repr", "");
    int1d out_of_region;
    for (size_t i = 0; i < gen_degs.size(); ++i)
        if (gen_degs[i].t > t_trunc /*|| gen_degs[i].stem()*/)
            out_of_region.push_back((int)i);

    std::map<int, int1d> map_h;
    load_map(dbResMap, table_map, map_h, fil);

    std::string db_out = fmt::format("map_AdamsSS_{}__{}.db", cw1, cw2);
    std::string table_out = fmt::format("map_AdamsE2_{}__{}", cw1, cw2);
    DbMapAdamsE2 dbMapAdamsE2(db_out);
    dbMapAdamsE2.recreate_tables(table_out);
    dbMapAdamsE2.begin_transaction();

    if (ut::has(RING_SPECTRA, to)) {
        auto basis2 = dbCw2.load_basis(table_cw2);
        auto basis_repr2 = dbCw2.load_basis_repr(table_cw2);
        Poly1d images;
        for (size_t i = 0; i < gen_reprs1.size(); ++i) {
            if (ut::has(out_of_region, (int)i))
                images.push_back(Poly::Gen(uint32_t(-1)));
            else {
                AdamsDeg d = gen_degs[i] + AdamsDeg(fil, fil - sus);
                if (!ut::has(basis2, d) || !ut::has(map_h, gen_reprs1[i])) {
                    images.push_back({});
                    continue;
                }
                const auto& basis2_d = basis2.at(d);
                const int2d& r = basis_repr2.at(d);
                int2d image, _k, g;
                lina::SetLinearMap(r, image, _k, g);
                images.push_back(Indices2Poly(lina::GetImage(image, g, map_h.at(gen_reprs1[i])), basis2_d));
            }
        }
        dbMapAdamsE2.save_map(table_out, images);
    }
    else {
        auto basis2 = dbCw2.load_basis_mod(table_cw2);
        auto basis_repr2 = dbCw2.load_basis_repr(table_cw2);
        Mod1d images;
        for (size_t i = 0; i < gen_reprs1.size(); ++i) {
            if (ut::has(out_of_region, (int)i))
                images.push_back(MMod(Mon(), uint32_t(-1)));
            else {
                AdamsDeg d = gen_degs[i] + AdamsDeg(fil, fil - sus);
                if (!ut::has(basis2, d) || !ut::has(map_h, gen_reprs1[i])) {
                    images.push_back({});
                    continue;
                }
                const auto& basis2_d = basis2.at(d);
                const int2d& r = basis_repr2.at(d);
                int2d image, _k, g;
                lina::SetLinearMap(r, image, _k, g);
                images.push_back(Indices2Mod(lina::GetImage(image, g, map_h.at(gen_reprs1[i])), basis2_d));
            }
        }
        dbMapAdamsE2.save_map(table_out, images);
    }
    {
        myio::Statement stmt(dbMapAdamsE2, "INSERT INTO version (id, name, value) VALUES (?1, ?2, ?3) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
        stmt.bind_and_step(651971502, std::string("filtration"), fil);
        stmt.bind_and_step(1585932889, std::string("suspension"), sus);
        stmt.bind_and_step(446174262, std::string("from"), from);
        stmt.bind_and_step(1713085477, std::string("to"), to);
    }
    set_db_t_max(dbMapAdamsE2, t_trunc);
    set_db_time(dbMapAdamsE2);
    dbMapAdamsE2.end_transaction();
}

void ExportMapSumAdamsE2(const std::string& cw1, const std::string& cw2, const std::string& db_map1, const std::string& db_map2, int v0, int v1, int sus, int fil, int t_trunc)
{
    using namespace alg2;
    std::regex is_map_regex("^map_AdamsSS_(\\w+?__\\w+?)(?:_t\\d+|).db$"); /* match example: map_AdamsSS_RP1_4__RP3_4_t169.db */
    std::smatch match;

    int t_max_map1 = 500, t_max_map2 = 500;
    Poly1d images1, images2;

    if (!db_map1.empty()) {
        myio::AssertFileExists(db_map1);
        myio::Database dbMap1(db_map1);
        t_max_map1 = get_db_t_max(dbMap1);

        std::string table1;
        if (std::regex_search(db_map1, match, is_map_regex); match[0].matched)
            table1 = fmt::format("map_AdamsE2_{}", match[1].str());
        else {
            fmt::print("filename={} not supported.\n", db_map1);
            throw MyException(0x248b04b, "File name is not supported.");
        }
        images1 = dbMap1.get_column_from_str<Poly>(table1, "map", "ORDER BY id", myio::Deserialize<Poly>);
    }

    if (!db_map2.empty()) {
        myio::AssertFileExists(db_map2);
        myio::Database dbMap2(db_map2);
        t_max_map2 = get_db_t_max(dbMap2);

        std::string table2;
        if (std::regex_search(db_map2, match, is_map_regex); match[0].matched)
            table2 = fmt::format("map_AdamsE2_{}", match[1].str());
        else {
            fmt::print("filename={} not supported.\n", db_map2);
            throw MyException(0xa833c7c4, "File name is not supported.");
        }
        images2 = dbMap2.get_column_from_str<Poly>(table2, "map", "ORDER BY id", myio::Deserialize<Poly>);
    }

    if (images1.empty())
        images1.resize(images2.size());
    else if (images2.empty())
        images2.resize(images1.size());
    Mod1d images3;
    for (size_t i = 0; i < images1.size(); ++i)
        images3.push_back(Mod(images1[i], v0) + Mod(images2[i], v1));

    std::string db_out = fmt::format("map_AdamsSS_{}__{}.db", cw1, cw2);
    std::string table_out = fmt::format("map_AdamsE2_{}__{}", cw1, cw2);
    DbMapAdamsE2 dbMapAdamsE2(db_out);
    dbMapAdamsE2.recreate_tables(table_out);
    dbMapAdamsE2.begin_transaction();
    dbMapAdamsE2.save_map(table_out, images3);

    {
        myio::Statement stmt(dbMapAdamsE2, "INSERT INTO version (id, name, value) VALUES (?1, ?2, ?3) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
        stmt.bind_and_step(1585932889, std::string("suspension"), sus);
        stmt.bind_and_step(651971502, std::string("filtration"), fil);
        stmt.bind_and_step(446174262, std::string("from"), cw1);
        stmt.bind_and_step(1713085477, std::string("to"), cw2);
    }
    int t_max_out = std::min({t_max_map1, t_max_map2});
    if (t_trunc > t_max_out) {
        t_trunc = t_max_out;
        fmt::print("t_max is truncated to {}\n", t_max_out);
    }
    set_db_t_max(dbMapAdamsE2, t_trunc);
    set_db_time(dbMapAdamsE2);
    dbMapAdamsE2.end_transaction();
}

void ExportMapFromFreeModAdamsE2(const std::string& cw1, const std::string& cw2, const nlohmann::json& map_json, int t_trunc)
{
    using namespace alg2;
    std::vector<std::string> RING_SPECTRA = {"S0", "tmf", "ko", "X2"};
    std::sort(RING_SPECTRA.begin(), RING_SPECTRA.end());

    const auto to = map_json.at("to").get<std::string>();
    const auto db_cw1 = fmt::format("{}_AdamsSS.db", cw1);
    const auto db_cw2 = fmt::format("{}_AdamsSS.db", to);
    const auto table_cw2 = fmt::format("{}_AdamsE2", to);
    myio::AssertFileExists(db_cw1);
    myio::AssertFileExists(db_cw2);

    AdamsDeg1d degs;
    int1d indices;
    for (auto& v : map_json.at("free_images")) {
        if (v.empty()) {
            degs.push_back(AdamsDeg(-1, -1));
            indices.push_back(-1);
        }
        else {
            int stem = v[0].get<int>();
            int s = v[1].get<int>();
            int i = v[2].get<int>();
            degs.push_back(AdamsDeg(s, stem + s));
            indices.push_back(i);
        }
    }

    MyDB dbCw1(db_cw1);
    MyDB dbCw2(db_cw2);
    int t_max1 = get_db_t_max(dbCw1);
    int t_max2 = get_db_t_max(dbCw2) - degs.front().t;
    int t_max_out = std::min(t_max1, t_max2);
    if (t_trunc > t_max_out) {
        t_trunc = t_max_out;
        fmt::print("t_max is truncated to {}\n", t_max_out);
    }

    Poly1d images;
    Mod1d images_mod;
    if (ut::has(RING_SPECTRA, to)) {
        auto basis2 = dbCw2.load_basis(table_cw2);
        for (size_t i = 0; i < degs.size(); ++i) {
            if (degs[i] != AdamsDeg(-1, -1))
                images.push_back(basis2[degs[i]][indices[i]]);
            else
                images.push_back({});
        }
    }
    else {
        auto basis2 = dbCw2.load_basis_mod(table_cw2);
        for (size_t i = 0; i < degs.size(); ++i) {
            if (degs[i] != AdamsDeg(-1, -1))
                images_mod.push_back(basis2[degs[i]][indices[i]]);
            else
                images_mod.push_back({});
        }
    }

    std::string db_out = fmt::format("map_AdamsSS_{}__{}.db", cw1, cw2);
    std::string table_out = fmt::format("map_AdamsE2_{}__{}", cw1, cw2);
    DbMapAdamsE2 dbMapAdamsE2(db_out);
    dbMapAdamsE2.recreate_tables(table_out);
    dbMapAdamsE2.begin_transaction();
    if (ut::has(RING_SPECTRA, to))
        dbMapAdamsE2.save_map(table_out, images);
    else
        dbMapAdamsE2.save_map(table_out, images_mod);
    {
        myio::Statement stmt(dbMapAdamsE2, "INSERT INTO version (id, name, value) VALUES (?1, ?2, ?3) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
        int sus = myio::get(map_json, "sus", 0);
        int fil = myio::get(map_json, "fil", 0);
        stmt.bind_and_step(1585932889, std::string("suspension"), sus);
        stmt.bind_and_step(651971502, std::string("filtration"), fil);
        stmt.bind_and_step(446174262, std::string("from"), cw1);
        stmt.bind_and_step(1713085477, std::string("to"), to);
    }
    set_db_t_max(dbMapAdamsE2, t_trunc);
    set_db_time(dbMapAdamsE2);
    dbMapAdamsE2.end_transaction();
}

void Export2Cell(std::string_view cw, std::string_view ring, int t_trunc, int stem_trunc)
{
    using namespace alg2;

    int t_cell = 0;
    if (cw == "C2")
        t_cell = 1;
    else if (cw == "Ceta")
        t_cell = 2;
    else if (cw == "Cnu")
        t_cell = 4;
    else if (cw == "Csigma")
        t_cell = 8;
    else {
        fmt::print("cw={} is not supported.\n", cw);
        return;
    }

    std::string db_in = fmt::format("{}_Adams_chain.db", cw);
    std::string table_in = fmt::format("{}_Adams", cw);
    if (ring != "S0") {
        db_in = fmt::format("{}_{}", ring, db_in);
        table_in = fmt::format("{}_{}", ring, table_in);
    }
    const std::string db_S0 = fmt::format("{}_AdamsSS.db", ring);
    const std::string table_S0 = fmt::format("{}_AdamsE2", ring);
    myio::AssertFileExists(db_in);
    myio::AssertFileExists(db_S0);

    std::string db_out = fmt::format("{}_AdamsSS.db", cw);
    std::string table_out = fmt::format("{}_AdamsE2", cw);
    if (ring != "S0") {
        db_out = fmt::format("{}_{}", ring, db_out);
        table_out = fmt::format("{}_{}", ring, table_out);
    }

    MyDB dbProd(db_in);
    AdamsDeg1d v_degs;
    int1d gen_reprs;
    int2d toS0res;
    dbProd.load_indecomposables_cell(table_in + "_E2", gen_reprs, toS0res, v_degs, t_trunc);
    auto map_h_dual = dbProd.load_cell_products_h(table_in, t_trunc);
    std::map<std::pair<int, int>, int1d> map_h;
    for (auto& [p, arr] : map_h_dual) {
        int id_ind_mod = p.first;
        int id_mod = p.second;
        for (int i : arr)
            map_h[std::make_pair(id_ind_mod, i)].push_back(id_mod);
    }

    MyDB dbS0(db_S0);
    auto S0_basis = dbS0.load_basis(table_S0);
    auto S0_basis_repr = dbS0.load_basis_repr(table_S0);

    size_t gen_size = gen_reprs.size();
    Poly1d toS0(gen_size);
    for (size_t i = 1; i < gen_size; ++i) {
        AdamsDeg d = v_degs[i] - AdamsDeg(0, t_cell);
        const Mon1d& b = S0_basis.at(d);
        const int2d& r = S0_basis_repr.at(d);
        int2d image, _k, g;
        lina::SetLinearMap(r, image, _k, g);
        toS0[i] = Indices2Poly(lina::GetImage(image, g, toS0res[i]), b);
    }

    std::map<AdamsDeg, Mod1d> gbm;
    Mon2d leads;
    std::map<AdamsDeg, MMod1d> basis;
    std::map<AdamsDeg, int2d> repr;

    /* Add new basis */
    for (int t = 0; t <= t_trunc; ++t) {
        std::map<AdamsDeg, MMod1d> basis_new;
        std::map<AdamsDeg, int2d> repr_new;

        /* Consider all possible basis in degree t */
        for (size_t gen_id = v_degs.size(); gen_id-- > 0;) {
            auto repr_g = gen_reprs[gen_id];
            int t1 = t - v_degs[gen_id].t;
            if (t1 >= 0) {
                auto p1 = S0_basis.lower_bound(AdamsDeg{0, t1});
                auto p2 = S0_basis.lower_bound(AdamsDeg{0, t1 + 1});
                for (auto p = p1; p != p2; ++p) {
                    auto deg_mon = p->first + v_degs[gen_id];
                    if (deg_mon.stem() > stem_trunc)
                        continue;
                    for (size_t i = 0; i < p->second.size(); ++i) {
                        MMod m = MMod(p->second[i], (uint32_t)gen_id);
                        if (size_t(gen_id) >= leads.size() || std::none_of(leads[gen_id].begin(), leads[gen_id].end(), [&m](const Mon& _m) { return divisible(_m, m.m); })) {
                            basis_new[deg_mon].push_back(m);
                            auto repr_mon = mul(map_h, repr_g, S0_basis_repr.at(p->first)[i]);
                            repr_new[deg_mon].push_back(std::move(repr_mon));
                        }
                    }
                }
            }
        }

        /* Compute groebner and basis in degree t */
        for (auto it = basis_new.begin(); it != basis_new.end(); ++it) {
            auto indices = ut::size_t_range(it->second.size());
            std::sort(indices.begin(), indices.end(), [&it](size_t a, size_t b) { return it->second[a] < it->second[b]; });
            MMod1d basis_new_d;
            int2d repr_new_d;
            for (size_t i = 0; i < indices.size(); ++i) {
                basis_new_d.push_back(it->second[indices[i]]);
                repr_new_d.push_back(repr_new.at(it->first)[indices[i]]);
            }

            int2d image, kernel, g;
            lina::SetLinearMap(repr_new_d, image, kernel, g);
            int1d lead_kernel;

            /* Add to groebner */
            for (const int1d& k : kernel) {
                lead_kernel.push_back(k[0]);
                gbm[it->first].push_back(Indices2Mod(k, basis_new_d));
                auto& x = gbm.at(it->first).back();
                size_t index = (size_t)x.GetLead().v;
                if (leads.size() <= index)
                    leads.resize(index + 1);
                leads[index].push_back(x.GetLead().m);
            }

            /* Add to basis_H */
            std::sort(lead_kernel.begin(), lead_kernel.end());
            int1d index_basis = lina::add(ut::int_range(int(basis_new_d.size())), lead_kernel);
            for (int i : index_basis) {
                basis[it->first].push_back(std::move(basis_new_d[i]));
                repr[it->first].push_back(std::move(repr_new_d[i]));
            }
        }
    }

    size_t size_basis = 0;
    for (auto it = basis.begin(); it != basis.end(); ++it) {
        size_basis += it->second.size();
    }
    int size_E2 = dbProd.get_int(fmt::format("SELECT count(*) FROM {}_E2 WHERE t<={} AND t-s<={}", table_in, t_trunc, stem_trunc));
    if (size_basis != size_t(size_E2)) {
        std::cout << "size_basis=" << size_basis << '\n';
        std::cout << "size_E2=" << size_E2 << '\n';
        throw MyException(0x571c8119U, "Error: sizes of bases do not match.");
    }

    /* Save the results */
    MyDB dbE2(db_out);
    dbE2.drop_and_create_generators_cell(table_out);
    dbE2.drop_and_create_relations(table_out);
    dbE2.drop_and_create_basis(table_out);

    dbE2.begin_transaction();

    dbE2.save_generators_cell(table_out, v_degs, toS0);
    dbE2.save_gb_mod(table_out, gbm);
    dbE2.save_basis_mod(table_out, basis, repr);

    dbE2.end_transaction();
}

int main_export(int argc, char** argv, int& index, const char* desc)
{
    std::string ring = "S0";
    int t_max = -1, stem_max = 383;

    myio::CmdArg1d args = {{"ring", &ring}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {{"stem_max", &stem_max}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    ExportRingAdamsE2(ring, t_max, stem_max);
    return 0;
}

int main_export_mod(int argc, char** argv, int& index, const char* desc)
{
    std::string mod, ring;
    int t_max = -1, stem_max = 383;

    myio::CmdArg1d args = {{"mod", &mod}, {"ring", &ring}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {{"stem_max", &stem_max}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;


    {
        auto js = myio::load_json("Adams.json");
        auto& cws = js.at("CW_complexes");
        if (cws.contains(mod)) {
            auto& cw_json = cws.at(mod);
            if (cw_json.contains("free")) {
                alg::int1d cells;
                for (auto& c : cw_json.at("cells_gen"))
                    cells.push_back(c.get<int>());
                ExportFreeModAdamsE2(mod, ring, cells, t_max, stem_max);
                return 0;
            }
            else if (cw_json.contains("summands")) {
				Export2Cell(mod, ring, t_max, stem_max);
				return 0;
            }
        }
    }
    
    ExportModAdamsE2(mod, ring, t_max, stem_max);
    return 0;
}

int main_export_d2(int argc, char** argv, int& index, const char* desc)
{
    std::string cw;

    myio::CmdArg1d args = {{"cw", &cw}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    ExportAdamsD2(cw);
    return 0;
}

int main_export_map(int argc, char** argv, int& index, const char* desc)
{
    std::string cw1, cw2;
    int t_max = -1;

    myio::CmdArg1d args = {{"cw1", &cw1}, {"cw2", &cw2}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    
    {
        auto js = myio::load_json("Adams.json");
        auto& maps = js.at("maps");
        auto name = fmt::format("{}__{}", cw1, cw2);
        if (maps.contains(name)) {
            auto& map_json = maps.at(name);
            if (map_json.contains("summands")) {
                auto summand0 = fmt::format("map_AdamsSS_{}.db", map_json.at("summands")[0].get<std::string>());
                auto summand1 = fmt::format("map_AdamsSS_{}.db", map_json.at("summands")[1].get<std::string>());
                ExportMapSumAdamsE2(cw1, cw2, summand0, summand1, 0, 1, myio::get(map_json, "sus", 0), myio::get(map_json, "fil", 0), t_max);
                return 0;
            }
            else if (map_json.contains("free_images")) {
                ExportMapFromFreeModAdamsE2(cw1, cw2, map_json, t_max);
				return 0;
            }
        }
    }
    ExportMapAdamsE2(cw1, cw2, t_max);
    return 0;
}

int main_2cell_export(int argc, char** argv, int& index, const char* desc)
{
    std::string mod;
    std::string ring = "S0";
    int t_max = -1, stem_max = 383;

    myio::CmdArg1d args = {{"mod:C2/Ceta/Cnu/Csigma", &mod}, {"ring", &ring}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {{"stem_max", &stem_max}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Export2Cell(mod, ring, t_max, stem_max);
    return 0;
}