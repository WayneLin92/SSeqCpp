#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/groebner.h"
#include "algebras/linalg.h"
#include <cstring>
#include <map>

class DbSteenrod : public myio::Database
{
    using Statement = myio::Statement;

public:
    DbSteenrod() = default;
    explicit DbSteenrod(const std::string& filename) : Database(filename) {}

public:
    void create_generators(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, repr SMALLINT, s SMALLINT, t SMALLINT);");
    }
    void create_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (rel BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_basis(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis (id INTEGER PRIMARY KEY, mon BLOB, s SMALLINT, t SMALLINT);");
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
    void create_basis_and_delete(const std::string& table_prefix) const
    {
        create_basis(table_prefix);
        delete_from(table_prefix + "_basis");
    }

    void load_indecomposables(const std::string& table_prefix, alg::array& array_id, alg::AdamsDeg1d& gen_degs) const
    {
        Statement stmt(*this, "SELECT id, s, t FROM " + table_prefix + "_generators_products WHERE indecomposable=1 ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0), s = stmt.column_int(1), t = stmt.column_int(2);
            array_id.push_back(id);
            gen_degs.push_back(alg::AdamsDeg{s, t});
        }
    }

    /**
     * loc2glo[s][i] = id
     * glo2loc[id] = (s, i)
     */
    void load_id_converter(const std::string& table_prefix, alg::array2d& loc2glo, alg::pairint1d& glo2loc) const
    {
        Statement stmt(*this, "SELECT id, s FROM " + table_prefix + "_generators_products ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int id = stmt.column_int(0), s = stmt.column_int(1);
            if (loc2glo.size() <= (size_t)s)
                loc2glo.resize(size_t(s + 1));
            loc2glo[s].push_back(id);
            glo2loc.push_back(std::make_pair(s, (int)(loc2glo[s].size() - 1)));
        }
    }

    /* id_ind * id = prod */
    std::map<std::pair<int, int>, alg::array> load_products_h(const std::string& table_prefix) const
    {
        std::map<std::pair<int, int>, alg::array> result;
        alg::array id_inds = get_column_int(table_prefix + "_generators_products", "id", "WHERE indecomposable=1 ORDER BY id");
        for (int id_ind : id_inds) {
            std::string id_ind_str = std::to_string(id_ind);
            Statement stmt(*this, "SELECT id, mh" + id_ind_str + " FROM " + table_prefix + "_generators_products WHERE mh" + id_ind_str + " is not NULL ORDER BY id;");
            while (stmt.step() == MYSQLITE_ROW) {
                int id = stmt.column_int(0);
                alg::array prod_h = stmt.column_blob_tpl<int>(1);
                result[std::make_pair(id_ind, id)] = std::move(prod_h);
            }
        }
        return result;
    }

    void save_generators(const std::string& table_prefix, const alg::AdamsDeg1d& gen_degs, alg::array& gen_repr) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (id, repr, s, t) VALUES (?1, ?2, ?3, ?4);");

        for (size_t i = 0; i < gen_degs.size(); ++i) {
            stmt.bind_int(1, (int)i);
            stmt.bind_int(2, gen_repr[i]);
            stmt.bind_int(3, gen_degs[i].s);
            stmt.bind_int(4, gen_degs[i].t);
            stmt.step_and_reset();
        }

        std::cout << gen_degs.size() << " generators are inserted into " + table_prefix + "_generators!\n";
    }

    void save_gb(const std::string& table_prefix, const alg::PolyRevlex1d& gb, const alg::AdamsDeg1d& gen_degs) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (rel, s, t) VALUES (?1, ?2, ?3);");

        for (size_t i = 0; i < gb.size(); ++i) {
            alg::AdamsDeg deg = alg::GetAdamsDeg(gb[i].GetLead(), gen_degs);
            stmt.bind_blob(1, gb[i].data);
            stmt.bind_int(2, deg.s);
            stmt.bind_int(3, deg.t);
            stmt.step_and_reset();
        }

        std::cout << gb.size() << " relations are inserted into " + table_prefix + "_relations!\n";
    }

    void save_basis(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::Mon1d>& basis) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (mon, s, t) VALUES (?1, ?2, ?3);");

        int count = 0;
        for (auto& [deg, basis_d] : basis) {
            for (auto& m : basis_d) {
                ++count;
                stmt.bind_blob(1, m);
                stmt.bind_int(2, deg.s);
                stmt.bind_int(3, deg.t);
                stmt.step_and_reset();
            }
        }

        std::cout << count << " bases are inserted into " + table_prefix + "_basis!\n";
    }
};

alg::array mul(const std::map<std::pair<int, int>, alg::array>& map_h, int id_ind, const alg::array& repr)
{
    alg::array result;
    for (int id : repr) {
        auto it = map_h.find(std::make_pair(id_ind, id));
        if (it != map_h.end())
            result = lina::AddVectors(result, it->second);
    }
    return result;
}

void AdamsE2Export()
{
    using namespace alg;
    DbSteenrod dbProd("AdamsE2Prod.db");
    AdamsDeg1d gen_degs;
    array gen_reprs;
    dbProd.load_indecomposables("SteenrodMRes", gen_reprs, gen_degs);
    auto map_h_dual = dbProd.load_products_h("SteenrodMRes");

    array2d loc2glo;
    pairint1d glo2loc;
    dbProd.load_id_converter("SteenrodMRes", loc2glo, glo2loc);
    std::map<std::pair<int, int>, array> map_h;
    for (auto& [p, arr] : map_h_dual) {
        int s_i = glo2loc[p.second].first - glo2loc[p.first].first;
        for (int i : arr)
            map_h[std::make_pair(p.first, loc2glo[s_i][i])].push_back(p.second);
    }

    PolyRevlex1d gb;
    Mon2d leads;
    std::map<AdamsDeg, Mon1d> basis;
    std::map<AdamsDeg, array2d> repr;

    /* Add new basis */
    basis[AdamsDeg{0, 0}].push_back(Mon());
    repr[AdamsDeg{0, 0}].push_back({0});

    /* Consider generators one by one */
    int t_trunc = dbProd.get_int("SELECT MAX(t) from SteenrodMRes_generators_products");
    /* Add new basis */
    for (int t = 1; t <= t_trunc; t++) {
        std::map<AdamsDeg, Mon1d> basis_new;
        std::map<AdamsDeg, array2d> repr_new;
        std::cout << "t=" << t << '/' << t_trunc << '\n';

        /* Consider all possible basis in degree t */
        for (size_t gen_id = gen_degs.size(); gen_id-- > 0;) {
            auto repr_g = gen_reprs[gen_id];
            int t1 = t - gen_degs[gen_id].t;
            if (t1 >= 0) {
                auto p1 = basis.lower_bound(AdamsDeg{0, t1});
                auto p2 = basis.lower_bound(AdamsDeg{0, t1 + 1});
                for (auto p = p1; p != p2; ++p) {
                    for (size_t i = 0; i < p->second.size(); ++i) {
                        auto& m = p->second[i];
                        auto& repr_m = repr.at(p->first)[i];
                        if (m.empty() || (int)gen_id <= m[0].gen) {
                            Mon mon(mul(m, GenPow::Mon((int)gen_id)));
                            auto repr_mon = mul(map_h, repr_g, repr_m);
                            auto deg_mon = p->first + gen_degs[gen_id];
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
        for (auto it = basis_new.begin(); it != basis_new.end(); ++it)
        {
            auto indices = ut::size_t_range(it->second.size());
            std::sort(indices.begin(), indices.end(), [&it](size_t a, size_t b) { return CmpRevlex::template cmp<Mon, Mon>(it->second[a], it->second[b]); });
            Mon1d basis_new_d;
            array2d repr_new_d;
            for (size_t i = 0; i < indices.size(); ++i) {
                basis_new_d.push_back(it->second[indices[i]]);
                repr_new_d.push_back(repr_new.at(it->first)[indices[i]]);
            }

            array2d image, kernel, g;
            lina::SetLinearMap(repr_new_d, image, kernel, g);
            array lead_kernel;

            /* Add to groebner */
            for (const array& k : kernel) {
                lead_kernel.push_back(k[0]);
                gb.push_back(PolyRevlex::Sort(Indices2Poly(k, basis_new_d)));
                int index = gb.back().data.front().front().gen;
                if (leads.size() <= (size_t)index)
                    leads.resize(size_t(index + 1));
                leads[index].push_back(gb.back().data.front());
            }

            /* Add to basis_H */
            std::sort(lead_kernel.begin(), lead_kernel.end());
            array index_basis = lina::AddVectors(ut::int_range(int(basis_new_d.size())), lead_kernel);
            for (int i : index_basis) {
                basis[it->first].push_back(std::move(basis_new_d[i]));
                repr[it->first].push_back(std::move(repr_new_d[i]));
            }
        }
    }

    /* Save the results */
    DbSteenrod dbE2("AdamsE2Export.db");
    dbE2.create_generators_and_delete("AdamsE2");
    dbE2.create_relations_and_delete("AdamsE2");
    dbE2.create_basis_and_delete("AdamsE2");

    dbE2.begin_transaction();

    dbE2.save_generators("AdamsE2", gen_degs, gen_reprs);
    dbE2.save_gb("AdamsE2", gb, gen_degs);
    dbE2.save_basis("AdamsE2", basis);

    dbE2.end_transaction();
}

int main()
{
    AdamsE2Export();
    return 0;
}