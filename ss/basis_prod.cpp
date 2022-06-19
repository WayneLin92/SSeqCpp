#include "algebras/dbAdamsSS.h"
#include "algebras/groebner.h"
#include "algebras/myhash.h"

using namespace alg;

class MyDB : public myio::DbAdamsSS
{
    using Statement = myio::Statement;

public:
    MyDB() = default;
    explicit MyDB(const std::string& filename) : DbAdamsSS(filename) {}

    void create_basis_products(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis_products (id1 INTEGER, id2 INTEGER, prod TEXT, PRIMARY KEY(id1, id2));");
    }

    void create_basis_products_and_delete(const std::string& table_prefix) const
    {
        create_basis_products(table_prefix);
        delete_from(table_prefix + "_basis_products");
    }

    void load_basis_v2(const std::string& table_prefix, int t_max, Mon1d& basis, array& t_basis) const
    {
        Statement stmt(*this, "SELECT mon, t FROM " + table_prefix + "_basis" + (t_max == alg::DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            basis.push_back(myio::Deserialize<Mon>(stmt.column_str(0)));
            t_basis.push_back(stmt.column_int(1));
        }
        std::clog << "basis loaded from " << table_prefix << "_basis, size=" << basis.size() << '\n';
    }

    /* void save_basis_products(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::Mon1d>& basis) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (id, mon, s, t) VALUES (?1, ?2, ?3, ?4);");

        int count = 0;
        for (auto& [deg, basis_d] : basis) {
            for (auto& m : basis_d) {
                stmt.bind_int(1, count);
                stmt.bind_str(2, Serialize(m));
                stmt.bind_int(3, deg.s);
                stmt.bind_int(4, deg.t);
                stmt.step_and_reset();
                ++count;
            }
        }

        std::clog << count << " bases are inserted into " + table_prefix + "_basis!\n";
    }*/
};

namespace ut {
template <>
uint64_t hash<GenPow>(const GenPow& p)
{
    return (uint64_t(p.gen) << 32) + uint64_t(p.exp);
}
}  // namespace ut

array PolyHash2indices(const PolyRevlex& poly, std::map<uint64_t, int>& hash2index)
{
    array result;
    for (const Mon& m : poly.data) {
        result.push_back(hash2index.at(ut::hash(m)));
    }
    return result;
}

std::ostream& operator<<(std::ostream& sout, const Mon& mon)
{
    return sout << myio::Serialize(mon);
}

int main()
{
    using namespace alg;
    MyDB db("AdamsE2Export.db");

    int t_max = 184;

    AdamsDeg1d gen_degs = db.load_gen_adamsdegs("AdamsE2");
    PolyRevlex1d polys = db.load_gb("AdamsE2", t_max);
    Mon1d basis;
    array t_basis;
    db.load_basis_v2("AdamsE2", t_max, basis, t_basis);

    GroebnerRevlex gb(184, polys);

    std::map<uint64_t, int> hash2index;
    for (size_t i = 0; i < basis.size(); ++i) {
        auto h = ut::hash(basis[i]);
        if (hash2index.find(h) != hash2index.end()) {
            std::cout << h << '\n';
            std::cout << basis[hash2index.at(h)] << '\n';
            std::cout << basis[i] << '\n';
            throw MyException(0xf0c20f8aU, "Hash collision!");
        }
        hash2index[ut::hash(basis[i])] = (int)i;
    }

    db.create_basis_products_and_delete("AdamsE2");
    db.begin_transaction();
    myio::Statement stmt(db, "INSERT INTO AdamsE2_basis_products (id1, id2, prod) VALUES (?1, ?2, ?3);");
    for (size_t i = 0; i < basis.size(); ++i) {
        std::cout << "i=" << i << '\n';
        for (size_t j = i; j < basis.size(); ++j) {
            if (t_basis[i] + t_basis[j] <= t_max) {
                PolyRevlex poly_prod = gb.Reduce(PolyRevlex::Mon_(basis[i]) * PolyRevlex::Mon_(basis[j]));

                //std::cout << "(" << basis[i] << ") * (" << basis[j] << ") = " << myio::Serialize(poly_prod.data) << '\n';
                array prod = PolyHash2indices(poly_prod, hash2index);
                stmt.bind_int(1, (int)i);
                stmt.bind_int(2, (int)j);
                stmt.bind_str(3, myio::Serialize(prod));
                stmt.step_and_reset();
            }
        }
    }
    db.end_transaction();
    return 0;
}