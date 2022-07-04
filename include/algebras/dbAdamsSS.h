#include "algebras.h"
#include "database.h"
#include "myio.h"
#include <map>

namespace myio {

using namespace alg;

inline std::string Serialize(const Mon& mon)
{
    return StrCont("", ",", "", "", mon, [](GenPow p) { return std::to_string(p.gen) + "," + std::to_string(p.exp); });
}

inline std::string Serialize(const Mon1d& obj)
{
    if (obj.size() == 1 && obj[0].empty())
        return ";";
    return StrCont("", ";", "", "", obj, [](const Mon& mon) { return Serialize(mon); });
}

template <>
Mon Deserialize<Mon>(const std::string& str);

template <>
Mon1d Deserialize<Mon1d>(const std::string& str);

/*****************************************************
 *             class DbAdamsSS
 *****************************************************/
class DbAdamsSS : public myio::Database
{
    using Statement = myio::Statement;

public:
    DbAdamsSS() = default;
    explicit DbAdamsSS(const std::string& filename) : Database(filename) {}

public:
    void create_generators(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, repr SMALLINT, s SMALLINT, t SMALLINT);");
    }
    void create_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (rel TEXT, s SMALLINT, t SMALLINT);");
    }
    void create_basis(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis (id INTEGER PRIMARY KEY, mon TEXT, s SMALLINT, t SMALLINT);");
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

public:
    void save_generators(const std::string& table_prefix, const alg::AdamsDeg1d& gen_degs, alg::int1d& gen_repr) const;
    void save_gb(const std::string& table_prefix, const alg::PolyRevlex1d& gb, const alg::AdamsDeg1d& gen_degs) const;
    void save_basis(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::Mon1d>& basis) const;

public:
    AdamsDeg1d load_gen_adamsdegs(const std::string& table_prefix) const;
    PolyRevlex1d load_gb(const std::string& table_prefix, int t_max) const;
    std::map<AdamsDeg, Mon1d> load_basis(const std::string& table_prefix) const;
};

}  // namespace myio