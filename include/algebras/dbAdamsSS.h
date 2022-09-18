#include "algebras.h"
#include "database.h"
#include "myio.h"
#include <map>

namespace myio {

using namespace alg;

inline std::string Serialize(const Mon& mon)
{
    return StrCont("", ",", "", "", mon, [](GE p) { return std::to_string(p.g()) + "," + std::to_string(p.e()); });
}

inline std::string Serialize(const MMod& mon)
{
    if (mon.m)
        return Serialize(mon.m) + "," + std::to_string(mon.v);
    else
        return std::to_string(mon.v);
}

inline std::string Serialize(const Mon1d& obj)
{
    if (obj.size() == 1 && !obj[0])
        return ";";
    return StrCont("", ";", "", "", obj, [](const Mon& mon) { return Serialize(mon); });
}

inline std::string Serialize(const MMod1d& obj)
{
    return StrCont("", ";", "", "", obj, [](const MMod& mon) { return Serialize(mon); });
}

template <>
Mon Deserialize<Mon>(const std::string& str);

template <>
Mon1d Deserialize<Mon1d>(const std::string& str);

template <>
MMod Deserialize<MMod>(const std::string& str);

template <>
MMod1d Deserialize<MMod1d>(const std::string& str);

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
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis (id INTEGER PRIMARY KEY, mon TEXT, repr TEXT, s SMALLINT, t SMALLINT);");
    }
    void drop_and_create_generators(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_generators");
        create_generators(table_prefix);
    }
    void drop_and_create_relations(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_relations");
        create_relations(table_prefix);
    }
    void drop_and_create_basis(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_basis");
        create_basis(table_prefix);
    }

public:
    void save_generators(const std::string& table_prefix, const alg::AdamsDeg1d& gen_degs, alg::int1d& gen_repr) const;
    void save_gb(const std::string& table_prefix, const std::map<AdamsDeg, Poly1d>& gb) const;
    void save_gb_mod(const std::string& table_prefix, const std::map<AdamsDeg, Mod1d>& gbm) const;
    void save_basis(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::Mon1d>& basis, const std::map<AdamsDeg, int2d>& repr) const;
    void save_basis_mod(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::MMod1d>& basis, const std::map<AdamsDeg, int2d>& repr) const;

public:
    AdamsDeg1d load_gen_adamsdegs(const std::string& table_prefix) const;
    Poly1d load_gb(const std::string& table_prefix, int t_max) const;
    Mod1d load_gb_mod(const std::string& table_prefix, int t_max) const;
    std::map<AdamsDeg, Mon1d> load_basis(const std::string& table_prefix) const;
    std::map<AdamsDeg, MMod1d> load_basis_mod(const std::string& table_prefix) const;
};

}  // namespace myio