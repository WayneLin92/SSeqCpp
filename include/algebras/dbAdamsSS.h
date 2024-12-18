#include "algebrasZ.h"
#include "database.h"
#include "myio.h"
#include <map>

namespace myio {

using namespace alg2;

/*****************************************************
 *             Serialization and Deserialization
 *****************************************************/
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

/*
 * Serializae algZ::Mon
 *
 * gen=0 means 2
 * gen=1 means fil
 * gen=2i means the i'th generator in even degree
 * gen=2i+1 means the i'th generator in odd degree
 */
std::string Serialize(const algZ::Mon& mon);

inline std::string Serialize(const algZ::MMod& mon)
{
    if (mon.IsUnKnown())
        return Serialize(mon.m) + ",-1";
    else
        return Serialize(mon.m) + "," + std::to_string(mon.v);
}

inline std::string Serialize(const Poly& poly)
{
    if (poly.data.size() == 1 && !poly.data[0])
        return ";";
    return StrCont("", ";", "", "", poly.data, [](const Mon& mon) { return Serialize(mon); });
}

inline std::string Serialize(const Mod& x)
{
    return StrCont("", ";", "", "", x.data, [](const MMod& mon) { return Serialize(mon); });
}

inline std::string Serialize(const algZ::Poly& poly)
{
    return StrCont("", ";", "", "", poly.data, [](const algZ::Mon& mon) { return Serialize(mon); });
}

inline std::string Serialize(const algZ::Mod& x)
{
    return StrCont("", ";", "", "", x.data, [](const algZ::MMod& mon) { return Serialize(mon); });
}

template <>
Mon Deserialize<Mon>(const std::string& str);

template <>
Poly Deserialize<Poly>(const std::string& str);

template <>
MMod Deserialize<MMod>(const std::string& str);

template <>
Mod Deserialize<Mod>(const std::string& str);

template <>
algZ::Mon Deserialize<algZ::Mon>(const std::string& str);

template <>
algZ::Poly Deserialize<algZ::Poly>(const std::string& str);

template <>
algZ::MMod Deserialize<algZ::MMod>(const std::string& str);

template <>
algZ::Mod Deserialize<algZ::Mod>(const std::string& str);



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
    void create_pi_generators(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, Einf TEXT, s SMALLINT, t SMALLINT);");
    }
    void create_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (rel TEXT, s SMALLINT, t SMALLINT);");
    }
    void create_pi_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_relations (rel TEXT, name TEXT, Einf TEXT, s SMALLINT, t SMALLINT);");
    }
    void create_basis(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis (id INTEGER PRIMARY KEY, mon TEXT, repr TEXT, s SMALLINT, t SMALLINT);");
    }
    void create_pi_basis(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_basis (id INTEGER PRIMARY KEY, mon TEXT, name TEXT, Einf TEXT, s SMALLINT, t SMALLINT);");
    }
    void drop_and_create_generators(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_generators");
        create_generators(table_prefix);
    }
    void drop_and_create_pi_generators(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_pi_generators");
        create_pi_generators(table_prefix);
    }
    void drop_and_create_relations(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_relations");
        create_relations(table_prefix);
    }
    void drop_and_create_pi_relations(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_pi_relations");
        create_pi_relations(table_prefix);
    }
    void drop_and_create_basis(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_basis");
        create_basis(table_prefix);
    }
    void drop_and_create_pi_basis(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_pi_basis");
        create_pi_basis(table_prefix);
    }

public:
    void save_generators(const std::string& table_prefix, const AdamsDeg1d& gen_degs, int1d& gen_repr) const;
    void save_gen_names(const std::string& table_prefix, const std::vector<std::string>& gen_names) const;
    void save_gb(const std::string& table_prefix, const std::map<AdamsDeg, Poly1d>& gb) const;
    void save_gb_mod(const std::string& table_prefix, const std::map<AdamsDeg, Mod1d>& gbm) const;
    void save_basis(const std::string& table_prefix, const std::map<AdamsDeg, Mon1d>& basis, const std::map<AdamsDeg, int2d>& repr) const;
    void save_basis_mod(const std::string& table_prefix, const std::map<AdamsDeg, MMod1d>& basis, const std::map<AdamsDeg, int2d>& repr) const;

public:
    AdamsDeg1d load_gen_adamsdegs(const std::string& table_prefix) const;
    std::vector<std::string> load_gen_names(const std::string& table_prefix) const;
    Poly1d load_gb(const std::string& table_prefix, int t_max) const;
    Mod1d load_gb_mod(const std::string& table_prefix, int t_max) const;
    BasisMon load_basis(const std::string& table_prefix) const;
    std::map<AdamsDeg, int2d> load_basis_d2(const std::string& table_prefix) const;
    BasisMMod load_basis_mod(const std::string& table_prefix) const;

public:
    void save_pi_generators(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Poly1d& gen_Einf) const;
    void save_pi_gb(const std::string& table_prefix, const std::map<AdamsDeg, algZ::Poly1d>& gb, const std::map<AdamsDeg, int2d>& gb_Einf) const;
    void save_pi_gb_mod(const std::string& table_prefix, const std::map<AdamsDeg, algZ::Mod1d>& gbm, const std::map<AdamsDeg, int2d>& gb_Einf) const;
    void save_pi_basis(const std::string& table_prefix, const std::map<AdamsDeg, algZ::Mon1d>& basis, const std::map<AdamsDeg, int2d>& basis_Einf) const;
    void save_pi_basis_mod(const std::string& table_prefix, const std::map<AdamsDeg, algZ::MMod1d>& basis, const std::map<AdamsDeg, int2d>& basis_Einf) const;

public:
    AdamsDeg1d load_pi_gen_adamsdegs(const std::string& table_prefix) const
    {
        return load_gen_adamsdegs(table_prefix + "_pi");
    }
    algZ::Poly1d load_pi_gb(const std::string& table_prefix, int t_max) const;
    algZ::Mod1d load_pi_gb_mod(const std::string& table_prefix, int t_max) const;
    std::map<AdamsDeg, algZ::Mon1d> load_pi_basis(const std::string& table_prefix) const;
    std::map<AdamsDeg, algZ::MMod1d> load_pi_basis_mod(const std::string& table_prefix) const;
};

}  // namespace myio
