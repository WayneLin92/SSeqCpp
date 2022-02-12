/** \file algdb.h
 * A module for interacting with a database.
 * TODO: treat FnCmp properly
 */
#ifndef ALGDB_H_INCLUDED
#define ALGDB_H_INCLUDED

#include "database.h"
#include "groebner.h"
#include "myio.h"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <map>
#include <sstream>

namespace alg {

constexpr auto kLevelMax = 10000;

struct BasisComplex
{
    alg::array2d boundaries;
    alg::array2d cycles;
};

struct Staircase
{
    alg::array2d basis_ind;
    alg::array2d diffs_ind; /* ={-1} means null */
    alg::array levels;
};

}  // namespace alg

namespace myio {

class DbAlg : public Database
{
    using MayDeg = alg::MayDeg;
    using MayDeg1d = alg::MayDeg1d;
    using array = alg::array;
    using array2d = alg::array2d;
    using Mon = alg::Mon;
    using Mon1d = alg::Mon1d;
    using Mon2d = alg::Mon2d;
    template <typename FnCmp>
    using Polynomial = alg::Polynomial<FnCmp>;
    template <typename FnCmp>
    using Poly1d = std::vector<Polynomial<FnCmp>>;
    template <typename FnCmp>
    using Groebner = alg::Groebner<FnCmp>;

public:
    bool bLogging_ = true;

private:
    const char* VERSION = "0.0"; /** Version of this class */
    std::string version_;        /** Version of the class who created the opened database */
    std::string mo_;             /** Monomial ordering */
    std::string date_;           /** Latest modification date */

public:
    DbAlg() = default;
    explicit DbAlg(const std::string& filename);

public:
    void create_generators(const std::string& table_prefix)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (gen_id INTEGER PRIMARY KEY, gen_name TEXT UNIQUE, gen_diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
    }
    void create_relations(const std::string& table_prefix)
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
    }
    void create_generators_and_delete(const std::string& table_prefix)
    {
        create_generators(table_prefix);
        delete_from(table_prefix + "_generators");
    }
    void create_relations_and_delete(const std::string& table_prefix)
    {
        create_relations(table_prefix);
        delete_from(table_prefix + "_relations");
    }

public:
    template <typename T>
    static T Deserialize(const std::string& str)
    {
        static_assert(false, "must use a specialized template.");
    }

    template <>
    static array Deserialize<array>(const std::string& str)
    {
        array result;
        if (str.empty())
            return result;
        std::stringstream ss(str);
        while (ss.good()) {
            int i;
            ss >> i;
            result.push_back(i);
            if (ss.peek() == ',')
                ss.ignore();
        }
        return result;
    }

    template <>
    static Mon Deserialize<Mon>(const std::string& str)
    {
        Mon result;
        if (str.empty())
            return result;
        std::stringstream ss(str);
        while (ss.good()) {
            int gen, exp;
            ss >> gen >> "," >> exp;
            result.emplace_back(gen, exp);
            if (ss.peek() == ',')
                ss.ignore();
        }
        return result;
    }

    template <>
    static Mon1d Deserialize<Mon1d>(const std::string& str)
    {
        Mon1d result;
        if (str.empty())
            return result; /* Return 0 as a polynomial */
        else if (str == ";") {
            result.push_back({});
            return result; /* Return 1 as a polynomial */
        }
        std::istringstream ss(str);
        while (ss.good()) {
            int gen, exp;
            ss >> gen >> "," >> exp;
            if (result.empty())
                result.emplace_back();
            result.back().emplace_back(gen, exp);
            if (ss.peek() == ',')
                ss.ignore();
            else if (ss.peek() == ';') {
                ss.ignore();
                result.emplace_back();
            }
        }
        return result;
    }

public:
    static std::string Serialize(array::const_iterator pbegin, array::const_iterator pend);
    static std::string Serialize(const array& obj)
    {
        return Serialize(obj.begin(), obj.end());
    }
    static std::string Serialize(const Mon& obj);
    static std::string Serialize(Mon1d::const_iterator pbegin, Mon1d::const_iterator pend);
    static std::string Serialize(const Mon1d& obj)
    {
        return Serialize(obj.begin(), obj.end());
    }

public:
    std::vector<array> get_column_array(const std::string& table_name, const std::string& column_name, const std::string& conditions) const;
    template <typename FnCmp>
    Poly1d<FnCmp> get_column_Poly(const std::string& table_name, const std::string& column_name, const std::string& conditions) const
    {
        return get_column_from_str<Polynomial<FnCmp>>(table_name, column_name, conditions, [](const std::string& s) { return Polynomial<FnCmp>{Deserialize<Mon1d>(s)}; });
    }
    template <typename FnCmp>
    Poly1d<FnCmp> get_column_Poly_with_null(const std::string& table_name, const std::string& column_name, const Polynomial<FnCmp>& null_poly, const std::string& conditions) const
    {
        return get_column_from_str_with_null<Polynomial<FnCmp>>(table_name, column_name, null_poly, conditions, [](const std::string& s) { return Polynomial<FnCmp>{Deserialize<Mon1d>(s)}; });
    }

public:
    MayDeg1d load_gen_maydegs(const std::string& table_prefix) const;
    std::vector<std::string> load_gen_names(const std::string& table_prefix) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0x1d983ee9U, "Should use prefix only.");
        return get_column_str(table_prefix + "_generators", "gen_name", "ORDER BY gen_id");
    }
    template <typename FnCmp>
    Poly1d<FnCmp> load_gen_diffs(const std::string& table_prefix) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0x323a1a78U, "Should use prefix only.");
        return get_column_Poly_with_null<FnCmp>(table_prefix + "_generators", "gen_diff", Polynomial<FnCmp>::Gen(-1), " ORDER BY gen_id");
    }
    template <typename FnCmp>
    Poly1d<FnCmp> load_gen_reprs(const std::string& table_prefix) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0x4c3e54bcU, "Should use prefix only.");
        return get_column_Poly<FnCmp>(table_prefix + "_generators", "repr", "ORDER BY gen_id");
    }
    template <typename FnCmp>
    Poly1d<FnCmp> load_gen_images(const std::string& table_prefix, const std::string& column_name) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0xc18ab9ceU, "Should use prefix only.");
        return get_column_Poly_with_null<FnCmp>(table_prefix + "_generators", column_name, Polynomial<FnCmp>::Gen(-1), " ORDER BY gen_id");
    }
    Mon2d load_leading_terms(const std::string& table_prefix, int t_max) const; /* The leading monomials are grouped by the first generator */
    template <typename FnCmp>
    Groebner<FnCmp> load_gb(const std::string& table_prefix, int t_max) const
    {
        if (table_prefix.find("_relations") != std::string::npos)  // Deprecated
            throw MyException(0xb4a5e736U, "Should use prefix only.");
        Groebner<FnCmp> gb(t_max);
        Statement stmt(*this, "SELECT leading_term, basis FROM " + table_prefix + "_relations" + (t_max == alg::DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
        while (stmt.step() == MYSQLITE_ROW) {
            Polynomial<FnCmp> lead = {{Deserialize<Mon>(stmt.column_str(0))}};
            Polynomial<FnCmp> basis = Polynomial<FnCmp>::Sort(Deserialize<Mon1d>(stmt.column_str(1)));  ////

            Polynomial<FnCmp> g = basis + lead;
            gb.push_back(std::move(g));
        }
        std::cout << "gb loaded from " << table_prefix + "_relations, size=" << gb.size() << '\n';
        return gb;
    }
    std::map<MayDeg, int> load_indices(const std::string& table_prefix, int t_max) const;
    template <typename T, typename FnMap>
    std::map<MayDeg, std::vector<T>> group_by_Maydeg(const std::string& table_name, const std::string& column_name, const std::string& conditions, FnMap map) const
    {
        std::map<MayDeg, std::vector<T>> result;
        Statement stmt(*this, "SELECT s, t, v, " + column_name + " FROM " + table_name + ' ' + conditions + ';');
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            ++count;
            MayDeg deg = {stmt.column_int(0), stmt.column_int(1), stmt.column_int(2)};
            result[deg].push_back(map(stmt.column_str(3)));
        }
        std::cout << column_name << "'s loaded from" << table_name << ", size = " << count << '\n';
        return result;
    }
    std::map<MayDeg, Mon1d> load_basis(const std::string& table_prefix, int t_max) const
    {
        if (table_prefix.find("_basis") != std::string::npos)  // Deprecated
            throw MyException(0x7c907e53U, "Should use prefix only.");
        return group_by_Maydeg<Mon>(table_prefix + "_basis", "mon", (t_max == alg::DEG_MAX ? std::string() : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY mon_id;", Deserialize<Mon>);
    }
    std::map<MayDeg, array2d> load_mon_diffs_ind(const std::string& table_prefix, int t_max) const
    {
        if (table_prefix.find("_basis") != std::string::npos)  // Deprecated
            throw MyException(0xb5f0285cU, "Should use prefix only.");
        return group_by_Maydeg<array>(table_prefix + "_basis", "diff", (t_max == alg::DEG_MAX ? " WHERE" : " WHERE t<=" + std::to_string(t_max) + " AND") + " diff IS NOT NULL ORDER BY mon_id;", Deserialize<array>);
    }
    std::map<MayDeg, array2d> load_mon_diffs_ind_with_null(const std::string& table_prefix, int t_max) const;
    template <typename FnCmp>
    std::map<MayDeg, Poly1d<FnCmp>> load_mon_diffs(const std::string& table_prefix, const std::map<MayDeg, Mon1d>& basis, int r, int t_max) const
    {
        if (table_prefix.find("_basis") != std::string::npos)  // Deprecated
            throw MyException(0x6e87250aU, "Should use prefix only.");
        std::map<MayDeg, Poly1d<FnCmp>> result;
        std::string sql = "SELECT s, t, v, diff FROM " + table_prefix + "_basis" + (t_max == alg::DEG_MAX ? " WHERE" : " WHERE t<=" + std::to_string(t_max) + " AND") + " diff IS NOT NULL ORDER BY mon_id;";
        Statement stmt(*this, sql);
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            ++count;
            MayDeg deg = {stmt.column_int(0), stmt.column_int(1), stmt.column_int(2)};
            array diff_index = Deserialize<array>(stmt.column_str(3));
            result[deg].push_back(diff_index.empty() ? Polynomial<FnCmp>{} : Polynomial<FnCmp>{alg::Indices2Poly(std::move(diff_index), basis.at(deg + MayDeg{1, 0, -r}))});
        }
        std::cout << "diffs loaded from " << table_prefix + "_basis, size=" << count << '\n';
        return result;
    }
    std::map<MayDeg, alg::BasisComplex> load_basis_ss(const std::string& table_prefix, int r, int t_max) const; /* load d_r-cycles and d_r-boundaries */
    std::map<MayDeg, alg::Staircase> load_basis_ss(const std::string& table_prefix, int t_max) const;

public:
    template <typename T, typename FnMap>
    void save_column(const std::string& table_name, const std::string& column_name, const std::string& index_name, const std::vector<T>& column, FnMap map, size_t i_start) const
    {
        Statement stmt(*this, "UPDATE " + table_name + " SET " + column_name + " = ?1 WHERE " + index_name + "= ?2;");
        for (size_t i = i_start; i < column.size(); ++i) {
            stmt.bind_str(1, map(column[i]));
            stmt.bind_int(2, (int)i);
            stmt.step_and_reset();
        }
        if (bLogging_)
            std::cout << column.size() - i_start << ' ' << column_name << "'s are inserted into " + table_name + "!\n";
    }
    void save_gen_maydegs(const std::string& table_prefix, const std::vector<MayDeg>& gen_degs) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0xc913403eU, "Should use prefix only.");

        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (s, t, v, gen_id) VALUES (?1, ?2, ?3, ?4);");

        for (size_t i = 0; i < gen_degs.size(); ++i) {
            stmt.bind_int(1, gen_degs[i].s);
            stmt.bind_int(2, gen_degs[i].t);
            stmt.bind_int(3, gen_degs[i].v);
            stmt.bind_int(4, (int)i);
            stmt.step_and_reset();
        }
        if (bLogging_)
            std::cout << gen_degs.size() << " generators are inserted into " + table_prefix + "_generators!\n";
    }
    void save_gen_names(const std::string& table_prefix, const std::vector<std::string>& gen_names) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0x8df979a9U, "Should use prefix only.");
        save_column(
            table_prefix + "_generators", "gen_name", "gen_id", gen_names, [](const std::string& c) { return c; }, 0);
    }
    template <typename FnCmp>
    void save_gen_reprs(const std::string& table_prefix, const std::vector<Polynomial<FnCmp>>& gen_reprs) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0x8df979a9U, "Should use prefix only.");
        save_column(
            table_prefix + "_generators", "repr", "gen_id", gen_reprs, [](const Polynomial<FnCmp>& p) { return Serialize(p.data); }, 0);
    }
    template <typename FnCmp>
    void save_gen_diffs(const std::string& table_prefix, const std::vector<Polynomial<FnCmp>>& gen_diffs) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0x8df979a9U, "Should use prefix only.");
        save_column(
            table_prefix + "_generators", "gen_diff", "gen_id", gen_diffs, [](const Polynomial<FnCmp>& p) { return Serialize(p.data); }, 0);
    }
    template <typename FnCmp>
    void save_gen_images(const std::string& table_prefix, const std::string& column_name, const std::vector<Polynomial<FnCmp>>& gen_images) const
    {
        if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
            throw MyException(0x8df979a9U, "Should use prefix only.");
        save_column(
            table_prefix + "_generators", column_name, "gen_id", gen_images, [](const Polynomial<FnCmp>& p) { return Serialize(p.data); }, 0);
    }
    template <typename FnCmp>
    void save_gb(const std::string& table_prefix, const Groebner<FnCmp>& gb, const std::vector<MayDeg>& gen_degs, size_t i_start = 0) const
    {
        if (table_prefix.find("_relations") != std::string::npos)  // Deprecated
            throw MyException(0x8df979a9U, "Should use prefix only.");

        Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (leading_term, basis, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5);");

        for (size_t i = i_start; i < gb.size(); ++i) {
            MayDeg deg = gb[i].GetMayDeg(gen_degs);
            stmt.bind_str(1, Serialize(gb[i].GetLead()));
            stmt.bind_str(2, Serialize(gb[i].data.begin() + 1, gb[i].data.end()));
            stmt.bind_int(3, deg.s);
            stmt.bind_int(4, deg.t);
            stmt.bind_int(5, deg.v);
            stmt.step_and_reset();
        }

        execute_cmd("CREATE TABLE IF NOT EXISTS info (key TEXT PRIMARY KEY, value TEXT);");
        execute_cmd("DELETE FROM info;");

        Statement stmt1(*this, "INSERT INTO info (key, value) VALUES (?1, ?2);");
        stmt1.bind_str(1, "version");
        stmt1.bind_str(2, VERSION);
        stmt1.step_and_reset();

        stmt1.bind_str(1, "mo");
        stmt1.bind_str(2, std::string(FnCmp::name));
        stmt1.step_and_reset();

        stmt1.bind_str(1, "date");
        stmt1.bind_str(2, ut::get_time());
        stmt1.step_and_reset();

        stmt1.bind_str(1, "grading");
        stmt1.bind_str(2, "t,s,v");
        stmt1.step_and_reset();

        if (bLogging_)
            std::cout << gb.size() - i_start << " relations are inserted into " + table_prefix + "_relations!\n";
    }

    /*
     * Save the values of a std::map
     */
    template <typename T, typename FnMap>
    void save_groups(const std::string& table_name, const std::string& column_name, const std::string& index_name, const std::map<MayDeg, std::vector<T>>& groups, FnMap map) const
    {
        Statement stmt(*this, "UPDATE " + table_name + " SET " + column_name + " = ?1 WHERE " + index_name + "= ?2;");

        int count = 0;
        for (auto& [deg, v] : groups) {
            for (auto& m : v) {
                ++count;
                stmt.bind_str(1, map(m));
                stmt.bind_int(2, count);
                stmt.step_and_reset();
            }
        }
        if (bLogging_)
            std::cout << count << ' ' << column_name << "'s are inserted into " + table_name + "!\n";
    }
    void save_basis(const std::string& table_prefix, const std::map<MayDeg, Mon1d>& basis) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_basis (mon_id INTEGER PRIMARY KEY, mon TEXT NOT NULL, diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (mon, s, t, v) VALUES (?1, ?2, ?3, ?4);");

        int count = 0;
        for (auto& [deg, basis_d] : basis) {
            for (auto& m : basis_d) {
                ++count;
                stmt.bind_str(1, Serialize(m));
                stmt.bind_int(2, deg.s);
                stmt.bind_int(3, deg.t);
                stmt.bind_int(4, deg.v);
                stmt.step_and_reset();
            }
        }
        if (bLogging_)
            std::cout << count << " bases are inserted into " + table_prefix + "_basis!\n";
    }
    void save_mon_reprs(const std::string& table_prefix, const std::map<MayDeg, array2d>& mon_reprs) const
    {
        if (table_prefix.find("_basis") != std::string::npos)  // Deprecated
            throw MyException(0x8df979a9U, "Should use prefix only.");
        save_groups(table_prefix + "_basis", "repr", "mon_id", mon_reprs, [](const alg::array& a) { return Serialize(a); });
    }
    template <typename FnCmp>
    void save_mon_reprs(const std::string& table_prefix, const std::map<MayDeg, Poly1d<FnCmp>>& mon_reprs) const
    {
        if (table_prefix.find("_basis") != std::string::npos)  // Deprecated
            throw MyException(0x8df979a9U, "Should use prefix only.");
        save_groups(table_prefix + "_basis", "repr", "mon_id", mon_reprs, [](const Polynomial<FnCmp>& a) { return Serialize(a); });
    }
    void save_ss(const std::string& table_prefix, const std::map<MayDeg, alg::Staircase>& basis_ss) const;
    void update_ss(const std::string& table_prefix, const std::map<MayDeg, alg::Staircase>& basis_ss) const;
};

#if defined(_MSC_VER)
constexpr std::string_view dir_db = "C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/";
#else
constexpr std::string_view dir_db = "/Users/weinanlin/MyData/Math_AlgTop/databases/";
#endif

}  // namespace myio

#endif /* ALGDB_H_INCLUDED */
