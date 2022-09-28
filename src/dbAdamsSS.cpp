#include "dbAdamsSS.h"

namespace myio {

template <>
Mon Deserialize<Mon>(const std::string& str)
{
    Mon result;
    if (str.empty())
        return result;
    std::stringstream ss(str);
    while (ss.good()) {
        int gen, exp;
        ss >> gen >> "," >> exp;
        result.push_back(GE(gen, exp));
        if (ss.peek() == ',')
            ss.ignore();
    }
    return result;
}

template <>
Poly Deserialize<Poly>(const std::string& str)
{
    Poly result;
    if (str.empty())
        return result; /* Return 0 as a polynomial */
    else if (str == ";") {
        result.data.push_back({});
        return result; /* Return 1 as a polynomial */
    }
    std::istringstream ss(str);
    while (ss.good()) {
        int gen, exp;
        ss >> gen >> "," >> exp;
        if (result.data.empty())
            result.data.push_back(Mon());
        result.data.back().push_back(GE(gen, exp));
        if (ss.peek() == ',')
            ss.ignore();
        else if (ss.peek() == ';') {
            ss.ignore();
            result.data.push_back(Mon());
        }
    }
    return result;
}

template <>
MMod Deserialize<MMod>(const std::string& str)
{
    Mon m;
    int v;
    if (str.empty())
        throw MyException(0x2fb9914bU, "Cannot initialize MMod with an empty string.");
    std::stringstream ss(str);
    while (ss.good()) {
        int gen, exp;
        ss >> gen >> "," >> exp;
        if (ss.good()) {
            m.push_back(GE(gen, exp));
            if (ss.peek() == ',')
                ss.ignore();
        }
        else
            v = gen;
    }
    return MMod(m, v);
}

template <>
Mod Deserialize<Mod>(const std::string& str)
{
    Mod result;
    if (str.empty())
        return result; /* Return 0 as a polynomial */

    Mon m;
    int v = -1;
    std::istringstream ss(str);
    while (ss.good()) {
        int gen, exp;
        ss >> gen >> "," >> exp;
        if (ss.good()) {
            m.push_back(GE(gen, exp));
            if (ss.peek() == ',')
                ss.ignore();
        }
        else {
            v = gen;
            ss.clear();
            result.data.push_back(MMod(m, v));
        }
        if (ss.peek() == ',')
            ss.ignore();
        else if (ss.peek() == ';') {
            ss.ignore();
            m.clear();
            v = -1;
        }
    }
    return result;
}

void DbAdamsSS::save_generators(const std::string& table_prefix, const alg::AdamsDeg1d& gen_degs, alg::int1d& gen_repr) const
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

void DbAdamsSS::save_gb(const std::string& table_prefix, const std::map<AdamsDeg, Poly1d>& gb) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (rel, s, t) VALUES (?1, ?2, ?3);");
    int count = 0;
    for (auto& [deg, polys] : gb) {
        for (auto& poly : polys) {
            stmt.bind_str(1, Serialize(poly));
            stmt.bind_int(2, deg.s);
            stmt.bind_int(3, deg.t);
            stmt.step_and_reset();
            ++count;
        }
    }

    std::cout << count << " relations are inserted into " + table_prefix + "_relations!\n";
}

void DbAdamsSS::save_gb_mod(const std::string& table_prefix, const std::map<AdamsDeg, Mod1d>& gbm) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (rel, s, t) VALUES (?1, ?2, ?3);");

    int count = 0;
    for (auto& [deg, xs] : gbm) {
        for (auto& x : xs) {
            stmt.bind_str(1, Serialize(x));
            stmt.bind_int(2, deg.s);
            stmt.bind_int(3, deg.t);
            stmt.step_and_reset();
            ++count;
        }
    }

    std::cout << gbm.size() << " relations are inserted into " + table_prefix + "_relations!\n";
}

template <typename T>
AdamsDeg1d OrderDegsV2(const T& cont)
{
    AdamsDeg1d result;
    for (auto& [d, _] : cont) {
        result.push_back(d);
    }
    std::sort(result.begin(), result.end(), [](const AdamsDeg& d1, const AdamsDeg& d2) { return d1.t < d2.t || (d1.t == d2.t && d1.s > d2.s); });
    return result;
}

void DbAdamsSS::save_basis(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::Mon1d>& basis, const std::map<AdamsDeg, int2d>& repr) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (id, mon, repr, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.size(); ++i) {
            stmt.bind_int(1, count);
            stmt.bind_str(2, Serialize(basis_d[i]));
            stmt.bind_str(3, Serialize(repr.at(deg)[i]));
            stmt.bind_int(4, deg.s);
            stmt.bind_int(5, deg.t);
            stmt.step_and_reset();
            ++count;
        }
    }

    std::cout << count << " bases are inserted into " + table_prefix + "_basis!\n";
}

void DbAdamsSS::save_basis_mod(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::MMod1d>& basis, const std::map<AdamsDeg, int2d>& repr) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (id, mon, repr, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.size(); ++i) {
            stmt.bind_int(1, count);
            stmt.bind_str(2, Serialize(basis_d[i]));
            stmt.bind_str(3, Serialize(repr.at(deg)[i]));
            stmt.bind_int(4, deg.s);
            stmt.bind_int(5, deg.t);
            stmt.step_and_reset();
            ++count;
        }
    }

    std::cout << count << " bases are inserted into " + table_prefix + "_basis!\n";
}

AdamsDeg1d DbAdamsSS::load_gen_adamsdegs(const std::string& table_prefix) const
{
    AdamsDeg1d result;
    Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators ORDER BY id;");
    while (stmt.step() == MYSQLITE_ROW)
        result.push_back({stmt.column_int(0), stmt.column_int(1)});
    std::cout << "gen_adamsdegs loaded from " << table_prefix + "_generators, size=" << result.size() << '\n';
    return result;
}

Poly1d DbAdamsSS::load_gb(const std::string& table_prefix, int t_max) const
{
    Poly1d result;
    Statement stmt(*this, "SELECT rel FROM " + table_prefix + "_relations" + (t_max == alg::DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    while (stmt.step() == MYSQLITE_ROW) {
        Poly g = Deserialize<Poly>(stmt.column_str(0));
        result.push_back(std::move(g));
    }
    std::cout << "gb loaded from " << table_prefix + "_relations, size=" << result.size() << '\n';
    return result;
}

Mod1d DbAdamsSS::load_gb_mod(const std::string& table_prefix, int t_max) const
{
    Mod1d result;
    Statement stmt(*this, "SELECT rel FROM " + table_prefix + "_relations" + (t_max == alg::DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    while (stmt.step() == MYSQLITE_ROW) {
        Mod g = Deserialize<Mod>(stmt.column_str(0));
        result.push_back(std::move(g));
    }
    std::cout << "gb loaded from " << table_prefix + "_relations, size=" << result.size() << '\n';
    return result;
}

std::map<AdamsDeg, Mon1d> DbAdamsSS::load_basis(const std::string& table_prefix) const
{
    std::map<AdamsDeg, Mon1d> result;
    Statement stmt(*this, "SELECT s, t, mon FROM " + table_prefix + "_basis ORDER BY id");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
        result[deg].push_back(Deserialize<Mon>(stmt.column_str(2)));
    }
    std::cout << "basis loaded from " << table_prefix << "_basis, size=" << count << '\n';
    return result;
}

std::map<AdamsDeg, MMod1d> DbAdamsSS::load_basis_mod(const std::string& table_prefix) const
{
    std::map<AdamsDeg, MMod1d> result;
    Statement stmt(*this, "SELECT s, t, mon FROM " + table_prefix + "_basis ORDER BY id");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
        result[deg].push_back(Deserialize<MMod>(stmt.column_str(2)));
    }
    std::cout << "basis loaded from " << table_prefix << "_basis, size=" << count << '\n';
    return result;
}

}  // namespace myio