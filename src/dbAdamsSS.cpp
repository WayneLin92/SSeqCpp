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
        result.emplace_back(gen, exp);
        if (ss.peek() == ',')
            ss.ignore();
    }
    return result;
}

template <>
Mon1d Deserialize<Mon1d>(const std::string& str)
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

    std::clog << gen_degs.size() << " generators are inserted into " + table_prefix + "_generators!\n";
}

void DbAdamsSS::save_gb(const std::string& table_prefix, const alg::PolyRevlex1d& gb, const alg::AdamsDeg1d& gen_degs) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (rel, s, t) VALUES (?1, ?2, ?3);");

    for (size_t i = 0; i < gb.size(); ++i) {
        alg::AdamsDeg deg = alg::GetAdamsDeg(gb[i].GetLead(), gen_degs);
        stmt.bind_str(1, Serialize(gb[i].data));
        stmt.bind_int(2, deg.s);
        stmt.bind_int(3, deg.t);
        stmt.step_and_reset();
    }

    std::clog << gb.size() << " relations are inserted into " + table_prefix + "_relations!\n";
}

void DbAdamsSS::save_basis(const std::string& table_prefix, const std::map<alg::AdamsDeg, alg::Mon1d>& basis) const
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
}

AdamsDeg1d DbAdamsSS::load_gen_adamsdegs(const std::string& table_prefix) const
{
    AdamsDeg1d result;
    Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators ORDER BY id;");
    while (stmt.step() == MYSQLITE_ROW)
        result.push_back({stmt.column_int(0), stmt.column_int(1)});
    std::clog << "gen_adamsdegs loaded from " << table_prefix + "_generators, size=" << result.size() << '\n';
    return result;
}

PolyRevlex1d DbAdamsSS::load_gb(const std::string& table_prefix, int t_max) const
{
    PolyRevlex1d result;
    Statement stmt(*this, "SELECT rel FROM " + table_prefix + "_relations" + (t_max == alg::DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    while (stmt.step() == MYSQLITE_ROW) {
        PolyRevlex g = PolyRevlex::Sort(Deserialize<Mon1d>(stmt.column_str(0)));
        result.push_back(std::move(g));
    }
    std::clog << "gb loaded from " << table_prefix + "_relations, size=" << result.size() << '\n';
    return result;
}

std::map<AdamsDeg, Mon1d> DbAdamsSS::load_basis(const std::string& table_prefix) const
{
    std::map<AdamsDeg, Mon1d> result;
    Statement stmt(*this, "SELECT s, t, mon FROM " + table_prefix + "_basis");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
        result[deg].push_back(Deserialize<Mon>(stmt.column_str(2)));
    }
    std::clog << "basis loaded from " << table_prefix << "_basis, size=" << count << '\n';
    return result;
}

}