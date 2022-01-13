#include "dbalg.h"
#include <iostream>
#include <sqlite3.h>

namespace myio {

DbAlg::DbAlg(const char* filename) : Database(filename), bLogging_(true)
{
    try {
        version_ = get_str(std::string("SELECT value FROM info WHERE key=\"version\""));
        mo_ = get_str(std::string("SELECT value FROM info WHERE key=\"mo\""));
        date_ = get_str(std::string("SELECT value FROM info WHERE key=\"date\""));
    }
    catch (MyException&) {
    }
}

std::string DbAlg::Serialize(array::const_iterator pbegin, array::const_iterator pend)
{
    std::stringstream ss;
    for (auto p = pbegin; p < pend; p++) {
        ss << *p;
        if (p + 1 != pend)
            ss << ",";
    }
    return ss.str();
}

std::string DbAlg::Serialize(const Mon& obj)
{
    std::ostringstream ss;
    for (auto p = obj.begin(); p != obj.end(); ++p) {
        ss << p->gen << ',' << p->exp;
        if (p + 1 != obj.end())
            ss << ',';
    }
    return ss.str();
}

std::string DbAlg::Serialize(Mon1d::const_iterator pbegin, Mon1d::const_iterator pend) /* Warning: the algebra should be connected */
{
    std::ostringstream ss;
    if (pend - pbegin == 1 && pbegin->empty()) {
        ss << ';';
        return ss.str();
    }
    for (auto pMon = pbegin; pMon < pend; ++pMon) {
        for (auto p = pMon->begin(); p != pMon->end(); ++p) {
            ss << p->gen << ',' << p->exp;
            if (p + 1 != pMon->end())
                ss << ',';
        }
        if (pMon + 1 != pend)
            ss << ';';
    }
    return ss.str();
}

alg::array2d DbAlg::get_column_array(const std::string& table_name, const std::string& column_name, const std::string& conditions) const
{
    return get_column_from_str<array>(table_name, column_name, conditions, Deserialize<array>);
}

alg::MayDeg1d DbAlg::load_gen_maydegs(const std::string& table_prefix) const
{
    if (table_prefix.find("_generators") != std::string::npos)  // Deprecated
        throw MyException(0xccd51941U, "Should use prefix only.");
    std::vector<MayDeg> gen_degs;
    Statement stmt(*this, "SELECT s, t, v FROM " + table_prefix + "_generators ORDER BY gen_id;");
    while (stmt.step() == SQLITE_ROW)
        gen_degs.push_back({stmt.column_int(0), stmt.column_int(1), stmt.column_int(2)});
    std::cout << "gen_degs loaded from " << table_prefix + "_generators, size=" << gen_degs.size() << '\n';
    return gen_degs;
}

alg::Mon2d DbAlg::load_leading_terms(const std::string& table_prefix, int t_max) const
{
    if (table_prefix.find("_relations") != std::string::npos)  // Deprecated
        throw MyException(0x3b6c138U, "Should use prefix only.");
    Mon2d leadings;
    Statement stmt(*this, "SELECT leading_term FROM " + table_prefix + "_relations" + (t_max == -1 ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    int count = 0;
    while (stmt.step() == SQLITE_ROW) {
        ++count;
        Mon mon = Deserialize<Mon>(stmt.column_str(0));
        if (size_t(mon[0].gen) >= leadings.size())
            leadings.resize(size_t(mon[0].gen) + 1);
        leadings[mon[0].gen].push_back(mon);
    }
    std::cout << "leading_term loaded from " << table_prefix + "_relations, size=" << count << '\n';
    return leadings;
}

std::map<alg::MayDeg, int> DbAlg::load_indices(const std::string& table_prefix, int t_max) const
{
    if (table_prefix.find("_basis") != std::string::npos)  // Deprecated
        throw MyException(0x2b4b52b6U, "Should use prefix only.");
    std::map<MayDeg, int> result;
    Statement stmt(*this, "SELECT s, t, v, min(base_id) FROM " + table_prefix + "_basis" + (t_max == -1 ? "" : " WHERE t<=" + std::to_string(t_max)) + " GROUP BY s, t, v;");
    int count = 0;
    while (stmt.step() == SQLITE_ROW) {
        ++count;
        MayDeg d = {stmt.column_int(0), stmt.column_int(1), stmt.column_int(2)};
        result[d] = stmt.column_int(3);
    }
    std::cout << "indices loaded from " << table_prefix + "_basis, size=" << count << '\n';
    return result;
}

std::map<alg::MayDeg, alg::array2d> DbAlg::load_mon_diffs_ind_with_null(const std::string& table_prefix, int t_max) const
{
    if (table_prefix.find("_basis") != std::string::npos)  // Deprecated
        throw MyException(0x2eae2501U, "Should use prefix only.");
    using T = array;
    std::string column_name = "diff";
    std::string conditions = (t_max == -1 ? std::string("") : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY mon_id;";

    std::map<MayDeg, std::vector<T>> result;
    Statement stmt(*this, "SELECT s, t, v, " + column_name + " FROM " + table_prefix + "_basis " + conditions + ';');
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        MayDeg deg = {stmt.column_int(0), stmt.column_int(1), stmt.column_int(2)};
        if (stmt.column_type(3) == SQLITE_TEXT)
            result[deg].push_back(Deserialize<array>(stmt.column_str(3)));
        else
            result[deg].push_back({-1});
    }
    std::cout << column_name << " loaded from" << table_prefix + "_basis, size = " << count << '\n';
    return result;
}

std::map<alg::MayDeg, alg::BasisComplex> DbAlg::load_basis_ss(const std::string& table_prefix, int r, int t_max) const
{
    if (table_prefix.find("_ss") != std::string::npos)  // Deprecated
        throw MyException(0x2db87425U, "Should use prefix only.");
    std::map<MayDeg, alg::BasisComplex> basis_ss;
    Statement stmt(*this, "SELECT s, t, v, level, base FROM " + table_prefix + "_ss" + (t_max == -1 ? "" : " WHERE t<=" + std::to_string(t_max)) + " ;");
    int count = 0;
    while (stmt.step() == SQLITE_ROW) {
        ++count;
        MayDeg deg = {stmt.column_int(0), stmt.column_int(1), stmt.column_int(2)};
        int level = stmt.column_int(3);
        if (level <= r)
            basis_ss[deg].boundaries.push_back(Deserialize<array>(stmt.column_str(4)));
        else if (level <= alg::kLevelMax - r - 2)
            basis_ss[deg].cycles.push_back(Deserialize<array>(stmt.column_str(4)));
    }
    std::cout << "basis_ss loaded from " << table_prefix + "_ss, size=" << count << '\n';
    return basis_ss;
}

std::map<alg::MayDeg, alg::Staircase> DbAlg::load_basis_ss(const std::string& table_prefix, int t_max) const
{
    if (table_prefix.find("_ss") != std::string::npos)  // Deprecated
        throw MyException(0xfca647e1U, "Should use prefix only.");
    std::map<MayDeg, alg::Staircase> basis_ss;
    Statement stmt(*this, "SELECT s, t, v, level, base, diff FROM " + table_prefix + "_ss" + (t_max == -1 ? "" : " WHERE t<=" + std::to_string(t_max)) + " ;");
    int count = 0;
    while (stmt.step() == SQLITE_ROW) {
        ++count;
        MayDeg deg = {stmt.column_int(0), stmt.column_int(1), stmt.column_int(2)};
        int level = stmt.column_int(3);
        array base = Deserialize<array>(stmt.column_str(4));

        basis_ss[deg].basis_ind.push_back(std::move(base));
        basis_ss[deg].levels.push_back(level);
        if (stmt.column_type(5) == SQLITE_TEXT) {
            array diff = Deserialize<array>(stmt.column_str(5));
            basis_ss[deg].diffs_ind.push_back(diff);
        }
        else
            basis_ss[deg].diffs_ind.push_back({-1});
    }
    std::cout << "basis_ss loaded from " << table_prefix + "_ss, size=" << count << '\n';
    return basis_ss;
}

void DbAlg::save_ss(const std::string& table_prefix, const std::map<alg::MayDeg, alg::Staircase>& basis_ss) const
{
    if (table_prefix.find("_ss") != std::string::npos)  // Deprecated
        throw MyException(0xfbb58c88U, "Should use prefix only.");
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_ss (base, diff, level, s, t, v, base_id) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7);");

    int count = 0;
    for (const auto& [deg, basis_ss_d] : basis_ss) {
        for (size_t i = 0; i < basis_ss_d.basis_ind.size(); ++i) {
            stmt.bind_str(1, Serialize(basis_ss_d.basis_ind[i]));
            if (basis_ss_d.diffs_ind[i] == array{-1})
                stmt.bind_null(2);
            else
                stmt.bind_str(2, Serialize(basis_ss_d.diffs_ind[i]));
            stmt.bind_int(3, basis_ss_d.levels[i]);
            stmt.bind_int(4, deg.s);
            stmt.bind_int(5, deg.t);
            stmt.bind_int(6, deg.v);
            stmt.bind_int(7, count);
            stmt.step_and_reset();
            ++count;
        }
    }
    if (bLogging_)
        std::cout << "basis_ss is inserted into " + table_prefix + "_ss, size=" << count << '\n';
}

void DbAlg::update_ss(const std::string& table_prefix, const std::map<alg::MayDeg, alg::Staircase>& basis_ss) const
{
    if (table_prefix.find("_ss") != std::string::npos)  // Deprecated
        throw MyException(0x13aa4307U, "Should use prefix only.");
    std::map<MayDeg, int> indices = load_indices(table_prefix + "_ss", -1);
    Statement stmt(*this, "UPDATE " + table_prefix + "_ss SET base=?1, diff=?2, level=?3 WHERE base_id=?4;");

    int count = 0;
    for (const auto& [deg, basis_ss_d] : basis_ss) {
        for (size_t i = 0; i < basis_ss_d.basis_ind.size(); ++i) {
            stmt.bind_str(1, Serialize(basis_ss_d.basis_ind[i]));
            if (basis_ss_d.diffs_ind[i] == array{-1})
                stmt.bind_null(2);
            else
                stmt.bind_str(2, Serialize(basis_ss_d.diffs_ind[i]));
            stmt.bind_int(3, basis_ss_d.levels[i]);
            stmt.bind_int(4, indices[deg] + (int)i);
            stmt.step_and_reset();
            ++count;
        }
    }
    if (bLogging_)
        std::cout << "basis_ss " + table_prefix + "_ss is updated, num_of_change=" << count << '\n';
}

}  // namespace myio
