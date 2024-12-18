#include "dbAdamsSS.h"

namespace myio {

std::string Serialize(const algZ::Mon& mon)
{
    static const uint32_t GEN_FIL = 0xffff;
    std::string result;
    if (mon.IsUnKnown())
        result = "-1," + std::to_string(mon.fil());
    else {
        if (mon.c() > 0) {
            result = "0," + std::to_string(mon.c());
        }
        if (mon.m()) {
            if (!result.empty())
                result += ',';
            Mon m0;
            for (auto& p : mon.m())
                m0.push_back(GE(p.g() * 2 + (p.e() & algZ::FLAG_ODD_EXP ? 1 : 0), p.e_masked()));
            result += Serialize(m0);
        }
        if (!result.empty())
            result += ',';
        result += "1," + std::to_string(mon.fil());
    }
    return result;
}

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
    int v = -1;
    if (str.empty())
        throw ErrorIdMsg(0x2fb9914bU, "Cannot initialize MMod with an empty string.");
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

/* Convert auxiliary input Mon to algZ::Mon */
algZ::Mon MonToMonZ(const Mon& m)
{
    algZ::Mon result;
    if (m.size() == 1 && m.begin()->g() == uint16_t(-1))
        result = algZ::Mon::O(m.begin()->e());
    else {
        int c = 0, fil = -1;
        Mon mm;
        for (GE p : m) {
            if (p.g() == 0)
                c = p.e();
            else if (p.g() == 1)
                fil = p.e();
            else if (p.g() % 2 == 0)
                mm = mm * GE(p.g() / 2, p.e());
            else
                mm = mm * GE(p.g() / 2, p.e() | algZ::FLAG_ODD_EXP);
        }
        if (fil == -1) {
            std::cout << "Incorrect input for algZ::Mon: " << m << '\n';
            throw ErrorIdMsg(0x3944e4aeU, "Incorrect input for algZ::Mon");
        }
        result = algZ::Mon(c, mm, fil);
    }
    return result;
}

/* Convert auxiliary input Mon to algZ::Mon */
algZ::Poly PolyToPolyZ(const Poly& p)
{
    algZ::Poly result;
    for (auto& m : p.data)
        result.data.push_back(MonToMonZ(m));
    return result;
}

/* Convert auxiliary input Mon to algZ::Mon */
algZ::MMod MModToMModZ(const MMod& p)
{
    return algZ::MMod(MonToMonZ(p.m), p.v, 0);
}

/* Convert auxiliary input Mon to algZ::Mon */
algZ::Mod ModToModZ(const Mod& p)
{
    algZ::Mod result;
    for (auto& m : p.data)
        result.data.push_back(MModToMModZ(m));
    return result;
}

template <>
algZ::Mon Deserialize<algZ::Mon>(const std::string& str)
{
    Mon m = Deserialize<Mon>(str);
    return MonToMonZ(m);
}

template <>
algZ::Poly Deserialize<algZ::Poly>(const std::string& str)
{
    Poly p = Deserialize<Poly>(str);
    return PolyToPolyZ(p);
}

template <>
algZ::MMod Deserialize<algZ::MMod>(const std::string& str)
{
    MMod m = Deserialize<MMod>(str);
    return MModToMModZ(m);
}

template <>
algZ::Mod Deserialize<algZ::Mod>(const std::string& str)
{
    Mod p = Deserialize<Mod>(str);
    return ModToModZ(p);
}

void DbAdamsSS::save_generators(const std::string& table_prefix, const alg2::AdamsDeg1d& gen_degs, alg2::int1d& gen_repr) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (id, repr, s, t) VALUES (?1, ?2, ?3, ?4);");

    for (size_t i = 0; i < gen_degs.size(); ++i) {
        stmt.bind_and_step((int)i, gen_repr[i], gen_degs[i].s, gen_degs[i].t);
    }
}

void DbAdamsSS::save_gen_names(const std::string& table_prefix, const std::vector<std::string>& gen_names) const
{
    Statement stmt(*this, "UPDATE " + table_prefix + "_generators SET name=?1 WHERE id=?2;");
    for (size_t i = 0; i < gen_names.size(); ++i)
        if (!gen_names[i].empty())
            stmt.bind_and_step(gen_names[i], (int)i);
}

void DbAdamsSS::save_gb(const std::string& table_prefix, const std::map<AdamsDeg, Poly1d>& gb) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (rel, s, t) VALUES (?1, ?2, ?3);");
    int count = 0;
    for (auto& [deg, polys] : gb) {
        for (auto& poly : polys) {
            stmt.bind_and_step(Serialize(poly), deg.s, deg.t);
            ++count;
        }
    }
}

void DbAdamsSS::save_gb_mod(const std::string& table_prefix, const std::map<AdamsDeg, Mod1d>& gbm) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (rel, s, t) VALUES (?1, ?2, ?3);");

    int count = 0;
    for (auto& [deg, xs] : gbm) {
        for (auto& x : xs) {
            stmt.bind_and_step(Serialize(x), deg.s, deg.t);
            ++count;
        }
    }
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

void DbAdamsSS::save_basis(const std::string& table_prefix, const std::map<alg2::AdamsDeg, alg2::Mon1d>& basis, const std::map<AdamsDeg, int2d>& repr) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (id, mon, repr, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.size(); ++i) {
            stmt.bind_and_step(count, Serialize(basis_d[i]), Serialize(repr.at(deg)[i]), deg.s, deg.t);
            ++count;
        }
    }
}

void DbAdamsSS::save_basis_mod(const std::string& table_prefix, const std::map<alg2::AdamsDeg, alg2::MMod1d>& basis, const std::map<AdamsDeg, int2d>& repr) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_basis (id, mon, repr, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.size(); ++i) {
            stmt.bind_and_step(count, Serialize(basis_d[i]), Serialize(repr.at(deg)[i]), deg.s, deg.t);
            ++count;
        }
    }
}

AdamsDeg1d DbAdamsSS::load_gen_adamsdegs(const std::string& table_prefix) const
{
    AdamsDeg1d result;
    Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators ORDER BY id;");
    while (stmt.step() == MYSQLITE_ROW)
        result.push_back({stmt.column_int(0), stmt.column_int(1)});
    return result;
}

std::vector<std::string> DbAdamsSS::load_gen_names(const std::string& table_prefix) const
{
    std::vector<std::string> result;
    Statement stmt(*this, "SELECT coalesce(name, \"\") FROM " + table_prefix + "_generators ORDER BY id;");
    while (stmt.step() == MYSQLITE_ROW)
        result.push_back(stmt.column_str(0));
    return result;
}

Poly1d DbAdamsSS::load_gb(const std::string& table_prefix, int t_max) const
{
    Poly1d result;
    Statement stmt(*this, "SELECT rel FROM " + table_prefix + "_relations" + (t_max == alg2::DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    while (stmt.step() == MYSQLITE_ROW) {
        Poly g = Deserialize<Poly>(stmt.column_str(0));
        result.push_back(std::move(g));
    }
    return result;
}

Mod1d DbAdamsSS::load_gb_mod(const std::string& table_prefix, int t_max) const
{
    Mod1d result;
    Statement stmt(*this, "SELECT rel FROM " + table_prefix + "_relations" + (t_max == alg2::DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    while (stmt.step() == MYSQLITE_ROW) {
        Mod g = Deserialize<Mod>(stmt.column_str(0));
        result.push_back(std::move(g));
    }
    return result;
}

BasisMon DbAdamsSS::load_basis(const std::string& table_prefix) const
{
    BasisMon result;
    Statement stmt(*this, "SELECT s, t, mon FROM " + table_prefix + "_basis ORDER BY id");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
        result[deg].push_back(Deserialize<Mon>(stmt.column_str(2)));
    }
    return result;
}

std::map<AdamsDeg, int2d> DbAdamsSS::load_basis_d2(const std::string& table_prefix) const
{
    std::map<AdamsDeg, int2d> result;
    if (has_column(table_prefix + "_basis", "d2")) {
        Statement stmt(*this, "SELECT s, t, d2 FROM " + table_prefix + "_basis WHERE d2 IS NOT NULL ORDER BY id");
        int count = 0;
        while (stmt.step() == MYSQLITE_ROW) {
            ++count;
            AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
            result[deg].push_back(Deserialize<int1d>(stmt.column_str(2)));
        }
    }
    return result;
}

BasisMMod DbAdamsSS::load_basis_mod(const std::string& table_prefix) const
{
    BasisMMod result;
    Statement stmt(*this, "SELECT s, t, mon FROM " + table_prefix + "_basis ORDER BY id");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
        result[deg].push_back(Deserialize<MMod>(stmt.column_str(2)));
    }
    return result;
}

void DbAdamsSS::save_pi_generators(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Poly1d& gen_Einf) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_generators (id, Einf, s, t) VALUES (?1, ?2, ?3, ?4);");

    for (size_t i = 0; i < gen_degs.size(); ++i) {
        stmt.bind_and_step((int)i, Serialize(gen_Einf[i]), gen_degs[i].s, gen_degs[i].t);
    }
}

void DbAdamsSS::save_pi_gb(const std::string& table_prefix, const std::map<AdamsDeg, algZ::Poly1d>& gb, const std::map<AdamsDeg, int2d>& gb_Einf) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_relations (rel, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");
    int count = 0;
    for (auto& [deg, polys] : gb) {
        auto& Einfs = gb_Einf.at(deg);
        for (size_t i = 0; i < polys.size(); ++i) {
            stmt.bind_and_step(Serialize(polys[i]), polys[i].Str(), Serialize(Einfs[i]), deg.s, deg.t);
            ++count;
        }
    }
}

void DbAdamsSS::save_pi_gb_mod(const std::string& table_prefix, const std::map<AdamsDeg, algZ::Mod1d>& gbm, const std::map<AdamsDeg, int2d>& gb_Einf) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_relations (rel, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");
    int count = 0;
    for (auto& [deg, polys] : gbm) {
        auto& Einfs = gb_Einf.at(deg);
        for (size_t i = 0; i < polys.size(); ++i) {
            stmt.bind_and_step(Serialize(polys[i]), polys[i].Str(), Serialize(Einfs[i]), deg.s, deg.t);
            ++count;
        }
    }
}

void DbAdamsSS::save_pi_basis(const std::string& table_prefix, const std::map<AdamsDeg, algZ::Mon1d>& basis, const std::map<AdamsDeg, int2d>& basis_Einf) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_basis (id, mon, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.size(); ++i) {
            stmt.bind_and_step(count, Serialize(basis_d[i]), basis_d[i].Str(), Serialize(basis_Einf.at(deg)[i]), deg.s, deg.t);
            ++count;
        }
    }
}

void DbAdamsSS::save_pi_basis_mod(const std::string& table_prefix, const std::map<AdamsDeg, algZ::MMod1d>& basis, const std::map<AdamsDeg, int2d>& basis_Einf) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_basis (id, mon, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.size(); ++i) {
            stmt.bind_and_step(count, Serialize(basis_d[i]), basis_d[i].Str(), Serialize(basis_Einf.at(deg)[i]), deg.s, deg.t);
            ++count;
        }
    }
}

algZ::Poly1d DbAdamsSS::load_pi_gb(const std::string& table_prefix, int t_max) const
{
    algZ::Poly1d result;
    Statement stmt(*this, "SELECT rel FROM " + table_prefix + "_pi_relations" + (t_max == DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    while (stmt.step() == MYSQLITE_ROW) {
        algZ::Poly g = Deserialize<algZ::Poly>(stmt.column_str(0));
        result.push_back(std::move(g));
    }
    return result;
}

algZ::Mod1d DbAdamsSS::load_pi_gb_mod(const std::string& table_prefix, int t_max) const
{
    algZ::Mod1d result;
    Statement stmt(*this, "SELECT rel FROM " + table_prefix + "_pi_relations" + (t_max == DEG_MAX ? "" : " WHERE t<=" + std::to_string(t_max)) + " ORDER BY t;");
    while (stmt.step() == MYSQLITE_ROW) {
        algZ::Mod g = Deserialize<algZ::Mod>(stmt.column_str(0));
        result.push_back(std::move(g));
    }
    return result;
}

std::map<AdamsDeg, algZ::Mon1d> DbAdamsSS::load_pi_basis(const std::string& table_prefix) const
{
    std::map<AdamsDeg, algZ::Mon1d> result;
    Statement stmt(*this, "SELECT s, t, mon FROM " + table_prefix + "_pi_basis ORDER BY id");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
        result[deg].push_back(Deserialize<algZ::Mon>(stmt.column_str(2)));
    }
    return result;
}

std::map<AdamsDeg, algZ::MMod1d> DbAdamsSS::load_pi_basis_mod(const std::string& table_prefix) const
{
    std::map<AdamsDeg, algZ::MMod1d> result;
    Statement stmt(*this, "SELECT s, t, mon FROM " + table_prefix + "_pi_basis ORDER BY id");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg deg = {stmt.column_int(0), stmt.column_int(1)};
        result[deg].push_back(Deserialize<algZ::MMod>(stmt.column_str(2)));
    }
    return result;
}

}  // namespace myio
