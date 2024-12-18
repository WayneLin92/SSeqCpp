#include "main.h"
#include "mylog.h"

using namespace alg2;

/* Order by (t, -s) */
template <typename T>
AdamsDeg1d OrderDegsV2(const T& cont)  //// To be deleted
{
    AdamsDeg1d result;
    for (auto [d, _] : cont) {
        result.push_back(d);
    }
    std::sort(result.begin(), result.end(), [](const AdamsDeg& d1, const AdamsDeg& d2) { return d1.t < d2.t || (d1.t == d2.t && d1.s > d2.s); });
    return result;
}

void DBSS::save_pi_generators_mod(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Mod1d& gen_Einf) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_generators (id, Einf, s, t) VALUES (?1, ?2, ?3, ?4);");
    for (size_t i = 0; i < gen_degs.size(); ++i)
        stmt.bind_and_step((int)i, myio::Serialize(gen_Einf[i]), gen_degs[i].s, gen_degs[i].t);
}

void DBSS::update_ss(const std::string& table_prefix, const SSNodes& nodes_ss) const
{
    MyException::Assert(nodes_ss.size() == 2, "update_ss(): Need nodes_ss.size() == 2");
    if (has_table(table_prefix + "_ss")) {
        auto indices = load_ss_indices(table_prefix);
        Statement stmt(*this, fmt::format("UPDATE {}_ss SET base=?1, diff=?2, level=?3 WHERE id=?4;", table_prefix));
        for (auto [deg, sc] : nodes_ss.back().items()) {
            for (size_t i = 0; i < sc.basis.size(); ++i)
                stmt.bind_and_step(myio::Serialize(sc.basis[i]), SerializeDiff(sc.diffs[i]), sc.levels[i], indices[deg] + (int)i);
        }
    }
    else
        save_ss(table_prefix, nodes_ss);
}

void DBSS::save_ss(const std::string& table_prefix, const SSNodes& nodes_ss) const
{
    drop_and_create_ss(table_prefix);
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_ss (id, base, diff, level, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");
    int count = 0;
    auto degs = nodes_ss.front().arr_degs();
    std::sort(degs.begin(), degs.end(), [](const AdamsDeg& d1, const AdamsDeg& d2) { return d1.t < d2.t || (d1.t == d2.t && d1.s > d2.s); });
    for (const auto& deg : degs) {
        const auto& sc = nodes_ss.GetRecentSc(deg);
        for (size_t i = 0; i < sc.basis.size(); ++i)
            stmt.bind_and_step(count++, myio::Serialize(sc.basis[i]), SerializeDiff(sc.diffs[i]), sc.levels[i], deg.s, deg.t);
    }
}

void DBSS::save_cofseq(const std::string& table, const CofSeq& cofseq, bool incremental) const
{
    Statement stmt(*this, "INSERT INTO " + table + " (iC, base, diff, level, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");
    if (incremental) {
        Statement stmt_del(*this, "DELETE FROM " + table + " WHERE iC=?1 AND s=?2 AND t=?3;");
        for (size_t iTri = 0; iTri < cofseq.nodes_cofseq.size(); ++iTri) {
            auto& ss = cofseq.nodes_cofseq[iTri].changes();
            auto degs = ss.arr_degs();
            for (const auto& deg : degs) {
                auto& sc = ss.at(deg);
                stmt_del.bind_and_step((int)iTri, deg.s, deg.t);
                for (size_t i = 0; i < sc.basis.size(); ++i)
                    stmt.bind_and_step((int)iTri, myio::Serialize(sc.basis[i]), SerializeDiff(sc.diffs[i]), sc.levels[i], deg.s, deg.t);
            }
        }
    }
    else {
        drop_and_create_cofseq(table);
        for (size_t iTri = 0; iTri < cofseq.nodes_cofseq.size(); ++iTri) {
            auto degs = cofseq.nodes_cofseq[iTri].front().arr_degs();
            for (const auto& deg : degs) {
                const auto& sc = cofseq.nodes_cofseq[iTri].GetRecentSc(deg);
                for (size_t i = 0; i < sc.basis.size(); ++i)
                    stmt.bind_and_step((int)iTri, myio::Serialize(sc.basis[i]), SerializeDiff(sc.diffs[i]), sc.levels[i], deg.s, deg.t);
            }
        }
    }
}

void DBSS::save_pi_basis(const std::string& table_prefix, const PiBasis& basis) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_basis (id, mon, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    int count = 0;
    auto degs = OrderDegsV2(basis);  ////
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.nodes_pi_basis.size(); ++i)
            stmt.bind_and_step(count++, myio::Serialize(basis_d.nodes_pi_basis[i]), basis_d.nodes_pi_basis[i].Str(), myio::Serialize(basis_d.Einf[i]), deg.s, deg.t);
    }
}

void DBSS::save_pi_basis_mod(const std::string& table_prefix, const PiBasisMod& basis) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_basis (id, mon, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.nodes_pi_basis.size(); ++i)
            stmt.bind_and_step(count++, myio::Serialize(basis_d.nodes_pi_basis[i]), basis_d.nodes_pi_basis[i].Str(), myio::Serialize(basis_d.Einf[i]), deg.s, deg.t);
    }
}

void DBSS::save_pi_def(const std::string& table_prefix, const std::vector<EnumDef>& pi_gen_defs, const std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_generators_def (id, def, map_ids, multipliers, mult_name, fils) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    for (size_t i = 0; i < pi_gen_defs.size(); ++i) {
        int1d map_ids;
        algZ::Poly mons;
        int1d fils;
        for (auto& gd : pi_gen_def_mons[i]) {
            map_ids.push_back(gd.map_index);
            mons.data.push_back(gd.m);
            fils.push_back(gd.O);
        }

        std::string name;
        if (mons) {
            name = mons.data[0].Str();
            for (size_t j = 1; j < mons.data.size(); ++j)
                name += "," + mons.data[j].Str();
        }
        stmt.bind_and_step((int)i, (int)pi_gen_defs[i], myio::Serialize(map_ids), myio::Serialize(mons), name, myio::Serialize(fils));
    }
}

std::map<AdamsDeg, int> DBSS::load_ss_indices(const std::string& table_prefix) const
{
    std::map<AdamsDeg, int> result;
    Statement stmt(*this, "SELECT s, t, min(id) FROM " + table_prefix + "_ss GROUP BY t, s;");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg d = {stmt.column_int(0), stmt.column_int(1)};
        result[d] = stmt.column_int(2);
    }
    return result;
}

SS DBSS::load_ss(const std::string& table_prefix) const
{
    SS nodes_ss;
    Statement stmt(*this, "SELECT base, COALESCE(diff, \"-1\"), level, s, t FROM " + table_prefix + "_ss;");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        int1d base = myio::Deserialize<int1d>(stmt.column_str(0));
        int1d diff = myio::Deserialize<int1d>(stmt.column_str(1));
        int level = stmt.column_int(2);
        AdamsDeg deg = {stmt.column_int(3), stmt.column_int(4)};
        auto& sc = nodes_ss[deg];
        sc.basis.push_back(std::move(base));
        sc.diffs.push_back(std::move(diff));
        sc.levels.push_back(level);
    }
    return nodes_ss;
}

std::array<SS, 3> DBSS::load_cofseq(const std::string& table) const
{
    std::array<SS, 3> node_cofseq;
    Statement stmt(*this, "SELECT iC, base, COALESCE(diff, \"-1\"), level, s, t FROM " + table);
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        size_t index = (size_t)stmt.column_int(0);
        int1d base = myio::Deserialize<int1d>(stmt.column_str(1));
        int1d diff = myio::Deserialize<int1d>(stmt.column_str(2));
        int level = stmt.column_int(3);
        AdamsDeg deg = {stmt.column_int(4), stmt.column_int(5)};

        node_cofseq[index][deg].diffs.push_back(std::move(diff));
        node_cofseq[index][deg].basis.push_back(std::move(base));
        node_cofseq[index][deg].levels.push_back(level);
    }
    return node_cofseq;
}

void DBSS::load_pi_def(const std::string& table_prefix, std::vector<EnumDef>& pi_gen_defs, std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const
{
    if (has_table(table_prefix + "_pi_generators_def")) {
        if (has_column(table_prefix + "_pi_generators_def", "map_ids")) {
            Statement stmt(*this, "SELECT def, map_ids, multipliers, fils FROM " + table_prefix + "_pi_generators_def order by id;");
            while (stmt.step() == MYSQLITE_ROW) {
                pi_gen_defs.push_back(EnumDef(stmt.column_int(0)));
                pi_gen_def_mons.push_back({});
                if (pi_gen_defs.back() == EnumDef::constraints) {
                    int1d map_ids = myio::Deserialize<int1d>(stmt.column_str(1));
                    algZ::Poly mons = myio::Deserialize<algZ::Poly>(stmt.column_str(2));
                    int1d fils = myio::Deserialize<int1d>(stmt.column_str(3));
                    for (size_t i = 0; i < map_ids.size(); ++i)
                        pi_gen_def_mons.back().push_back(GenConstraint{map_ids[i], mons.data[i], fils[i]});
                }
            }
        }
        else {  //// Compatibility
            std::string c;
            if (has_column(table_prefix + "_pi_generators_def", "property"))
                c = "property";
            else
                c = "mons";
            Statement stmt(*this, "SELECT def, " + c + " FROM " + table_prefix + "_pi_generators_def order by id;");
            while (stmt.step() == MYSQLITE_ROW) {
                pi_gen_defs.push_back(EnumDef(stmt.column_int(0)));
                pi_gen_def_mons.push_back({});
                if (pi_gen_defs.back() == EnumDef::constraints) {
                    algZ::Poly mons = myio::Deserialize<algZ::Poly>(stmt.column_str(1));
                    for (size_t i = 0; i < mons.data.size(); i += 2)
                        pi_gen_def_mons.back().push_back(GenConstraint{0, mons.data[i], mons.data[i + 1].fil()});
                }
            }
        }
    }
}

/* generate the table of the spectral sequence */
void generate_ss(const std::string& name, const std::string& path, bool isRing, int r = 2)
{
    using namespace alg2;

    DBSS db(path);
    std::string table_prefix = fmt::format("{}_AdamsE2", name);
    BasisMon basis = db.load_basis(table_prefix);
    SSNodes nodes_ss;
    nodes_ss.push_back({});
    nodes_ss.push_back({});
    SS& ss = nodes_ss[0];

    /* fill ss */
    for (AdamsDeg d : basis.degs()) {
        auto& basis_d = basis.at(d);
        for (int i = 0; i < (int)basis_d.size(); ++i) {
            ss[d].basis.push_back({i});
            ss[d].diffs.push_back({-1});
            ss[d].levels.push_back(LEVEL_MAX - r);
        }
    }

    if (isRing) {
        ss[AdamsDeg{0, 0}].diffs = {{}};
        ss[AdamsDeg{0, 0}].levels = {LEVEL_PERM};
    }

    /* insert into the database */
    db.begin_transaction();
    db.drop_and_create_ss(table_prefix);
    db.save_ss(table_prefix, nodes_ss);

    db.drop_and_create_pi_relations(name);
    db.drop_and_create_pi_basis(name);
    db.drop_and_create_pi_definitions(name);
    if (isRing)
        db.drop_and_create_pi_generators(name);
    else
        db.drop_and_create_pi_generators_mod(name);
    db.end_transaction();
}

int main_reset(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    bool confirm = 0;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {{"confirm", &confirm}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    fmt::print("Confirm to reset {}\n", cat_name);
    if (confirm || myio::UserConfirm()) {
        std::vector<std::string> names, paths;
        std::vector<int> isRing;
        GetAllDbNames(cat_name, names, paths, isRing);
        for (size_t i = 0; i < names.size(); ++i)
            generate_ss(names[i], paths[i], isRing[i]);
        for (const auto& entry : std::filesystem::directory_iterator(cat_name)) {
            auto filename = entry.path().filename().string();
            if (filename.starts_with("cofseq_") && filename.ends_with(".db"))
                std::filesystem::remove(entry.path());
        }
        auto path_log = fmt::format("{}/log.db", cat_name);
        Logger::SetOutDeduce(path_log.c_str());
        Logger::Reset();
        auto flag = SSFlag::cofseq;
        Category category(cat_name, "", flag, false, true);

        if (auto iCw = category.GetIndexCwByName("S0")) {
            auto& nodes_ss = category.GetRings()[iCw.index].nodes_ss;
            if (auto deg = AdamsDeg(25, 63 + 25); nodes_ss.front().has(deg)) {
                if (auto rt = category.SetCwDiffGlobal(iCw, deg, int1d{1}, int1d{0}, 5, false, flag))
                    return -127;
            }
            if (auto deg = AdamsDeg(56, 127 + 56); nodes_ss.front().has(deg)) {
                if (auto rt = category.SetCwDiffGlobal(iCw, deg, int1d{1}, int1d{0}, 6, false, flag))
                    return -253;
            }
        }
        if (auto iCw = category.GetIndexCwByName("tmf")) {
            auto& nodes_ss = category.GetRings()[iCw.index].nodes_ss;
            if (auto deg = AdamsDeg(16, 96 + 16); nodes_ss.front().has(deg)) {
                if (auto rt = category.SetCwDiffGlobal(iCw, deg, int1d{0}, int1d{0}, 3, false, flag))
                    return -325;
            }
        }
        if (auto rt = category.DeduceTrivialDiffs(flag))
            throw RunTimeError("DeduceTrivialDiffs() failed");
        category.SaveNodes(cat_name, "", true, flag);
    }

    return 0;
}

int main_reset_cofseq(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    fmt::print("Confirm to reset_cofseq {}\n", cat_name);
    if (myio::UserConfirm()) {
        for (const auto& entry : std::filesystem::directory_iterator(cat_name)) {
            auto filename = entry.path().filename().string();
            if (filename.starts_with("cofseq_") && filename.ends_with(".db"))
                std::filesystem::remove(entry.path());
        }
    }

    return 0;
}

int main_reset_pi(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    fmt::print("Confirm to reset_pi {}\n", cat_name);
    if (myio::UserConfirm()) {
        std::vector<std::string> names, paths;
        std::vector<int> isRing;
        GetAllDbNames(cat_name, names, paths, isRing);
        for (size_t i = 0; i < names.size(); ++i) {
            DBSS db(paths[i]);
            db.begin_transaction();
            db.drop_and_create_pi_relations(names[i]);
            db.drop_and_create_pi_basis(names[i]);
            db.drop_and_create_pi_definitions(names[i]);
            if (isRing[i])
                db.drop_and_create_pi_generators(names[i]);
            else
                db.drop_and_create_pi_generators_mod(names[i]);
            db.end_transaction();
        }
    }

    return 0;
}

int main_checkpoint(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name, ckpt_name;
    myio::CmdArg1d args = {{"category", &cat_name}, {"checkpoint", &ckpt_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::cofseq;
    Category category(cat_name, "", flag, true, true);
    category.SaveNodes(cat_name, ckpt_name, false, flag);

    return 0;
}

int main_rollback(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name, ckpt_name;
    myio::CmdArg1d args = {{"category", &cat_name}, {"checkpoint", &ckpt_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::cofseq;
    Category category(cat_name, ckpt_name , flag, true, true);
    category.SaveNodes(cat_name, "", false, flag); //// TODO: deal with exclusions
    return 0;
}
