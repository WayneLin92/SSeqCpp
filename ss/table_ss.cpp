#include "main.h"
#include <fmt/ranges.h>  ///

using namespace alg2;

void DBSS::save_pi_generators_mod(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Mod1d& gen_Einf, const algZ::Poly1d& to_S0) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_generators (id, Einf, to_S0, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");

    for (size_t i = 0; i < gen_degs.size(); ++i) {
        stmt.bind_int(1, (int)i);
        stmt.bind_str(2, myio::Serialize(gen_Einf[i]));
        stmt.bind_str(3, myio::Serialize(to_S0[i]));
        stmt.bind_int(4, gen_degs[i].s);
        stmt.bind_int(5, gen_degs[i].t);
        stmt.step_and_reset();
    }

    // myio::Logger::out() << "Generators inserted into " + table_prefix + "_pi_generators, size=" << gen_degs.size() << '\n';
}

void DBSS::save_basis_ss(const std::string& table_prefix, const Staircases& nodes_ss) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_ss (id, base, diff, level, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");
    int count = 0;
    auto degs = OrderDegsV2(nodes_ss);
    for (const auto& deg : degs) {
        auto& basis_ss_d = nodes_ss.at(deg);
        for (size_t i = 0; i < basis_ss_d.basis.size(); ++i) {
            stmt.bind_int(1, count++);
            stmt.bind_str(2, myio::Serialize(basis_ss_d.basis[i]));
            if (basis_ss_d.diffs[i] == NULL_DIFF)
                stmt.bind_null(3);
            else
                stmt.bind_str(3, myio::Serialize(basis_ss_d.diffs[i]));
            stmt.bind_int(4, basis_ss_d.levels[i]);
            stmt.bind_int(5, deg.s);
            stmt.bind_int(6, deg.t);
            stmt.step_and_reset();
        }
    }
    // myio::Logger::out() << "basis_ss are inserted into " << table_prefix << "_ss, size=" << count << '\n';
}

void DBSS::save_pi_basis(const std::string& table_prefix, const PiBasis& basis) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_basis (id, mon, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.pi_basis.size(); ++i) {
            stmt.bind_int(1, count);
            stmt.bind_str(2, myio::Serialize(basis_d.pi_basis[i]));
            stmt.bind_str(3, basis_d.pi_basis[i].Str());
            stmt.bind_str(4, myio::Serialize(basis_d.Einf[i]));
            stmt.bind_int(5, deg.s);
            stmt.bind_int(6, deg.t);
            stmt.step_and_reset();
            ++count;
        }
    }

    // myio::Logger::out() << count << " bases are inserted into " + table_prefix + "_pi_basis!\n";
}

void DBSS::save_pi_basis_mod(const std::string& table_prefix, const PiBasisMod& basis) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_basis (id, mon, name, Einf, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    int count = 0;
    auto degs = OrderDegsV2(basis);
    for (AdamsDeg deg : degs) {
        auto& basis_d = basis.at(deg);
        for (size_t i = 0; i < basis_d.pi_basis.size(); ++i) {
            stmt.bind_int(1, count);
            stmt.bind_str(2, myio::Serialize(basis_d.pi_basis[i]));
            stmt.bind_str(3, basis_d.pi_basis[i].Str());
            stmt.bind_str(4, myio::Serialize(basis_d.Einf[i]));
            stmt.bind_int(5, deg.s);
            stmt.bind_int(6, deg.t);
            stmt.step_and_reset();
            ++count;
        }
    }

    // myio::Logger::out() << count << " bases are inserted into " + table_prefix + "_pi_basis!\n";
}

void DBSS::save_pi_def(const std::string& table_prefix, const std::vector<EnumDef>& pi_gen_defs, const std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_generators_def (id, def, map_ids, multipliers, mult_name, fils) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");

    for (size_t i = 0; i < pi_gen_defs.size(); ++i) {
        stmt.bind_int(1, (int)i);
        stmt.bind_int(2, (int)pi_gen_defs[i]);

        int1d map_ids;
        algZ::Poly mons;
        int1d fils;
        for (auto& gd : pi_gen_def_mons[i]) {
            map_ids.push_back(gd.map_id);
            mons.data.push_back(gd.m);
            fils.push_back(gd.O);
        }
        stmt.bind_str(3, myio::Serialize(map_ids));
        stmt.bind_str(4, myio::Serialize(mons));
        if (mons) {
            std::string name = mons.data[0].Str();
            for (size_t i = 1; i < mons.data.size(); ++i)
                name += "," + mons.data[i].Str();
            stmt.bind_str(5, name);
        }
        else
            stmt.bind_str(5, "");
        stmt.bind_str(6, myio::Serialize(fils));
        stmt.step_and_reset();
    }

    // myio::Logger::out() << pi_gen_defs.size() << " definitions are inserted into " + table_prefix + "_pi_generators_def!\n";
}

std::map<AdamsDeg, int> DBSS::load_basis_indices(const std::string& table_prefix) const
{
    std::map<AdamsDeg, int> result;
    Statement stmt(*this, "SELECT s, t, min(id) FROM " + table_prefix + "_basis GROUP BY t, s;");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg d = {stmt.column_int(0), stmt.column_int(1)};
        result[d] = stmt.column_int(2);
    }
    // myio::Logger::out() << "Indices loaded from " << table_prefix << "_ss, size=" << count << '\n';
    return result;
}

void DBSS::update_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& nodes_ss) const
{
    std::map<AdamsDeg, int> indices = load_basis_indices(table_prefix);
    Statement stmt(*this, "UPDATE " + table_prefix + "_ss SET base=?1, diff=?2, level=?3 WHERE id=?4;");

    int count = 0;
    for (const auto& [deg, basis_ss_d] : nodes_ss) {
        for (size_t i = 0; i < basis_ss_d.basis.size(); ++i) {
            stmt.bind_str(1, myio::Serialize(basis_ss_d.basis[i]));
            if (basis_ss_d.diffs[i] == NULL_DIFF)
                stmt.bind_null(2);
            else
                stmt.bind_str(2, myio::Serialize(basis_ss_d.diffs[i]));
            stmt.bind_int(3, basis_ss_d.levels[i]);
            stmt.bind_int(4, indices[deg] + (int)i);
            stmt.step_and_reset();
            ++count;
        }
    }
    // myio::Logger::out() << table_prefix + "_ss is updated, num_of_change=" << count << '\n';
}

Staircases DBSS::load_basis_ss(const std::string& table_prefix) const
{
    Staircases nodes_ss;
    Statement stmt(*this, "SELECT base, COALESCE(diff, \"-1\"), level, s, t FROM " + table_prefix + "_ss;");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        int1d base = myio::Deserialize<int1d>(stmt.column_str(0));
        int level = stmt.column_int(2);
        AdamsDeg deg = {stmt.column_int(3), stmt.column_int(4)};

        nodes_ss[deg].basis.push_back(std::move(base));
        nodes_ss[deg].levels.push_back(level);

        int1d diff = myio::Deserialize<int1d>(stmt.column_str(1));
        nodes_ss[deg].diffs.push_back(std::move(diff));
    }
    // myio::Logger::out() << "basis_ss loaded from " << table_prefix << "_ss, size=" << count << '\n';
    return nodes_ss;
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
        // myio::Logger::out() << "Definitions loaded from " << table_prefix << "_pi_generators_def, size=" << pi_gen_defs.size() << '\n';
    }
}

/* generate the table of the spectral sequence */
void generate_ss(const std::string& db_filename, int r)
{
    using namespace alg2;

    DBSS db(db_filename);
    std::string table_prefix = GetE2TablePrefix(db_filename);
    std::map<AdamsDeg, Mon1d> basis = db.load_basis(table_prefix);
    std::map<AdamsDeg, Staircase> nodes_ss;

    /* fill nodes_ss */
    int prev_t = 0;
    for (auto& [d, basis_d] : basis) {
        for (int i = 0; i < (int)basis_d.size(); ++i) {
            nodes_ss[d].basis.push_back({i});
            nodes_ss[d].diffs.push_back({-1});
            nodes_ss[d].levels.push_back(LEVEL_MAX - r);
        }
    }

    if (GetComplexName(db_filename) == "S0") {  ////if it is a ring
        nodes_ss[AdamsDeg{0, 0}].diffs = {{}};
        nodes_ss[AdamsDeg{0, 0}].levels = {LEVEL_MAX / 2};
    }

    /* insert into the database */
    db.begin_transaction();
    db.drop_and_create_basis_ss(table_prefix);
    db.save_basis_ss(table_prefix, nodes_ss);

    auto pi_table = GetComplexName(db_filename);
    db.begin_transaction();
    db.drop_and_create_pi_relations(pi_table);
    db.drop_and_create_pi_basis(pi_table);
    db.drop_and_create_pi_definitions(pi_table);
    if (pi_table == "S0")
        db.drop_and_create_pi_generators(pi_table);
    else
        db.drop_and_create_pi_generators_mod(pi_table);
    db.end_transaction();
}

int main_reset(int argc, char** argv, int index)
{
    std::string selector = "default";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Initialize the ss table\n";
        std::cout << "Usage:\n  ss reset [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    std::cout << "Confirm to reset " << selector << '\n';
    if (myio::UserConfirm()) {
        for (size_t k = 0; k < dbnames.size(); ++k)
            generate_ss(dbnames[k], 2);

        Diagram diagram(dbnames, DeduceFlag::no_op);

        int count = diagram.DeduceTrivialDiffs();
        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(GetE2TablePrefix(dbnames[k]), diagram.GetBasisSSChanges(k));
            db.end_transaction();
        }
    }

    return 0;
}

int main_resetpi(int argc, char** argv, int index)
{
    std::string selector = "default";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Initialize the homotopy data\n";
        std::cout << "Usage:\n  ss resetpi [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    std::cout << "Confirm to resetpi " << selector << '\n';
    if (myio::UserConfirm()) {
        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            auto pi_table = GetComplexName(dbnames[k]);
            db.begin_transaction();
            db.drop_and_create_pi_relations(pi_table);
            db.drop_and_create_pi_basis(pi_table);
            db.drop_and_create_pi_definitions(pi_table);
            if (k == 0)
                db.drop_and_create_pi_generators(pi_table);
            else
                db.drop_and_create_pi_generators_mod(pi_table);
            db.end_transaction();
        }
    }

    return 0;
}
