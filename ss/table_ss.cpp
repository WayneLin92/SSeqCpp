#include "main.h"

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

    myio::Logger::out() << "Generators inserted into " + table_prefix + "_pi_generators, size=" << gen_degs.size() << '\n';
}

void DBSS::save_basis_ss(const std::string& table_prefix, const Staircases& basis_ss) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_ss (id, base, diff, level, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");
    int count = 0;
    auto degs = OrderDegsV2(basis_ss);
    for (const auto& deg : degs) {
        auto& basis_ss_d = basis_ss.at(deg);
        for (size_t i = 0; i < basis_ss_d.basis_ind.size(); ++i) {
            stmt.bind_int(1, count++);
            stmt.bind_str(2, myio::Serialize(basis_ss_d.basis_ind[i]));
            if (basis_ss_d.diffs_ind[i] == int1d{-1})
                stmt.bind_null(3);
            else
                stmt.bind_str(3, myio::Serialize(basis_ss_d.diffs_ind[i]));
            stmt.bind_int(4, basis_ss_d.levels[i]);
            stmt.bind_int(5, deg.s);
            stmt.bind_int(6, deg.t);
            stmt.step_and_reset();
        }
    }
    myio::Logger::out() << "basis_ss are inserted into " << table_prefix << "_ss, size=" << count << '\n';
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

    myio::Logger::out() << count << " bases are inserted into " + table_prefix + "_pi_basis!\n";
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

    myio::Logger::out() << count << " bases are inserted into " + table_prefix + "_pi_basis!\n";
}

void DBSS::save_pi_def(const std::string& table_prefix, const std::vector<DefFlag>& pi_gen_defs, const std::vector<std::set<algZ::Mon>>& pi_gen_def_mons) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_generators_def (id, def, mons) VALUES (?1, ?2, ?3);");

    for (size_t i = 0; i < pi_gen_defs.size(); ++i) {
        stmt.bind_int(1, (int)i);
        stmt.bind_int(2, (int)pi_gen_defs[i]);
        algZ::Poly mons;
        mons.data.insert(mons.data.end(), pi_gen_def_mons[i].begin(), pi_gen_def_mons[i].end());
        stmt.bind_str(3, myio::Serialize(mons));
        stmt.step_and_reset();
    }

    myio::Logger::out() << pi_gen_defs.size() << " definitions are inserted into " + table_prefix + "_pi_generators_def!\n";
}

void DBSS::save_pi_def_mod(const std::string& table_prefix, const std::vector<DefFlag>& pi_gen_defs, const std::vector<std::set<algZ::MMod>>& pi_gen_def_mons) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_pi_generators_def (id, def, mons) VALUES (?1, ?2, ?3);");

    for (size_t i = 0; i < pi_gen_defs.size(); ++i) {
        stmt.bind_int(1, (int)i);
        stmt.bind_int(2, (int)pi_gen_defs[i]);
        algZ::Mod mons;
        mons.data.insert(mons.data.end(), pi_gen_def_mons[i].begin(), pi_gen_def_mons[i].end());
        stmt.bind_str(3, myio::Serialize(mons));
        stmt.step_and_reset();
    }

    myio::Logger::out() << pi_gen_defs.size() << " definitions are inserted into " + table_prefix + "_pi_generators_def!\n";
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
    myio::Logger::out() << "Indices loaded from " << table_prefix << "_ss, size=" << count << '\n';
    return result;
}

void DBSS::update_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& basis_ss) const
{
    std::map<AdamsDeg, int> indices = load_basis_indices(table_prefix);
    Statement stmt(*this, "UPDATE " + table_prefix + "_ss SET base=?1, diff=?2, level=?3 WHERE id=?4;");

    int count = 0;
    for (const auto& [deg, basis_ss_d] : basis_ss) {
        for (size_t i = 0; i < basis_ss_d.basis_ind.size(); ++i) {
            stmt.bind_str(1, myio::Serialize(basis_ss_d.basis_ind[i]));
            if (basis_ss_d.diffs_ind[i] == int1d{-1})
                stmt.bind_null(2);
            else
                stmt.bind_str(2, myio::Serialize(basis_ss_d.diffs_ind[i]));
            stmt.bind_int(3, basis_ss_d.levels[i]);
            stmt.bind_int(4, indices[deg] + (int)i);
            stmt.step_and_reset();
            ++count;
        }
    }
    myio::Logger::out() << table_prefix + "_ss is updated, num_of_change=" << count << '\n';
}

Staircases DBSS::load_basis_ss(const std::string& table_prefix) const
{
    Staircases basis_ss;
    Statement stmt(*this, "SELECT base, COALESCE(diff, \"-1\"), level, s, t FROM " + table_prefix + "_ss;");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        int1d base = myio::Deserialize<int1d>(stmt.column_str(0));
        int level = stmt.column_int(2);
        AdamsDeg deg = {stmt.column_int(3), stmt.column_int(4)};

        basis_ss[deg].basis_ind.push_back(std::move(base));
        basis_ss[deg].levels.push_back(level);

        int1d diff = myio::Deserialize<int1d>(stmt.column_str(1));
        basis_ss[deg].diffs_ind.push_back(std::move(diff));
    }
    myio::Logger::out() << "basis_ss loaded from " << table_prefix << "_ss, size=" << count << '\n';
    return basis_ss;
}

void DBSS::load_pi_def(const std::string& table_prefix, std::vector<DefFlag>& pi_gen_defs, std::vector<std::set<algZ::Mon>>& pi_gen_def_mons) const
{
    Statement stmt(*this, "SELECT def, mons FROM " + table_prefix + "_pi_generators_def order by id;");
    while (stmt.step() == MYSQLITE_ROW) {
        pi_gen_defs.push_back(DefFlag(stmt.column_int(0)));
        pi_gen_def_mons.push_back({});
        if (pi_gen_defs.back() == DefFlag::mon) {
            algZ::Poly mons = myio::Deserialize<algZ::Poly>(stmt.column_str(1));
            pi_gen_def_mons.back().insert(mons.data.begin(), mons.data.end());
        }
    }
    myio::Logger::out() << "Definitions loaded from " << table_prefix << "_pi_generators_def, size=" << pi_gen_defs.size() << '\n';
}

void DBSS::load_pi_def_mod(const std::string& table_prefix, std::vector<DefFlag>& pi_gen_defs, std::vector<std::set<algZ::MMod>>& pi_gen_def_mons) const
{
    Statement stmt(*this, "SELECT def, mons FROM " + table_prefix + "_pi_generators_def order by id;");
    while (stmt.step() == MYSQLITE_ROW) {
        pi_gen_defs.push_back(DefFlag(stmt.column_int(0)));
        pi_gen_def_mons.push_back({});
        if (pi_gen_defs.back() == DefFlag::mon) {
            algZ::Mod mons = myio::Deserialize<algZ::Mod>(stmt.column_str(1));
            pi_gen_def_mons.back().insert(mons.data.begin(), mons.data.end());
        }
    }
    myio::Logger::out() << "Definitions loaded from " << table_prefix << "_pi_generators_def, size=" << pi_gen_defs.size() << '\n';
}

/* generate the table of the spectral sequence */
void generate_ss(const std::string& db_filename, int r)
{
    using namespace alg2;

    DBSS db(db_filename);
    std::string table_prefix = GetE2TablePrefix(db_filename);
    std::map<AdamsDeg, Mon1d> basis = db.load_basis(table_prefix);
    std::map<AdamsDeg, Staircase> basis_ss;

    /* fill basis_ss */
    int prev_t = 0;
    for (auto& [d, basis_d] : basis) {
        for (int i = 0; i < (int)basis_d.size(); ++i) {
            basis_ss[d].basis_ind.push_back({i});
            basis_ss[d].diffs_ind.push_back({-1});
            basis_ss[d].levels.push_back(kLevelMax - r);
        }
    }

    basis_ss[AdamsDeg{0, 0}].diffs_ind = {{}};
    basis_ss[AdamsDeg{0, 0}].levels = {kLevelMax / 2};

    /* insert into the database */
    db.begin_transaction();
    db.drop_and_create_basis_ss(table_prefix);
    db.save_basis_ss(table_prefix, basis_ss);

    auto pi_table = GetComplexName(db_filename);
    db.begin_transaction();
    db.drop_and_create_pi_relations(pi_table);
    db.drop_and_create_pi_basis(pi_table);
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

    for (size_t k = 0; k < dbnames.size(); ++k) {
        generate_ss(dbnames[k], 2);
    }

    Diagram diagram(dbnames);

    std::cout << "Confirm to reset\n";
    if (myio::UserConfirm()) {
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

    std::cout << "Confirm to resetpi\n";
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

int main_resetfrom(int argc, char** argv, int index)
{
    std::string selector = "debug";
    std::string selector_from = "debug-reset";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Initialize the homotopy data\n";
        std::cout << "Usage:\n  ss resetfrom [selector] [selector_from]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";
        std::cout << "  selector_from = " << selector_from << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);
    auto dbnames_from = GetDbNames(selector_from);

    for (size_t k = 0; k < dbnames.size(); ++k) {
        std::ifstream src(dbnames_from[k], std::ios::binary);
        std::ofstream dst(dbnames[k], std::ios::binary);
        dst << src.rdbuf();
    }

    return 0;
}

int main_truncate(int argc, char** argv, int index)
{
    std::string selector = "default";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Debugging\n";
        std::cout << "Usage:\n  ss truncate [selector]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  selector = " << selector << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    for (size_t k = 0; k < dbnames.size(); ++k) {
        DBSS db(dbnames[k]);
        auto table = GetE2TablePrefix(dbnames[k]);
        auto basis_ss = db.load_basis_ss(table);
        int t_max = basis_ss.rbegin()->first.t;

        for (auto& [d, sc] : basis_ss) {
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] > kLevelPC) {
                    int r = kLevelMax - sc.levels[i];
                    if (d.t + r - 1 > t_max && sc.diffs_ind[i] != int1d{-1})
                        sc.diffs_ind[i] = int1d{-1};
                }
            }
        }

        db.begin_transaction();
        db.update_basis_ss(table, basis_ss);
        db.end_transaction();
    }

    return 0;
}
