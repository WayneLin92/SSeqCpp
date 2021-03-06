#include "main.h"

using namespace alg;

void SSDB::save_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& basis_ss) const
{
    Statement stmt(*this, "INSERT INTO " + table_prefix + "_ss (id, base, diff, level, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");
    int count = 0;
    for (const auto& [deg, basis_ss_d] : basis_ss) {
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
    std::cout << count << " basis_ss are inserted into " + table_prefix + "_ss.\n";
}

std::map<AdamsDeg, int> SSDB::load_indices(const std::string& table_prefix) const
{
    std::map<AdamsDeg, int> result;
    Statement stmt(*this, "SELECT s, t, min(id) FROM " + table_prefix + "_ss GROUP BY t, s;");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        AdamsDeg d = {stmt.column_int(0), stmt.column_int(1)};
        result[d] = stmt.column_int(2);
    }
    std::cout << count << " indices loaded from " << table_prefix << "_ss.\n";
    return result;
}

void SSDB::update_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& basis_ss) const
{
    std::map<AdamsDeg, int> indices = load_indices(table_prefix);
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
    std::cout << table_prefix + "_ss is updated, num_of_change=" << count << '\n';
}

Staircases SSDB::load_basis_ss(const std::string& table_prefix) const
{
    Staircases basis_ss;
    Statement stmt(*this, "SELECT base, diff, level, s, t FROM " + table_prefix + "_ss;");
    int count = 0;
    while (stmt.step() == MYSQLITE_ROW) {
        ++count;
        int1d base = myio::Deserialize<int1d>(stmt.column_str(0));
        int level = stmt.column_int(2);
        AdamsDeg deg = {stmt.column_int(3), stmt.column_int(4)};

        basis_ss[deg].basis_ind.push_back(std::move(base));
        basis_ss[deg].levels.push_back(level);
        if (stmt.column_type(1) == MYSQLITE3_TEXT) {
            int1d diff = myio::Deserialize<int1d>(stmt.column_str(1));
            basis_ss[deg].diffs_ind.push_back(diff);
        }
        else
            basis_ss[deg].diffs_ind.push_back({-1});
    }
    std::cout << "basis_ss loaded from " << table_prefix << "_ss, size=" << count << '\n';
    return basis_ss;
}

SS SSDB::load_ss(const std::string& table_prefix) const
{
    Staircases basis_ss = load_basis_ss("AdamsE2");
    Poly1d polys = load_gb("AdamsE2", DEG_MAX);
    auto basis = load_basis("AdamsE2");
    Groebner gb(basis.rbegin()->first.t, {}, std::move(polys));
    return SS(std::move(gb), std::move(basis), std::move(basis_ss));
}

/* generate the table of the spectral sequence */
void generate_ss(const std::string& db_filename, const std::string& table_prefix, int r)
{
    using namespace alg;

    SSDB db(db_filename);
    std::map<AdamsDeg, Mon1d> basis = db.load_basis(table_prefix);

    int s_diff = 2;
    int t_diff = 0;
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
    db.create_basis_ss_and_delete(table_prefix);
    db.save_basis_ss(table_prefix, basis_ss);
    db.end_transaction();
}

int main_generate_ss(int argc, char** argv, int index)
{
    std::string db_filename = "AdamsE2Export_t220.db";
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Generate the ss table\n";
        std::cout << "Usage:\n  ss ss [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.0 (2022-7-11)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_prefix", table_prefix))
        return index;

    generate_ss(db_filename, table_prefix, 2);
    return 0;
}
