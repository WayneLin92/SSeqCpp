#include "main.h"

class DbMiagrate : public DBSS
{
public:
    DbMiagrate() = default;
    explicit DbMiagrate(const std::string& filename) : DBSS(filename) {}

    void create_ss_primitives(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss_primitives (id INTEGER PRIMARY KEY, base TEXT, diff TEXT, level SMALLINT, s SMALLINT, t SMALLINT);");
    }

    void drop_and_create_ss_primitives(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss_primitives");
        create_ss_primitives(table_prefix);
    }

    void save_ss_primitives(const std::string& table_prefix, const Staircases& basis_ss) const
    {
        myio::Statement stmt(*this, "INSERT INTO " + table_prefix + "_ss_primitives (id, base, diff, level, s, t) VALUES (?1, ?2, ?3, ?4, ?5, ?6);");
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
        std::cout << count << " basis_ss are inserted into " + table_prefix + "_ss_primitives.\n";
    }
};

/* Copy */
void migrate(const Diagram& ss1, Diagram& ss2, Staircases& primitives, int t_max_zero)
{
    const auto& basis_ss1 = ss1.GetS0().basis_ss;

    Timer timer(3600);
    ss2.DeduceDiffs(10, 10, 1, 1, timer);
    ss2.DeduceImageJ();

    AdamsDeg1d degs;
    for (auto& [d, _] : basis_ss1.front())
        degs.push_back(d);
    std::sort(degs.begin(), degs.end(), [](AdamsDeg d1, AdamsDeg d2) { return d1.stem() < d2.stem() || (d1.stem() == d2.stem() && d1.s < d2.s); });
    for (AdamsDeg deg : degs) {
        const auto& sc = ss1.GetRecentStaircase(ss1.GetS0().basis_ss, deg);
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            if (sc.levels[i] > kLevelMax / 2) {
                int rt = 0;
                if (sc.diffs_ind[i] != int1d{-1})
                    rt = ss2.SetS0DiffLeibnizV2(deg, sc.basis_ind[i], sc.diffs_ind[i], kLevelMax - sc.levels[i]);
                else if (sc.levels[i] > kLevelPC || deg.stem() <= t_max_zero)
                    rt = ss2.SetS0DiffLeibnizV2(deg, sc.basis_ind[i], int1d{}, kLevelMax - sc.levels[i] - 1);
                if (rt) {
                    primitives[deg].basis_ind.push_back(sc.basis_ind[i]);
                    primitives[deg].diffs_ind.push_back(sc.diffs_ind[i]);
                    primitives[deg].levels.push_back(sc.levels[i]);

                    Timer timer(3600);
                    ss2.DeduceDiffs(10, 10, 1, 1, timer);
                    ss2.DeduceImageJ();
                }
            }
        }
    }
}

int main_deduce_migrate(int argc, char** argv, int index)
{
    std::string db_in = "S0_AdamsSS_t245.db";
    std::string table_in = "AdamsE2";
    std::string db_out = DB_DEFAULT;
    std::string table_out = TABLE_DEFAULT;
    int t_max_zero = 381;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Deduce trivial differentials for degree reason\n";
        std::cout << "Usage:\n  ss deduce migrate <db_in> <table_in> <db_out> <table_out> [t_max_zero]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  db_in = " << db_in << "\n";
        std::cout << "  table_in = " << table_in << "\n";
        std::cout << "  db_out = " << db_out << "\n";
        std::cout << "  table_out = " << table_out << "\n";
        std::cout << "  t_max_zero = " << t_max_zero << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "db_in", db_in))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_in", table_in))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_out", db_out))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_out", table_out))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max_zero", t_max_zero))
        return index;

    bench::Timer timer;
    DBSS db1(db_in);
    DbMiagrate db2(db_out);

    Diagram diagram1({db_in}), diagram2({db_out});

    try {
        Staircases primitives;
        migrate(diagram1, diagram2, primitives, t_max_zero);

        db2.begin_transaction();
        db2.update_basis_ss(table_out, diagram2.GetChanges(0));
        db2.drop_and_create_ss_primitives(table_out);
        db2.save_ss_primitives(table_out, primitives);
        db2.end_transaction();
    }
    /*catch (SSException& e) {
        std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
    }*/
    catch (NoException&) {
        ;
    }

    return 0;
}