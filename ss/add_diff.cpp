#include "algebras/linalg.h"
#include "main.h"

using namespace alg2;

std::ostream& operator<<(std::ostream& sout, const int1d& arr)
{
    return sout << myio::TplStrCont("(", ",", ")", "()", arr.begin(), arr.end(), [](int i) { return std::to_string(i); });
}

std::ostream& operator<<(std::ostream& sout, const Staircase& sc)
{
    sout << "Staircase:\n";
    for (size_t i = 0; i < sc.levels.size(); ++i) {
        sout << sc.levels[i] << "  " << sc.basis_ind[i] << "  " << sc.diffs_ind[i] << '\n';
    }
    return sout;
}

int main_add_diff(int argc, char** argv, int index)
{
    int stem = 0, s = 0, r = 0;
    std::string x_str, dx_str;
    size_t iDb = 0;
    std::string mode = "add";
    std::string db_S0 = DB_S0;
    std::vector<std::string> dbnames = {
        DB_C2,
        DB_Ceta,
        DB_Cnu,
        DB_Csigma,
    };

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Manually input a differential into the ss table\n";
        std::cout << "Usage:\n  ss add_diff <stem> <s> <r> <x> <dx> [iDb] [mode:add/deduce/try] [db_S0] [db_Cofibs ...]\n\n";

        std::cout << "mode:\n";
        std::cout << "add: add differential and save\n";
        std::cout << "deduce: add differential, deduce and save\n";
        std::cout << "try: add differential, deduce and do not save\n";

        std::cout << "Default values:\n";
        std::cout << "  iDb = " << iDb << "\n";
        std::cout << "  mode = " << mode << "\n";
        std::cout << "  db_S0 = " << db_S0 << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "stem", stem))
        return index;
    if (myio::load_arg(argc, argv, ++index, "s", s))
        return index;
    if (myio::load_arg(argc, argv, ++index, "r", r))
        return index;
    if (myio::load_arg(argc, argv, ++index, "x", x_str))
        return index;
    if (myio::load_arg(argc, argv, ++index, "dx", dx_str))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "iDb", iDb))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "mode", mode))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_S0", db_S0))
        return index;
    if (myio::load_args(argc, argv, ++index, "db_Cofibs", dbnames))
        return index;
    dbnames.insert(dbnames.begin(), db_S0);

    AdamsDeg deg_x(s, stem + s);
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    Diagram diagram(dbnames);

    /* Check if x, dx are valid */
    std::sort(x.begin(), x.end());
    std::sort(dx.begin(), dx.end());
    auto& basis_ss = diagram.GetAllBasisSs()[iDb]->front();
    if (!x.empty()) {
        if (basis_ss.find(deg_x) == basis_ss.end()) {
            std::cout << "deg_x not found" << std::endl;
            return 101;
        }
        if (x.front() < 0 || basis_ss.at(deg_x).levels.size() <= (size_t)x.back()) {
            std::cout << "Invalid x";
            return 102;
        }
    }
    if (!dx.empty()) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
        if (basis_ss.find(deg_dx) == basis_ss.end()) {
            std::cout << "deg_dx not found" << std::endl;
            return 103;
        }
        if (dx.front() < 0 || basis_ss.at(deg_dx).levels.size() <= (size_t)dx.back()) {
            std::cout << "Invalid dx" << std::endl;
            return 104;
        }
    }

    int count = diagram.SetDiffLeibnizV2(iDb, deg_x, x, dx, r);
    if (count > 0 && (mode == "try" || mode == "deduce")) {
        Timer timer(300);
        diagram.DeduceDiffs(0, 1, timer);
    }

    if (mode == "add" || mode == "deduce") {
        diagram.ApplyChanges(1);
        for (size_t k = 0; k < dbnames.size(); ++k) {
            DBSS db(dbnames[k]);
            db.begin_transaction();
            db.update_basis_ss(GetE2TablePrefix(dbnames[k]), diagram.GetChanges(k));
            db.end_transaction();
        }
    }
    std::cout << "Done" << std::endl;
    return 0;
}
