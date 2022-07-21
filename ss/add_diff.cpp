#include "algebras/linalg.h"
#include "main.h"

using namespace alg;

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
    std::string db_filename = "AdamsE2Export_t220.db";
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Manually input a differential into the ss table\n";
        std::cout << "Usage:\n  ss add_diff <stem> <s> <r> <x> <dx> [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.0 (2022-7-11)" << std::endl;
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
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "table_prefix", table_prefix))
        return index;

    AdamsDeg deg_x(s, stem + s);
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    SSDB db(db_filename);
    SS ss = db.load_ss("AdamsE2");

    // TODO: add sanitizer check.
    ss.SetDiffLeibnizV2(deg_x, x, dx, r);

    db.begin_transaction();
    db.update_basis_ss("AdamsE2", ss.GetChanges());
    db.end_transaction();
    return 0;
}
