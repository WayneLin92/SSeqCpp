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
    std::string db_filename = db_ss_default;
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Manually input a differential into the ss table\n";
        std::cout << "Usage:\n  ss add_diff <stem> <s> <r> <x> <dx> [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.0 (2022-7-21)" << std::endl;
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

    DBSS db(db_filename);
    SS ss = db.LoadSS("AdamsE2");

    /* Check if x, dx are valid */
    std::sort(x.begin(), x.end());
    std::sort(dx.begin(), dx.end());
    auto& basis = ss.GetBasis();
    if (!x.empty()) {
        if (basis.find(deg_x) == basis.end()) {
            std::cout << "deg_x not found";
            return 1;
        }
        if (x.front() < 0 || basis.at(deg_x).size() <= (size_t)x.back()) {
            std::cout << "Invalid x";
            return 2;
        }
    }
    if (!dx.empty()) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
        if (basis.find(deg_dx) == basis.end()) {
            std::cout << "deg_dx not found";
            return 3;
        }
        if (dx.front() < 0 || basis.at(deg_dx).size() <= (size_t)dx.back()) {
            std::cout << "Invalid dx";
            return 4;
        }
    }

    ss.SetDiffLeibnizV2(deg_x, x, dx, r);

    db.begin_transaction();
    db.update_basis_ss("AdamsE2", ss.GetChanges());
    db.end_transaction();
    return 0;
}

int main_try_add_diff(int argc, char** argv, int index)
{
    int stem = 0, s = 0, r = 0;
    std::string x_str, dx_str;
    std::string db_filename = db_ss_default;
    std::string table_prefix = "AdamsE2";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Try to input a differential into the ss table and detect contradictions without changing the database\n";
        std::cout << "Usage:\n  ss try_add_diff <stem> <s> <r> <x> <dx> [db_filename] [table_prefix]\n\n";

        std::cout << "Default values:\n";

        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  table_prefix = " << table_prefix << "\n\n";

        std::cout << "Version:\n  1.0 (2022-7-21)" << std::endl;
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

    DBSS db(db_filename);
    SS ss = db.LoadSS("AdamsE2");

    /* Check if x, dx are valid */
    std::sort(x.begin(), x.end());
    std::sort(dx.begin(), dx.end());
    auto& basis = ss.GetBasis();
    if (!x.empty()) {
        if (basis.find(deg_x) == basis.end()) {
            std::cout << "deg_x not found";
            return 1;
        }
        if (x.front() < 0 || basis.at(deg_x).size() <= (size_t)x.back()) {
            std::cout << "Invalid x";
            return 2;
        }
    }
    if (!dx.empty()) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
        if (basis.find(deg_dx) == basis.end()) {
            std::cout << "deg_dx not found";
            return 3;
        }
        if (dx.front() < 0 || basis.at(deg_dx).size() <= (size_t)dx.back()) {
            std::cout << "Invalid dx";
            return 4;
        }
    }

    try {
        ss.SetDiffLeibnizV2(deg_x, x, dx, r);
        Timer timer(600);
        int count = ss.DeduceDiffs(10, 10, 1, 1, timer);
        std::cout << "count=" << count << "                \n";
    }
    catch (SSException& e) {
        std::cerr << "SSException " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    catch (MyException& e) {
        std::cerr << "MyException " << std::hex << e.id() << ": " << e.what() << '\n';
    }
    /*catch (NoException&) {
        ;
    }*/
    return 0;
}