#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"

using namespace alg2;

std::ostream& operator<<(std::ostream& sout, const int1d& arr)
{
    return sout << myio::TplStrCont("(", ",", ")", "()", arr.begin(), arr.end(), [](int i) { return std::to_string(i); });
}

std::ostream& operator<<(std::ostream& sout, const Staircase& sc)
{
    sout << "Staircase:\n";
    for (size_t i = 0; i < sc.levels.size(); ++i) {
        sout << sc.levels[i] << "  " << sc.basis[i] << "  " << sc.diffs[i] << '\n';
    }
    return sout;
}

int main_add_diff(int argc, char** argv, int index)
{
    int stem = 0, s = 0, r = 0;
    std::string x_str, dx_str;
    size_t iDb = 0;
    std::string mode = "add";
    std::string selector = "default";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Manually input a differential into the ss table\n";
        std::cout << "Usage:\n  ss add_diff <stem> <s> <r> <x> <dx> [iDb] [mode:add/deduce/try] [selector]\n\n";

        std::cout << "mode:\n";
        std::cout << "add: add differential and save\n";
        std::cout << "deduce: add differential, deduce and save\n";
        std::cout << "try: add differential, deduce and do not save\n";

        std::cout << "Default values:\n";
        std::cout << "  iDb = " << iDb << "\n";
        std::cout << "  mode = " << mode << "\n";
        std::cout << "  selector = " << selector << "\n\n";

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
    if (myio::load_op_arg(argc, argv, ++index, "selector", selector))
        return index;
    auto dbnames = GetDbNames(selector);

    AdamsDeg deg_x(s, stem + s);
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    Diagram diagram(dbnames, DeduceFlag::no_op);

    /* Check if x, dx are valid */
    std::sort(x.begin(), x.end());
    std::sort(dx.begin(), dx.end());
    auto& nodes_ss = diagram.GetAllBasisSs()[iDb]->front();
    if (!x.empty()) {
        if (nodes_ss.find(deg_x) == nodes_ss.end()) {
            std::cout << "deg_x not found" << std::endl;
            return 101;
        }
        if (x.front() < 0 || nodes_ss.at(deg_x).levels.size() <= (size_t)x.back()) {
            std::cout << "Invalid x";
            return 102;
        }
    }
    if (!dx.empty()) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
        if (nodes_ss.find(deg_dx) == nodes_ss.end()) {
            std::cout << "deg_dx not found" << std::endl;
            return 103;
        }
        if (dx.front() < 0 || nodes_ss.at(deg_dx).levels.size() <= (size_t)dx.back()) {
            std::cout << "Invalid dx" << std::endl;
            return 104;
        }
    }

    int count = diagram.SetDiffGlobal(iDb, deg_x, x, dx, r);
    if (count > 0 && (mode == "try" || mode == "deduce")) {
        Timer timer(300);
        diagram.DeduceDiffs(0, 500, 0, DeduceFlag::no_op);////
    }

    if (mode == "add" || mode == "deduce") {
        Logger::LogSummary("Changed differentials", count);
        diagram.save(dbnames, DeduceFlag::no_op);
    }
    std::cout << "Done" << std::endl;
    return 0;
}
