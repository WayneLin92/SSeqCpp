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
    std::string cw;
    std::string diagram_name = "default";
    std::string mode = "add";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Manually input a differential into the ss table\n";
        std::cout << "Usage:\n  ss add_diff <cw> <stem> <s> <r> <x> <dx> [diagram] [mode:add/deduce/try]\n\n";

        std::cout << "mode:\n";
        std::cout << "add: add differential and save\n";
        std::cout << "deduce: add differential, deduce and save\n";
        std::cout << "try: add differential, deduce and do not save\n";

        std::cout << "Default values:\n";
        std::cout << "  mode = " << mode << "\n";
        std::cout << "  diagram = " << diagram_name << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cw", cw))
        return index;
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
    if (myio::load_op_arg(argc, argv, ++index, "selector", diagram_name))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "mode", mode))
        return index;

    AdamsDeg deg_x(s, stem + s);
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    Diagram diagram(diagram_name, DeduceFlag::no_op);

    /* #Check if x, dx are valid */
    std::sort(x.begin(), x.end());
    std::sort(dx.begin(), dx.end());
    bool isRing = diagram.GetRingIndexByName(cw) != -1;
    auto& ss = isRing ? diagram.GetRingByName(cw).nodes_ss.front() : diagram.GetModuleByName(cw).nodes_ss.front();
    if (!x.empty()) {
        if (!ut::has(ss, deg_x)) {
            fmt::print("deg_x not found\n");
            return -101;
        }
        if (x.front() < 0 || ss.at(deg_x).levels.size() <= (size_t)x.back()) {
            fmt::print("Invalid x\n");
            return -102;
        }
    }
    if (!dx.empty()) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
        if (ss.find(deg_dx) == ss.end()) {
            fmt::print("deg_dx not found\n");
            return -103;
        }
        if (dx.front() < 0 || ss.at(deg_dx).levels.size() <= (size_t)dx.back()) {
            fmt::print("Invalid dx\n");
            return -104;
        }
    }

    /* #Add diff */
    int count = 0;
    if (isRing) {
        size_t iRing = (size_t)diagram.GetRingIndexByName(cw);
        count += diagram.SetRingDiffGlobal(iRing, deg_x, x, dx, r);
    }
    else {
        size_t iMod = (size_t)diagram.GetModuleIndexByName(cw);
        count += diagram.SetModuleDiffGlobal(iMod, deg_x, x, dx, r);
    }
    if (count > 0 && (mode == "try" || mode == "deduce")) {
        Logger::LogSummary("Changed differentials", count);
        diagram.DeduceDiffs(0, 500, 0, DeduceFlag::no_op);////
    }
    if (mode == "add" || mode == "deduce") {
        Logger::LogSummary("Changed differentials", count);
        diagram.save(diagram_name, DeduceFlag::no_op);
    }
    fmt::print("Done\n");
    return 0;
}
