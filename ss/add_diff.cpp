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

int main_add_diff(int argc, char** argv, int& index, const char* desc)
{
    int stem = 0, s = 0, r = 0;
    std::string x_str, dx_str;
    std::string cw;
    std::string diagram_name = "default";
    std::string mode = "add";

    myio::CmdArg1d args = {{"cw", &cw}, {"stem", &stem}, {"s", &s}, {"r", &r}, {"x", &x_str}, {"dx", &dx_str}};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}, {"mode:add/deduce/try", &mode}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

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
        count += diagram.SetRingDiffGlobal(iRing, deg_x, x, dx, r, false);
    }
    else {
        size_t iMod = (size_t)diagram.GetModuleIndexByName(cw);
        count += diagram.SetModuleDiffGlobal(iMod, deg_x, x, dx, r, false);
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
