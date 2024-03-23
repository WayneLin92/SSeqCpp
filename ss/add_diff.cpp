#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"
#include <regex>

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
    std::string diagram_name;
    std::string mode = "add";

    myio::CmdArg1d args = {{"cw", &cw}, {"stem", &stem}, {"s", &s}, {"r", &r}, {"x", &x_str}, {"dx", &dx_str}, {"diagram", &diagram_name}};
    myio::CmdArg1d op_args = {{"mode:add/deduce/try", &mode}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::no_op;

    AdamsDeg deg_x(s, stem + s);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    Diagram diagram(diagram_name, flag);

    /* #Check if x, dx are valid */
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
        if (!ut::has(ss, deg_dx)) {
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
        count += diagram.SetRingDiffGlobal(iRing, deg_x, x, dx, r, false, flag);
    }
    else {
        size_t iMod = (size_t)diagram.GetModuleIndexByName(cw);
        count += diagram.SetModuleDiffGlobal(iMod, deg_x, x, dx, r, false, flag);
    }
    if (count > 0 && (mode == "try" || mode == "deduce")) {
        Logger::LogSummary("Changed differentials", count);
        diagram.DeduceDiffs(0, 500, 0, SSFlag::no_op);  ////
    }
    if (mode == "add" || mode == "deduce") {
        Logger::LogSummary("Changed differentials", count);
        diagram.save(diagram_name, SSFlag::no_op);
    }
    return 0;
}

int main_add_diff_from_file(int argc, char** argv, int& index, const char* desc)
{
    int lineNum = 0;
    std::string filenameLog;
    std::string diagram_name = "mix-hopf";

    myio::CmdArg1d args = {{"filenameLog", &filenameLog}, {"lineNum", &lineNum}};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::no_op;

    myio::AssertFileExists(filenameLog);
    std::ifstream fileLog(filenameLog);
    std::string line;
    int count_lines = 0, count_diffs = 0;
    std::regex is_duduce_regex("^(?:deduce - |)(\\w+) \\((\\d+), (\\d+)\\) d_(\\d+)\\[((?:\\d|\\s|,)*)\\]=\\[((?:\\d|\\s|,)*)\\]"); /* match example: deduce - S0 (66, 6) d_5[0]=[] */
    std::smatch match;

    std::string cw;
    int stem = -1, s = -1, r = 0;
    int1d x, dx;
    Diagram diagram(diagram_name, flag);

    try {
        while (std::getline(fileLog, line) && count_lines++ < lineNum) {
            if (std::regex_search(line, match, is_duduce_regex); match[0].matched) {
                cw = match[1].str();
                stem = std::stoi(match[2].str());
                s = std::stoi(match[3].str());
                r = std::stoi(match[4].str());

                AdamsDeg deg_x(s, stem + s);
                x = myio::Deserialize<int1d>(match[5].str());
                dx = myio::Deserialize<int1d>(match[6].str());

                bool isRing = diagram.GetRingIndexByName(cw) != -1;
                if (isRing) {
                    size_t iRing = (size_t)diagram.GetRingIndexByName(cw);
                    Logger::LogDiff(0, EnumReason::manual, diagram.GetRings()[iRing].name, deg_x, x, dx, r);
                    count_diffs += diagram.SetRingDiffGlobal(iRing, deg_x, x, dx, r, false, flag);
                }
                else {
                    size_t iMod = (size_t)diagram.GetModuleIndexByName(cw);
                    MyException::Assert(iMod != -1, "iMod != -1");
                    Logger::LogDiff(0, EnumReason::manual, diagram.GetModules()[iMod].name, deg_x, x, dx, r);
                    count_diffs += diagram.SetModuleDiffGlobal(iMod, deg_x, x, dx, r, false, flag);
                }
            }
            if (count_lines % 10000 == 0) {
                fmt::print("Deduce trivial diffs\n\n");
                diagram.DeduceTrivialDiffs(flag);
            }
        }
    }
    catch (SSException&) {
    }

    diagram.save(diagram_name, flag);
    return 0;
}

int main_add_diff_from_log(int argc, char** argv, int& index, const char* desc)
{
    int lineNum = 0;
    std::string filenameLog;
    std::string diagram_name = "mix-hopf";

    myio::CmdArg1d args = {{"filenameLog", &filenameLog}, {"lineNum", &lineNum}};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::no_op;

    myio::AssertFileExists(filenameLog);
    DbLog dbLog(filenameLog);
    int count_lines = 0, count_diffs = 0;

    std::string cw;
    int stem = -1, s = -1, r = 0;
    int1d x, dx;
    Diagram diagram(diagram_name, flag);

    try {
        myio::Statement stmt(dbLog, fmt::format("SELECT name, stem, s, r, x, dx FROM log WHERE id<={} AND depth=0 AND reason!=\"Error\" AND reason!=\"nat\"", lineNum));
        while (stmt.step() == MYSQLITE_ROW && count_lines++ < lineNum) {
            cw = stmt.column_str(0);
            stem = stmt.column_int(1);
            s = stmt.column_int(2);
            r = stmt.column_int(3);

            AdamsDeg deg_x(s, stem + s);
            x = myio::Deserialize<int1d>(stmt.column_str(4));
            dx = myio::Deserialize<int1d>(stmt.column_str(5));

            bool isRing = diagram.GetRingIndexByName(cw) != -1;
            if (isRing) {
                size_t iRing = (size_t)diagram.GetRingIndexByName(cw);
                Logger::LogDiff(0, EnumReason::manual, diagram.GetRings()[iRing].name, deg_x, x, dx, r);
                count_diffs += diagram.SetRingDiffGlobal(iRing, deg_x, x, dx, r, false, flag);
            }
            else {
                size_t iMod = (size_t)diagram.GetModuleIndexByName(cw);
                MyException::Assert(iMod != -1, "iMod != -1");
                Logger::LogDiff(0, EnumReason::manual, diagram.GetModules()[iMod].name, deg_x, x, dx, r);
                count_diffs += diagram.SetModuleDiffGlobal(iMod, deg_x, x, dx, r, false, flag);
            }
        }
        if (count_lines % 10000 == 0) {
            fmt::print("Deduce trivial diffs\n\n");
            diagram.DeduceTrivialDiffs(flag);
        }
    }
    catch (SSException&) {
    }

    diagram.save(diagram_name, flag);
    return 0;
}

int main_add_cofseq_diff(int argc, char** argv, int& index, const char* desc)
{
    int stem = 0, s = 0, r = 0;
    std::string x_str, dx_str;
    std::string cofseq_name;
    int iCs = 0;
    std::string diagram_name = "mix-exact";

    myio::CmdArg1d args = {{"cofseq", &cofseq_name}, {"index", &iCs}, {"stem", &stem}, {"s", &s}, {"r", &r}, {"x", &x_str}, {"dx", &dx_str}};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    AdamsDeg deg_x(s, stem + s);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    auto flag = SSFlag::cofseq;

    Diagram diagram(diagram_name, flag);

    /* #Check if x, dx are valid */
    int index_cofseq = diagram.GetCofSeqIndexByName(cofseq_name);
    auto& cofseq = diagram.GetCofSeqs()[index_cofseq];  //// const

    std::sort(x.begin(), x.end());
    std::sort(dx.begin(), dx.end());
    auto& ss1 = cofseq.nodes_ss[iCs]->front();
    auto& ss2 = cofseq.nodes_ss[(iCs + 1) % 3]->front();
    if (!x.empty()) {
        if (!ut::has(ss1, deg_x)) {
            fmt::print("deg_x not found\n");
            return -101;
        }
        if (x.front() < 0 || ss1.at(deg_x).levels.size() <= (size_t)x.back()) {
            fmt::print("Invalid x\n");
            return -102;
        }
    }
    if (!dx.empty()) {
        int stem_map = cofseq.degMap[iCs].stem();
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, stem_map + r);
        if (!ut::has(ss2, deg_dx)) {
            fmt::print("deg_dx not found\n");
            return -103;
        }
        if (dx.front() < 0 || ss2.at(deg_dx).levels.size() <= (size_t)dx.back()) {
            fmt::print("Invalid dx\n");
            return -104;
        }
    }

    /* #Add diff */
    int count = diagram.SetDiffGlobalCofseq(cofseq, iCs, deg_x, x, dx, r, false, flag);
    diagram.save(diagram_name, SSFlag::cofseq);
    Logger::LogSummary("Changed differentials", count);
    return 0;
}
