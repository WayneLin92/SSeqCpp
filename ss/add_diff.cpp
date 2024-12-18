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
    std::string cat_name;
    std::string mode = "add";

    myio::CmdArg1d args = {{"cw", &cw}, {"stem", &stem}, {"s", &s}, {"r", &r}, {"x", &x_str}, {"dx", &dx_str}, {"category", &cat_name}};
    myio::CmdArg1d op_args = {{"mode:add/deduce/try", &mode}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::no_op;

    AdamsDeg deg_x(s, stem + s);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    Category category(cat_name, "", flag);

    /* #Check if x, dx are valid */
    auto iCw = category.GetIndexCwByName(cw);
    auto& ss = category.GetNodesSS(iCw).front();
    if (!x.empty()) {
        if (!ss.has(deg_x)) {
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
        if (!ss.has(deg_dx)) {
            fmt::print("deg_dx not found\n");
            return -103;
        }
        if (dx.front() < 0 || ss.at(deg_dx).levels.size() <= (size_t)dx.back()) {
            fmt::print("Invalid dx\n");
            return -104;
        }
    }

    /* #Add diff */
    SSRet rt = category.SetCwDiffGlobal(iCw, deg_x, x, dx, r, false, flag);
    if (rt)
        throw ErrorIdMsg(rt.code, "Failed to SetCwDiffGlobal()");
    if (rt.IsChanged() && (mode == "try" || mode == "deduce")) {
        if (rt += category.DeduceDiffs(0, 500, 1, 0, SSFlag::no_op))
            throw ErrorIdMsg(rt.code, "Failed to DeduceDiffs()");
    }
    if (mode == "add" || mode == "deduce")
        category.SaveNodes(cat_name, "", true, SSFlag::no_op);
    category.PrintSummary();
    return 0;
}

int main_add_diff_from_file(int argc, char** argv, int& index, const char* desc)
{
    int lineNum = 0;
    std::string filenameLog;
    std::string cat_name;

    myio::CmdArg1d args = {{"filenameLog", &filenameLog}, {"lineNum", &lineNum}, {"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::no_op;

    myio::AssertFileExists(filenameLog);
    std::ifstream fileLog(filenameLog);
    std::string line;
    int count_lines = 0, count_diffs = 0;
    std::regex regex_deduce("^(?:deduce - |)((?:\\w|_|)+) \\((\\d+),(?:\\s|)(\\d+)\\) d_(\\d+)\\[((?:\\d|\\s|,)*)\\]=\\[((?:\\d|\\s|,)*)\\]"); /* match example: deduce - S0 (66, 6) d_5[0]=[] */
    std::smatch match;

    std::string cw;
    int stem = -1, s = -1, r = 0;
    int1d x, dx;
    Category category(cat_name, "", flag);

    while (std::getline(fileLog, line) && count_lines++ < lineNum) {
        if (std::regex_search(line, match, regex_deduce); match[0].matched) {
            cw = match[1].str();
            stem = std::stoi(match[2].str());
            s = std::stoi(match[3].str());
            r = std::stoi(match[4].str());

            AdamsDeg deg_x(s, stem + s);
            x = myio::Deserialize<int1d>(match[5].str());
            dx = myio::Deserialize<int1d>(match[6].str());

            auto iCw = category.GetIndexCwByName(cw);
            if (category.IsNewDiff(category.GetNodesSS(iCw), deg_x, x, dx, r)) {
                Logger::LogDiff(0, EnumReason::manual, category.GetCwName(iCw), deg_x, r, x, dx, "", flag);
                if (auto rt = category.SetCwDiffGlobal(iCw, deg_x, x, dx, r, false, flag))
                    throw ErrorIdMsg(rt.code, "Failed to SetCwDiffGlobal()");
                else if (rt.IsChanged())
                    ++count_diffs;
            }
        }
        if (count_diffs % 10000 == 0) {
            fmt::print("Deduce trivial diffs\n\n");
            if (auto rt = category.DeduceTrivialDiffs(flag))
                throw ErrorIdMsg(rt.code, "Failed to deduce trivial diffs");
        }
    }

    category.SaveNodes(cat_name, "", true, flag);
    return 0;
}

int main_add_diff_from_log(int argc, char** argv, int& index, const char* desc)
{
    std::string cat_name;
    int id_max = 0;
    std::string filenameLog;

    myio::CmdArg1d args = {{"category", &cat_name}, {"filenameLog", &filenameLog}, {"id_max", &id_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    auto flag = SSFlag::cofseq | SSFlag::no_exclusions;

    myio::AssertFileExists(filenameLog);
    DbLog dbLog(filenameLog);
    Category category(cat_name, "", flag);

    try {
        myio::Statement stmt(dbLog, fmt::format("SELECT id, name, s, t, r, x, dx, name like '%' || ':' || '%', reason like '%' || 'I' FROM log WHERE id<={} AND depth=0 AND reason!=\"E\" AND reason!=\"N\"", id_max));
        while (stmt.step() == MYSQLITE_ROW) {
            auto [id, name, deg, r, x, dx, isCs, isDInv] = columns_diff(stmt);
            fmt::print("id={} ", id);
            if (!isDInv)
                Logger::LogDiff(0, EnumReason::manual, name, deg, r, x, dx, "", flag);
            else
                Logger::LogDiffInv(0, EnumReason::manual, name, deg, r, x, dx, "", flag);
            category.SetUnivDiffGlobal(name, deg, r, x, dx, isCs, isDInv, flag);
        }
    }
    catch (NoException&) {
    }

    category.SaveNodes(cat_name, "", true, flag);
    return 0;
}

int main_add_cofseq_diff(int argc, char** argv, int& index, const char* desc)
{
    int stem = 0, s = 0, r = 0;
    std::string x_str, dx_str;
    std::string cofseq_name;
    int iTri = 0;
    std::string cat_name = "mix-exact";

    myio::CmdArg1d args = {{"cofseq", &cofseq_name}, {"index", &iTri}, {"stem", &stem}, {"s", &s}, {"r", &r}, {"x", &x_str}, {"dx", &dx_str}};
    myio::CmdArg1d op_args = {{"category", &cat_name}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    AdamsDeg deg_x(s, stem + s);
    int1d x = myio::Deserialize<int1d>(x_str);
    int1d dx = myio::Deserialize<int1d>(dx_str);

    auto flag = SSFlag::cofseq;

    Category category(cat_name, "", flag);

    /* #Check if x, dx are valid */
    int index_cofseq = category.GetCofSeqIndexByName(cofseq_name);
    auto& cofseq = category.GetCofSeqs()[index_cofseq];  //// const

    std::sort(x.begin(), x.end());
    std::sort(dx.begin(), dx.end());
    auto& ss1 = cofseq.nodes_ss[iTri]->front();
    auto& ss2 = cofseq.nodes_ss[NextiTri(iTri)]->front();
    if (!x.empty()) {
        if (!ss1.has(deg_x)) {
            fmt::print("deg_x not found\n");
            return -101;
        }
        if (x.front() < 0 || ss1.at(deg_x).levels.size() <= (size_t)x.back()) {
            fmt::print("Invalid x\n");
            return -102;
        }
    }
    if (!dx.empty()) {
        int stem_map = cofseq.degMap[iTri].stem();
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, stem_map + r);
        if (!ss2.has(deg_dx)) {
            fmt::print("deg_dx not found\n");
            return -103;
        }
        if (dx.front() < 0 || ss2.at(deg_dx).levels.size() <= (size_t)dx.back()) {
            fmt::print("Invalid dx\n");
            return -104;
        }
    }

    /* #Add diff */
    if (auto rt = category.SetDiffGlobalCofseq(cofseq, iTri, deg_x, x, dx, r, false, flag))
        throw ErrorIdMsg(rt.code, "Failed to SetDiffGlobalCofseq()");
    category.SaveNodes(cat_name, "", true, SSFlag::cofseq);
    category.PrintSummary();
    return 0;
}
