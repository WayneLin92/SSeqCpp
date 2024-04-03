#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "groebner_res.h"
#include "main.h"

#include <cstring>
#include <regex>
#include <sstream>

void ResolveV2(const Mod1d& rels, const int1d& v_degs, int t_trunc, int stem_trunc, const std::string& db_filename, const std::string& tablename)
{
    using namespace steenrod;

    auto gb = AdamsRes::load(db_filename, tablename, t_trunc, stem_trunc);
    Resolve(gb, rels, v_degs, t_trunc, stem_trunc, db_filename, tablename);

    std::cout << "gb.dim_Ext()=" << gb.dim_Ext() << '\n';
    std::cout << "gb.dim_Gb()=" << gb.dim_Gb() << '\n';
}

int GetCoh(int1d& v_degs, Mod1d& rels, int t_max, const std::string& name);

int main_res(int argc, char** argv, int& index, const char* desc)
{
    std::string cw = "S0";
    int t_max = 100, stem_max = DEG_MAX;

    myio::CmdArg1d args = {{"cw", &cw}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {{"stem_max", &stem_max}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
        
/* Prevent double run on linux */
#ifdef __linux__
    if (IsAdamsRunning(fmt::format("./Adams res {} ", cw))) {
        fmt::print("Error: ./Adams res {} is already running.\n", cw);
        return -1;
    }
#endif

    if (stem_max > DEG_MAX) {
        stem_max = DEG_MAX;
        fmt::print("stem_max is truncated to {}.", stem_max);
    }

    int1d v_degs;
    Mod1d rels;
    int d_max = std::min(t_max, stem_max);
    if (int error = GetCoh(v_degs, rels, d_max, cw))
        return -2;

    std::string db_filename = cw + "_Adams_res.db";
    std::string tablename = cw + "_Adams_res";
    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}