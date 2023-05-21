#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "groebner_res.h"
#include "main.h"

#include <cstring>
#include <sstream>

void ResolveV2(const Mod1d& rels, const int1d& v_degs, int t_trunc, int stem_trunc, const std::string& db_filename, const std::string& tablename)
{
    using namespace steenrod;

    auto gb = AdamsRes::load(db_filename, tablename, t_trunc, stem_trunc);
    Resolve(gb, rels, v_degs, t_trunc, stem_trunc, db_filename, tablename);

    std::cout << "gb.dim_Ext()=" << gb.dim_Ext() << '\n';
    std::cout << "gb.dim_Gb()=" << gb.dim_Gb() << '\n';
}

Mod1d GetS0Rels(int t_trunc)
{
    using namespace steenrod;
    Mod1d rels;
    for (int i = 0; (1 << i) <= t_trunc; ++i)
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
    return rels;
}

void Coh_S0(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_tmf(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_ko(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_X2(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_Chn(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_tmf_Chn(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_C2hn(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_j(int1d& v_degs, Mod1d& rels, int t_max);

int main_res(int argc, char** argv, int& index, const char* desc)
{
    std::string ring = "S0";
    int t_max = 100, stem_max = DEG_MAX;

    myio::CmdArg1d args = {{"ring", &ring}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {{"stem_max", &stem_max}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    if (stem_max > DEG_MAX) {
        stem_max = DEG_MAX;
        fmt::print("stem_max is truncated to {}.", stem_max);
    }
    std::string db_filename = ring + "_Adams_res.db";
    std::string tablename = ring + "_Adams_res";

    int d_max = std::min(t_max, stem_max);
    int1d v_degs;
    Mod1d rels;
    if (ring == "S0")
        Coh_S0(v_degs, rels, d_max);
    else if (ring == "X2")
        Coh_X2(v_degs, rels, d_max);
    else if (ring == "tmf")
        Coh_tmf(v_degs, rels, d_max);
    else if (ring == "ko")
        Coh_ko(v_degs, rels, d_max);
    else if (ring == "C2")
        Coh_Chn(v_degs, rels, 0, d_max);
    else if (ring == "Ceta")
        Coh_Chn(v_degs, rels, 1, d_max);
    else if (ring == "Cnu")
        Coh_Chn(v_degs, rels, 2, d_max);
    else if (ring == "Csigma")
        Coh_Chn(v_degs, rels, 3, d_max);
    else if (ring == "C2h4")
        Coh_C2hn(v_degs, rels, 4, d_max);
    else if (ring == "C2h5")
        Coh_C2hn(v_degs, rels, 5, d_max);
    else if (ring == "C2h6")
        Coh_C2hn(v_degs, rels, 6, d_max);
    else if (ring == "j")
        Coh_j(v_degs, rels, d_max);
    else if (ring == "tmf_C2")
        Coh_tmf_Chn(v_degs, rels, 0, d_max);
    else if (ring == "tmf_Ceta")
        Coh_tmf_Chn(v_degs, rels, 1, d_max);
    else if (ring == "tmf_Cnu")
        Coh_tmf_Chn(v_degs, rels, 2, d_max);
    else {
        fmt::print("Unsupported arugment ring:{}\n", ring);
        return 0;
    }

    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void normalize_RP(int n1, int n2, int& n1_, int& n2_);
void Coh_RP(int1d& v_degs, Mod1d& rels, Mod1d& convert_v, int n1, int n2, int t_max);
void Coh_X2_RP(int1d& v_degs, Mod1d& rels, Mod1d& convert_v, int n1, int n2, int t_max);

int main_res_RP(int argc, char** argv, int& index, const char* desc)
{
    int n1 = 0, n2 = 0;
    int t_max = 100, stem_max = DEG_MAX;
    std::string ring = "S0";

    myio::CmdArg1d args = {{"n1", &n1}, {"n2", &n2}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {{"stem_max", &stem_max}, {"ring", &ring}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    if (stem_max > DEG_MAX) {
        stem_max = DEG_MAX;
        fmt::print("stem_max is truncated to {}.", stem_max);
    }
    if (!(n1 % 2 == 1 && n1 >= 1 && n1 + 1 < n2)) {
        fmt::print("We need n1 % 2 == 1 && n1 >= 1 && n1 + 1 < n2");
        return -index;
    }
    int n1_ = -1, n2_ = -1;
    normalize_RP(n1, n2, n1_, n2_);
    if (n1 != n1_)
        fmt::print("We use n1={} n2={} instead because of James periodicity", n1_, n2_);

    int1d v_degs;
    Mod1d rels, convert_v;
    if (ring == "S0")
        Coh_RP(v_degs, rels, convert_v, n1_, n2_, t_max);
    else if (ring == "X2")
        Coh_X2_RP(v_degs, rels, convert_v, n1_, n2_, t_max);
    else {
        fmt::print("ring={} not supported\n", ring);
        return 1;
    }
    std::string db_filename, tablename;
    if (ring != "S0") {
        db_filename = fmt::format("{}_RP{}_{}_Adams_res.db", ring, n1_, n2_);
        tablename = fmt::format("{}_RP{}_{}_Adams_res", ring, n1_, n2_);
    }
    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void SetCohMapImages(std::string& cw1, std::string& cw2, Mod1d& images, int& sus);
void SetCohMap(const std::string& db_map, const std::string& table_map, const Mod1d& images, int sus);

int main_map_coh(int argc, char** argv, int& index, const char* desc)
{
    std::string cw1, cw2;

    myio::CmdArg1d args = {{"cw1", &cw1}, {"cw2", &cw2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Mod1d images;
    int sus = 0;
    SetCohMapImages(cw1, cw2, images, sus);
    std::string db_filename = fmt::format("map_Adams_res_{}_to_{}.db", cw1, cw2);
    std::string tablename = fmt::format("map_Adams_res_{}_to_{}", cw1, cw2);
    SetCohMap(db_filename, tablename, images, sus);

    return 0;
}