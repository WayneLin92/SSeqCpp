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

void Coh_S0(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_tmf(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_ko(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_X2(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_Chn(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_tmf_Chn(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_C2hn(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_DC2hn(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_three_cell(int1d& v_degs, Mod1d& rels, int n1, int n2, int t_max);
void Coh_smash_2cell(int1d& v_degs, Mod1d& rels, int n1, int n2, int t_max);
void Coh_j(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_j_C2(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_Fphi(int1d& v_degs, Mod1d& rels, Mod1d& cell_reduced, int1d& min_rels, int t_max);
int res_P(const std::string& cw, int t_max, int stem_max);

int main_res(int argc, char** argv, int& index, const char* desc)
{
    std::string cw = "S0";
    int t_max = 100, stem_max = DEG_MAX;

    myio::CmdArg1d args = {{"cw", &cw}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {{"stem_max", &stem_max}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    if (stem_max > DEG_MAX) {
        stem_max = DEG_MAX;
        fmt::print("stem_max is truncated to {}.", stem_max);
    }

    std::regex is_P_regex("^(?:tmf_|X2_|)(?:R|C|H)P(?:m|)\\d+_(?:m|)\\d+$"); /* match example: RP1_4, CPm1_10 */
    std::smatch match;
    if (std::regex_search(cw, match, is_P_regex); match[0].matched)
        return res_P(cw, t_max, stem_max);

    int d_max = std::min(t_max, stem_max);
    int1d v_degs;
    Mod1d rels;
    if (cw == "S0")
        Coh_S0(v_degs, rels, d_max);
    else if (cw == "X2")
        Coh_X2(v_degs, rels, d_max);
    else if (cw == "tmf")
        Coh_tmf(v_degs, rels, d_max);
    else if (cw == "ko")
        Coh_ko(v_degs, rels, d_max);
    else if (cw == "C2")
        Coh_Chn(v_degs, rels, 0, d_max);
    else if (cw == "Ceta")
        Coh_Chn(v_degs, rels, 1, d_max);
    else if (cw == "Cnu")
        Coh_Chn(v_degs, rels, 2, d_max);
    else if (cw == "Csigma")
        Coh_Chn(v_degs, rels, 3, d_max);
    else if (cw == "C2h4")
        Coh_C2hn(v_degs, rels, 4, d_max);
    else if (cw == "C2h5")
        Coh_C2hn(v_degs, rels, 5, d_max);
    else if (cw == "C2h6")
        Coh_C2hn(v_degs, rels, 6, d_max);
    else if (cw == "DC2h4")
        Coh_DC2hn(v_degs, rels, 4, d_max);
    else if (cw == "DC2h5")
        Coh_DC2hn(v_degs, rels, 5, d_max);
    else if (cw == "DC2h6")
        Coh_DC2hn(v_degs, rels, 6, d_max);
    else if (cw == "CW_2_eta")
        Coh_three_cell(v_degs, rels, 0, 1, d_max);
    else if (cw == "CW_eta_nu")
        Coh_three_cell(v_degs, rels, 1, 2, d_max);
    else if (cw == "CW_nu_sigma")
        Coh_three_cell(v_degs, rels, 2, 3, d_max);
    else if (cw == "CW_eta_2")
        Coh_three_cell(v_degs, rels, 1, 0, d_max);
    else if (cw == "CW_nu_eta")
        Coh_three_cell(v_degs, rels, 2, 1, d_max);
    else if (cw == "CW_sigma_nu")
        Coh_three_cell(v_degs, rels, 3, 2, d_max);
    else if (cw == "C2_Ceta")
        Coh_smash_2cell(v_degs, rels, 0, 1, d_max);
    else if (cw == "Ceta_nu")
        Coh_smash_2cell(v_degs, rels, 1, 2, d_max);
    else if (cw == "Cnu_Csigam")
        Coh_smash_2cell(v_degs, rels, 2, 3, d_max);
    else if (cw == "j")
        Coh_j(v_degs, rels, d_max);
    else if (cw == "j_C2")
        Coh_j_C2(v_degs, rels, d_max);
    else if (cw == "tmf_C2")
        Coh_tmf_Chn(v_degs, rels, 0, d_max);
    else if (cw == "tmf_Ceta")
        Coh_tmf_Chn(v_degs, rels, 1, d_max);
    else if (cw == "tmf_Cnu")
        Coh_tmf_Chn(v_degs, rels, 2, d_max);
    else if (cw == "Fphi") {
        Mod1d tmp;
        int1d tmp_int1d;
        Coh_Fphi(v_degs, rels, tmp, tmp_int1d, d_max);
    }
    else {
        fmt::print("Unsupported arugment cw={}\n", cw);
        return -1;
    }

    std::string db_filename = cw + "_Adams_res.db";
    std::string tablename = cw + "_Adams_res";
    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void ParseFP(const std::string& cw, int& hopf, int& n1, int& n2);
void normalize_RP(int n1, int n2, int& n1_, int& n2_);
void Coh_P(int1d& v_degs, Mod1d& rels, Mod1d& cell_reduced, int1d& min_rels, int n1, int n2, int t_max, const std::string& over, int hopf);
void Coh_X2_RP(int1d& v_degs, Mod1d& rels, Mod1d& cell_reduced, int1d& min_rels, int n1, int n2, int t_max);

int res_P(const std::string& cw, int t_max, int stem_max)
{
    std::regex is_P_regex("^(tmf_|X2_|)((R|C|H)P(?:m|)\\d+_(?:m|)\\d+)$"); /* match example: RP1_4, CPm1_10 */
    std::smatch match;
    std::regex_search(cw, match, is_P_regex);

    std::string ring = match[1].str();
    int hopf, n1, n2;
    ParseFP(match[2].str(), hopf, n1, n2);
    std::string field = match[3].str();

    if (!(n1 + 1 < n2)) {
        fmt::print("We need n1 + 1 < n2");
        return -1;
    }
    int n1_ = n1, n2_ = n2;

    if (field == "R") {
        normalize_RP(n1, n2, n1_, n2_);
        if (n1 != n1_)
            fmt::print("We use n1={} n2={} because of James periodicity\n", n1_, n2_);
    }

    int1d v_degs;
    Mod1d rels, tmp_Mod1d;
    int1d tmp_ind1d;
    if (ring == "")
        Coh_P(v_degs, rels, tmp_Mod1d, tmp_ind1d, n1_, n2_, t_max, "S0", hopf);
    else if (ring == "tmf_")
        Coh_P(v_degs, rels, tmp_Mod1d, tmp_ind1d, n1_, n2_, t_max, "tmf", hopf);
    else if (ring == "X2_")
        Coh_X2_RP(v_degs, rels, tmp_Mod1d, tmp_ind1d, n1_, n2_, t_max);  ////
    else {
        fmt::print("ring={} not supported\n", ring);
        return 1;
    }
    std::string db_filename, tablename;
    std::string str_n1 = std::to_string(abs(n1_));
    std::string str_n2 = std::to_string(abs(n2_));
    if (n1_ < 0)
        str_n1 = "m" + str_n1;
    if (n2_ < 0)
        str_n2 = "m" + str_n2;
    db_filename = fmt::format("{}{}P{}_{}_Adams_res.db", ring, field, str_n1, str_n2);
    tablename = fmt::format("{}{}P{}_{}_Adams_res", ring, field, str_n1, str_n2);

    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void SetCohMap(const std::string& cw1, const std::string& cw2, std::string& from, std::string& to, Mod1d& images, int& sus, int& fil);
void SetDbCohMap(const std::string& db_map, const std::string& table_map, const std::string& from, const std::string& to, const Mod1d& images, int sus, int fil);

int main_map_coh(int argc, char** argv, int& index, const char* desc)
{
    std::string cw1, cw2;

    myio::CmdArg1d args = {{"cw1", &cw1}, {"cw2", &cw2}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Mod1d images;
    int sus = 0, fil = 0;
    std::string from, to;
    SetCohMap(cw1, cw2, from, to, images, sus, fil);
    std::string db_filename = fmt::format("map_Adams_res_{}_to_{}.db", cw1, cw2);
    std::string tablename = fmt::format("map_Adams_res_{}_to_{}", cw1, cw2);
    SetDbCohMap(db_filename, tablename, from, to, images, sus, fil);

    return 0;
}