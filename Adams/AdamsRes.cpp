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
void Coh_Chopf(int1d& v_degs, Mod1d& rels, int n, int t_max);
void Coh_C2h4(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_j(int1d& v_degs, Mod1d& rels, int t_max);

int main_res(int argc, char** argv, int index)
{
    std::string cw = "S0";
    int t_max = 100, stem_max = DEG_MAX;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Calculate the minimal resolution for the Adams spectral sequence\n");
        fmt::print("Usage:\n  Adams res <cw:S0/tmf/X2> <t_max> [stem_max]\n\n");

        fmt::print("Default values:\n");
        fmt::print("  stem_max = {}\n", stem_max);

        fmt::print("{}\n", VERSION);
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cw", cw))
        return -index;
    if (myio::load_arg(argc, argv, ++index, "t_max", t_max))
        return -index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return -index;

    if (stem_max > DEG_MAX) {
        stem_max = DEG_MAX;
        fmt::print("stem_max is truncated to {}.", stem_max);
    }
    std::string db_filename = cw + "_Adams_res.db";
    std::string tablename = cw + "_Adams_res";

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
        Coh_Chopf(v_degs, rels, 0, d_max);
    else if (cw == "Ceta")
        Coh_Chopf(v_degs, rels, 1, d_max);
    else if (cw == "Cnu")
        Coh_Chopf(v_degs, rels, 2, d_max);
    else if (cw == "Csigma")
        Coh_Chopf(v_degs, rels, 3, d_max);
    else if (cw == "C2h4")
        Coh_C2h4(v_degs, rels, d_max);
    else if (cw == "j")
        Coh_j(v_degs, rels, d_max);
    else {
        fmt::print("Unsupported arugment cw:{}\n", cw);
        return 0;
    }

    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void normalize_RP(int n1, int n2, int& n1_, int& n2_);
void Coh_RP(int1d& v_degs, Mod1d& rels, Mod1d& convert_v, int n1, int n2, int t_max);

int main_res_RP(int argc, char** argv, int index)
{
    int t_max = 100, stem_max = DEG_MAX;
    int n1 = 0, n2 = 0;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Calculate the minimal resolution of P_{{n1}}^{{n2}} for the Adams spectral sequence\n");
        fmt::print("Usage:\n  Adams res_RP <n1> <n2> <t_max> [stem_max]\n\n");

        fmt::print("Default values:\n");
        fmt::print("  stem_max = {}\n", stem_max);

        fmt::print("{}\n", VERSION);
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "n1", n1))
        return -index;
    if (myio::load_arg(argc, argv, ++index, "n2", n2))
        return -index;
    if (myio::load_arg(argc, argv, ++index, "t_max", t_max))
        return -index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return -index;
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
    std::string db_filename = fmt::format("RP{}_{}_Adams_res.db", n1_, n2_);
    std::string tablename = fmt::format("RP{}_{}_Adams_res", n1_, n2_);

    int1d v_degs;
    Mod1d rels, convert_v;
    Coh_RP(v_degs, rels, convert_v, n1_, n2_, t_max);
    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void SetCohMapImages(std::string& cw1, std::string& cw2, Mod1d& images, int& sus);
void SetCohMap(const std::string& db_map, const std::string& table_map, const Mod1d& images, int sus);

int main_map_coh(int argc, char** argv, int index)
{
    std::string cw1, cw2;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Set the map H^*(cw2) -> H^*(cw1) \n");
        fmt::print("Usage:\n  Adams map_coh <cw1> <cw2>\n\n");

        fmt::print("{}\n", VERSION);
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cw1", cw1))
        return -index;
    if (myio::load_arg(argc, argv, ++index, "cw2", cw2))
        return -index;

    Mod1d images;
    int sus = 0;
    SetCohMapImages(cw1, cw2, images, sus);
    std::string db_filename = fmt::format("map_Adams_res_{}_to_{}.db", cw1, cw2);
    std::string tablename = fmt::format("map_Adams_res_{}_to_{}", cw1, cw2);
    SetCohMap(db_filename, tablename, images, sus);

    return 0;
}