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
void Coh_j(int1d& v_degs, Mod1d& rels, int t_max);
void Coh_X2(int1d& v_degs, Mod1d& rels, int t_max);

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
    else if (cw == "j")
        Coh_j(v_degs, rels, d_max);
    else {
        fmt::print("Unsupported arugment cw:{}\n", cw);
        return 0;
    }

    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void Coh_RP(int1d& v_degs, Mod1d& rels, int n1, int n2, int t_max);

int main_res_RP(int argc, char** argv, int index)
{
    int t_max = 100, stem_max = DEG_MAX;
    int n1 = 0, n2 = 0;

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Calculate the minimal resolution of P_{n1}^{n2} for the Adams spectral sequence");
        fmt::print("Usage:\n  Adams res_RP <n1> <n2> <t_max> [stem_max]\n\n");

        fmt::print("Default values:\n");
        fmt::print("  stem_max = {}\n", stem_max);

        std::cout << VERSION << std::endl;
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
    if (n1 < 1 || n1 > n2) {
        fmt::print("Inappropriate n1, n2 argument");
        return -(++index);
    }
    std::string db_filename = fmt::format("RP{}_{}_Adams_res.db", n1, n2);
    std::string tablename = fmt::format("RP{}_{}_Adams_res", n1, n2);

    int1d v_degs;
    Mod1d rels;
    Coh_RP(v_degs, rels, n1, n2, t_max);
    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}