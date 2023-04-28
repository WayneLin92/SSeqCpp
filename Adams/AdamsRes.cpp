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
    int t_max = 100, stem_max = 261;
    std::string db_filename = "<cw>_Adams_res.db";
    std::string tablename = "<cw>_Adams_res";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        fmt::print("Calculate the minimal resolution for the Adams spectral sequence\n");
        fmt::print("Usage:\n  Adams res <cw:S0/tmf/X2> [t_max] [stem_max] [db_filename] [tablename]\n\n");

        fmt::print("Default values:\n");
        fmt::print("  t_max = {}\n", t_max);
        fmt::print("  stem_max = {}\n", stem_max);
        fmt::print("  db_filename = {}\n", db_filename);
        fmt::print("  tablename = {}\n\n", tablename);

        fmt::print("{}\n", VERSION);
        return 0;
    }
    if (myio::load_arg(argc, argv, ++index, "cw", cw))
        return -index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return -index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return -index;
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return -index;
    if (myio::load_op_arg(argc, argv, ++index, "tablename", tablename))
        return -index;

    if (db_filename == "<cw>_Adams_res.db")
        db_filename = cw + "_Adams_res.db";
    if (tablename == "<cw>_Adams_res")
        tablename = cw + "_Adams_res";

    int1d v_degs;
    Mod1d rels;
    if (cw == "S0")
        Coh_S0(v_degs, rels, t_max);
    else if (cw == "X2")
        Coh_X2(v_degs, rels, t_max);
    else if (cw == "tmf")
        Coh_tmf(v_degs, rels, t_max);
    else if (cw == "j")
        Coh_j(v_degs, rels, t_max);
    else {
        fmt::print("Unsupported arugment cw:{}\n", cw);
        return 0;
    }

    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);
    return 0;
}

void Coh_RPn(int1d& v_degs, Mod1d& rels, int n, int t_max);

int main_res_RP(int argc, char** argv, int index)
{
    int t_max = 100, stem_max = 261;
    int n = -1;
    std::string db_filename = "RPn_Adams_res.db";
    std::string tablename = "RPn_Adams_res";

    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "-h") == 0) {
        std::cout << "Calculate the minimal resolution of RP^n for the Adams spectral sequence\n";
        std::cout << "Usage:\n  Adams res_RP [n] [t_max] [stem_max] [db_filename] [tablename]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  n = inf\n";
        std::cout << "  t_max = " << t_max << "\n";
        std::cout << "  stem_max = " << stem_max << "\n";
        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  tablename = " << tablename << "\n\n";

        std::cout << VERSION << std::endl;
        return 0;
    }
    if (argc > index + 1 && strcmp(argv[size_t(index + 1)], "inf") == 0) {
        n = -1;
        ++index;
    }
    else if (myio::load_op_arg(argc, argv, ++index, "n", n))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "tablename", tablename))
        return index;
    std::string strN = n == -1 ? "inf" : std::to_string(n);
    if (n == -1)
        n = t_max;
    if (db_filename == "RPn_Adams_res.db")
        db_filename = "RP" + strN + "_Adams_res.db";
    if (tablename == "RPn_Adams_res")
        tablename = "RP" + strN + "_Adams_res";

    int1d v_degs;
    Mod1d rels;
    Coh_RPn(v_degs, rels, n, t_max);
    ResolveV2(rels, v_degs, t_max, stem_max, db_filename, tablename);

#ifndef MYDEPLOY
    bench::Counter::print();
#endif
    return 0;
}