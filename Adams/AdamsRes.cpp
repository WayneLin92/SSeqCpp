#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "groebner_steenrod.h"

#include <cstring>
#include <sstream>

void ResolveV2(const std::string& X, int t_trunc, int stem_trunc, std::string& db_filename, const std::string& tablename)
{
    using namespace steenrod;

    Mod1d rels;
    if (X == "S0") {
        for (int i = 0; 1; ++i) {
            rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
            if (MMilnor::P(i, i + 1).deg() > t_trunc)
                break;
        }
    }
    else {
        std::cout << "X=" << X << " is not supported yet\n";
        return;
    }

    auto gb = AdamsRes::load(db_filename, tablename, t_trunc, stem_trunc);

    Resolve(gb, rels, t_trunc, stem_trunc, db_filename);

    std::cout << "gb.dim_Ext()=" << gb.dim_Ext() << '\n';
    std::cout << "gb.dim_Gb()=" << gb.dim_Gb() << '\n';
}

int main_res(int argc, char** argv, int index)
{
    std::string X = "S0";
    int t_max = 392, stem_max = 261;
    std::string db_filename = "S0_Adams_res.db";
    std::string tablename = "S0_Adams_res";
    int num_threads = 128;

    if (argc >= 2 && strcmp(argv[1], "-h") == 0) {
        std::cout << "Calculate the minimal resolution for the Adams spectral sequence\n";
        std::cout << "Usage:\n  Adams res [X] [t_max] [stem_max] [db_filename] [tablename] [num_threads]\n\n";

        std::cout << "Default values:\n";
        std::cout << "  X = " << X << "\n";
        std::cout << "  t_max = " << t_max << "\n";
        std::cout << "  stem_max = " << stem_max << "\n";
        std::cout << "  db_filename =" << db_filename << "\n";
        std::cout << "  tablename =" << tablename << "\n";
        std::cout << "  num_threads = " << num_threads << "\n\n";

        std::cout << "Version:\n  2.1 (2022-08-07)" << std::endl;
        return 0;
    }
    if (myio::load_op_arg(argc, argv, ++index, "X", X))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "tablename", tablename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "num_threads", num_threads))
        return index;

    ut::FUTURE_NUM_THREADS = num_threads;
    bench::Timer timer;

    ResolveV2(X, t_max, stem_max, db_filename, tablename);

#ifndef MYDEPLOY
    bench::Counter::print();
#endif
    return 0;
}