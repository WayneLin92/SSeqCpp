#include "algebras/benchmark.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "groebner_steenrod.h"

#include <cstring>
#include <sstream>

#ifndef MYDEPLOY
std::vector<int> bench::Counter::counts_ = {0, 0, 0};
#endif

void AdamsE2(const std::string& X, int t_trunc, int stem_trunc, std::string& db_filename)
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

#ifdef MYDEPLOY
    // ResetDb(db_filename);
    auto gb = SteenrodMRes::load(db_filename, t_trunc, stem_trunc);
#else
    auto gb = SteenrodMRes(t_trunc, {}, {}, {}, {}, 0, 0);
#endif

    ResolveMRes(gb, rels, t_trunc, stem_trunc, db_filename);

    std::cout << "gb.dim_Ext()=" << gb.dim_Ext() << '\n';
    std::cout << "gb.dim_Gb()=" << gb.dim_Gb() << '\n';
}

int main(int argc, char** argv)
{
    std::string X = "S0";
    int t_max = 392, stem_max = 261;
    std::string db_filename = "S0_Adams_res.db";
    auto num_threads = 128;

    if (argc >= 2 && strcmp(argv[1], "-h") == 0) {
        std::cout << "Calculate the minimal resolution for the Adams spectral sequence\n";
        std::cout << "Usage:\n  AdamsRes <X> [t_max] [stem_max] [db_filename] [num_threads]\n\n";

        std::cout << "<cmd> can be one of the following:\n";
        std::cout << "  S0: the sphere\n";
        std::cout << "  C2: the cofiber of 2\n";
        std::cout << "  Ceta: the cofiber of eta\n";
        std::cout << "  Cnu: the cofiber of nu\n";
        std::cout << "  Csigma: the cofiber of sigma\n\n";

        std::cout << "Default values:\n";
        std::cout << "  t_max = " << t_max << "\n";
        std::cout << "  stem_max = " << stem_max << "\n";
        std::cout << "  db_filename = <X>_Adams_res.db\n";
        std::cout << "  num_threads = " << num_threads << "\n\n ";

        std::cout << "Version:\n  2.1 (2022-8-6)" << std::endl;
        return 0;
    }
    int index = 0;
    if (myio::load_arg(argc, argv, ++index, "X", X))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "t_max", t_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "stem_max", stem_max))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "db_filename", db_filename))
        return index;
    if (myio::load_op_arg(argc, argv, ++index, "num_threads", num_threads))
        return index;

    ut::FUTURE_NUM_THREADS = num_threads;
    bench::Timer timer;

    AdamsE2(X, t_max, stem_max, db_filename);

#ifndef MYDEPLOY
    bench::Counter::print();
#endif
    return 0;
}