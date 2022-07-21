#include "algebras/benchmark.h"
#include "algebras/utility.h"
#include "algebras/myio.h"
#include "groebner_steenrod.h"

#include <sstream>
#include <cstring>

#ifndef MYDEPLOY
std::vector<int> bench::Counter::counts_ = {0, 0, 0};
#endif

void AdamsE2(int t_trunc, int stem_trunc, std::string& db_filename)
{
    using namespace steenrod;

    Mod1d rels;
    for (int i = 0; i < 10; ++i) {
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
        if (MMilnor::P(i, i + 1).deg() > t_trunc)
            break;
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
    int t_max = 392, stem_max = 261;
    std::string db_filename = "AdamsE2.db";
    auto num_threads = 128;

    if (argc >= 2 && strcmp(argv[1], "-h") == 0) {
        std::cout << "Usage:\n  AdamsE2 <t_max> <stem_max> <db_filename> <num_threads>\n\n";
        std::cout << "Default values:\n";

        std::cout << "  t_max = " << t_max << "\n";
        std::cout << "  stem_max = " << stem_max << "\n";
        std::cout << "  db_filename = " << db_filename << "\n";
        std::cout << "  num_threads = " << num_threads << "\n\n ";

        std::cout << "Version:\n  2.1 (2022-7-10)" << std::endl;
        return 0;
    }
    int index = 0;
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

    AdamsE2(t_max, stem_max, db_filename);

#ifndef MYDEPLOY
    bench::Counter::print();
#endif
    return 0;
}