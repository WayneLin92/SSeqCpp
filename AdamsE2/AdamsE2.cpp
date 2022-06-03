#include "algebras/benchmark.h"
#include "algebras/utility.h"
#include "groebner_steenrod.h"

#if !defined(MYDEPLOY) || defined(MYDEPLOY_TEST_FILE)
std::vector<int> bench::Counter::counts_ = {0, 0, 0};
#endif

void AdamsE2()
{
    using namespace steenrod;
    bench::Timer timer;

#ifdef MYDEPLOY
    int t_trunc = DEG_MAX_MULT;
#else
    int t_trunc = 100;
#endif
    Mod1d rels;
    for (int i = 0; i < 10; ++i) {
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
        if (MMilnor::P(i, i + 1).deg() > t_trunc)
            break;
    }
#ifdef MYDEPLOY
    auto gb = GroebnerMRes::load("AdamsE2.db", t_trunc);
#else
    auto gb = GroebnerMRes(t_trunc, {}, {}, {}, {}, 0, 0);
#endif
    AddRelsMRes(gb, rels, t_trunc);

    std::cout << "gb.dim_Ext()=" << gb.dim_Ext() << '\n';
    std::cout << "gb.dim_Gb()=" << gb.dim_Gb() << '\n';
}

#ifdef MYDEPLOY_TEST_FILE
/* Test the file saving and reading */
void AdamsE2_TestFile()
{
    using namespace steenrod;

    // int t_trunc1 = 61;
    int t_trunc = TEST_DEG;

    /*Mod1d rels;
    for (int i = 0; i < 10; ++i) {
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
        if (MMilnor::P(i, i + 1).deg() > t_trunc)
            break;
    }*/

    // GroebnerMRes::reset("AdamsE2.db");
    auto gb = GroebnerMRes::load("AdamsE2_right.db", t_trunc);
    // AddRelsMRes(gb, rels, t_trunc1);

    auto gb1 = GroebnerMRes::load("AdamsE2_wrong.db", t_trunc);
    std::cout << "gb.hash=" << gb.hash() << '\n';
    std::cout << "gb1.hash=" << gb1.hash() << '\n';

    /*AddRelsMRes(gb, rels, t_trunc);
    AddRelsMRes(gb1, rels, t_trunc);*/

    std::cout << "gb.dim_Ext()=" << gb.dim_Ext() << '\n';
    std::cout << "gb.dim_Gb()=" << gb.dim_Gb() << '\n';

    std::cout << "gb1.dim_Ext()=" << gb1.dim_Ext() << '\n';
    std::cout << "gb1.dim_Gb()=" << gb1.dim_Gb() << '\n';
}
#endif

int main()
{
#ifdef MYDEPLOY_TEST_FILE
    AdamsE2_TestFile();
#else
    AdamsE2();
#endif

#if !defined(MYDEPLOY) || defined(MYDEPLOY_TEST_FILE)
    bench::Counter::print();
#endif
    return 0;
}