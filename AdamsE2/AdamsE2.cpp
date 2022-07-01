#include "algebras/benchmark.h"
#include "algebras/utility.h"
#include "groebner_steenrod.h"

#ifndef TO_GUOZHEN
std::vector<int> bench::Counter::counts_ = {0, 0, 0};
#endif

void AdamsE2()
{
    using namespace steenrod;
    bench::Timer timer;

#ifdef TO_GUOZHEN
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
#ifdef TO_GUOZHEN
    auto gb = SteenrodMRes::load("AdamsE2.db", t_trunc);
#else
    auto gb = SteenrodMRes(t_trunc, {}, {}, {}, {}, 0, 0);
#endif
    ResolveMRes(gb, rels, t_trunc);

    std::cout << "gb.dim_Ext()=" << gb.dim_Ext() << '\n';
    std::cout << "gb.dim_Gb()=" << gb.dim_Gb() << '\n';
}

int main()
{
    AdamsE2();

#ifndef TO_GUOZHEN
    bench::Counter::print();
#endif
    return 0;
}