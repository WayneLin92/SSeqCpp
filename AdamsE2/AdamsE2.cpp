#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/utility.h"
#include "groebner_steenrod.h"

std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0, 0};

void AdamsE2()
{
    using namespace steenrod;

    bench::Timer timer;

    int t_trunc = 70;
    Mod1d rels;
    for (int i = 0; i < 10; ++i) {
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
        if (MMilnor::P(i, i + 1).deg() > t_trunc)
            break;
    }
    auto gb = GroebnerMRes::load("AdamsE2.db", t_trunc);
    size_t dim = AddRelsMRes(gb, rels, t_trunc);
    std::cout << "dim=" << dim << '\n';

    auto& data = gb.data();
    size_t size_gb = 0;
    for (size_t s = 0; s < data.size(); ++s)
        size_gb += data[s].size();
    std::cout << "size_gb=" << size_gb << '\n';

    //for (size_t s = 0; s < data.size(); ++s) {
    //    std::cout << "s=" << s << '\n';
    //    for (size_t i = 0; i < data[s].size(); ++i) {
    //        if (data[s][i].x1.LF().GetLead().v() == 3 /*|| data[s][i].x2.data.size() == 1*/)
    //        //if (data[s][i].x1.LF().data.size() > 1)
    //            std::cout << data[s][i].x1.LF() /*<< " = " << data[s][i].x2*/ /*<< " = " << data[s][i].x2m*/ << '\n';
    //    }
    //}

    /* std::cout << "\ngb_x2m:\n";
    auto& data_x2m = gb.data_x2m();
    for (size_t s = 0; s < data_x2m.size(); ++s) {
        std::cout << "s=" << s << '\n';
        for (size_t i = 0; i < data_x2m[s].size(); ++i) {
            std::cout << data_x2m[s][i] << '\n';
        }
    }*/
}

int main()
{
    // steenrod::generate_basis("AdamsE2.db", 250);
    AdamsE2();
    bench::Counter::print();
    return 0;
}