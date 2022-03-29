#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/utility.h"
#include "groebner_steenrod.h"

std::vector<double> bench::Timer::counts_ = {};
std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0, 0};

void AdamsE2()
{
    using namespace steenrod;

    bench::Timer timer;

    int t_trunc = 100;
    ModCpt1d rels;
    for (int i = 0; i < 10; ++i) {
        rels.push_back(MModCpt(MMay::P(i, i + 1), 0));
        if (MMay::P(i, i + 1).deg() > t_trunc)
            break;
    }
    auto gb = GroebnerMRes::load("AdamsE2.db", t_trunc);
    AddRelsMRes(gb, rels, t_trunc);
}
struct A
{
    short a, b;
};
constexpr size_t n = sizeof(A) * 8;

int main()
{
    // steenrod::generate_basis("AdamsE2.db", 250);
    AdamsE2();
    std::cout << "counts_=\n";
    bench::Counter::print_counts_();
    std::cout << "timers_=\n";
    bench::Timer::print_counts_();
    return 0;
}