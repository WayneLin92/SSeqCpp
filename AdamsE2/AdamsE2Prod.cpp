#include "algebras/benchmark.h"
#include "algebras/utility.h"
#include "groebner_steenrod_const.h"

// std::vector<int> bench::Counter::counts_ = {0, 0, 0, 0};
int main()
{
    using namespace steenrod;
    bench::Timer timer;

#ifdef MYDEPLOY
    compute_products(100);
#else
    compute_products(50);
#endif

    compute_products_ind();
    compute_products(DEG_MAX);

    // bench::Counter::print();
    return 0;
}