#include "algebras/benchmark.h"
#include "algebras/utility.h"
#include "groebner_steenrod_const.h"

int main()
{
    using namespace steenrod;
    bench::Timer timer;

    compute_products();

    return 0;
}