#include "algebras/algebras.h"
#include "algebras/myio.h"


int main()
{
    alg::Mon m1 = {{1, 1}};
    alg::Mon m2 = {{2, 1}};
    std::cout << mul(mul(m1, m2), m2) << '\n';
    return 0;
}