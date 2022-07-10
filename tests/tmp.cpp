#include "algebras/groebner.h"
#include "algebras/myio.h"
#include <iostream>

void test()
{
    using namespace alg;
    auto x0 = Poly::Gen(0);
    auto x1 = Poly::Gen(1);
    auto f = x0 * x0 + x1 * x1;
    auto g = pow(x1, 3);

    Poly1d rels = {f, g};

    auto gb = Groebner(100, {1, 1});
    gb.AddRels(rels, 100);
    for (auto& rel : gb.data())
        std::cout << rel << '\n';
}

int main()
{
    test();

    return 0;
}