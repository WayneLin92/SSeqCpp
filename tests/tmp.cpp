#include "algebras/groebner.h"
#include "algebras/myio.h"
#include <iostream>

void test()
{
    using namespace alg;
    auto h0 = Poly::Gen(0);
    auto h1 = Poly::Gen(1);
    auto h2 = Poly::Gen(2);
    auto h0_1 = Poly::Gen(3);
    auto b02 = Poly::Gen(4);
    auto b13 = Poly::Gen(5);

    Poly1d rels = {h0 * h1, h1 * h2, h0 * h0_1 + b02 * h2, h2 * h0_1 + b13 * h0, b02 * h2 * h2 + b13 * h0 * h0};

    auto gb = Groebner(100, {1, 1, 1, 1, 2, 2});
    gb.AddRels(rels, 100);
    
    auto gbm = GroebnerMod(&gb, 100, {0});
    auto x1 = h0;
    auto x2 = h2;
    auto x1m = Mod(x1, 0);
    auto x2m = Mod(x2, 0);

    int1d indices;
    gbm.ToSubMod({x1m, x2m}, 100, indices);

    for (auto& g : gbm.data())
        std::cout << g.Str() << '\n';
}

int main()
{
    test();

    return 0;
}