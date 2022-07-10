#include "algebras/benchmark.h"
#include "algebras/utility.h"
#include "algebras/groebner.h"

void benchmark_B9_Revlex()
{
    using namespace alg;
    std::cout << "benchmark_B9_Revlex: \n";
    bench::Timer timer;
    constexpr auto Gen = Poly::Gen;

    int n_max = 9;
    int1d gen_degs;
    for (int d = 1; d <= n_max; d++) {
        for (int i = 0; i <= n_max - d; i++) {
            int j = i + d;
            gen_degs.push_back((1 << j) - (1 << i));
        }
    }
    Poly1d rels;
    for (int d = 2; d <= n_max; d++) {
        for (int i = 0; i <= n_max - d; i++) {
            int j = i + d;
            Poly rel;
            for (int k = i + 1; k < j; k++) {
                int a = (1 << k) - (1 << i);
                int b = (1 << j) - (1 << k);
                auto p1 = std::find(gen_degs.begin(), gen_degs.end(), a);
                auto p2 = std::find(gen_degs.begin(), gen_degs.end(), b);
                int index1 = int(p1 - gen_degs.begin());
                int index2 = int(p2 - gen_degs.begin());
                rel += Gen(index1, 1) * Gen(index2, 1);
            }
            rels.push_back(std::move(rel));
        }
    }

    Groebner gb(DEG_MAX, gen_degs);
    gb.AddRels(rels, DEG_MAX);
    size_t answer = 462;
    std::cout << "gb.size()=" << gb.size() << "==" << answer << '\n';
}

int main()
{
    benchmark_B9_Revlex();
    return 0;
}