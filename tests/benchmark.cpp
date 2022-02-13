#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/utility.h"

void benchmark_B9_Lex()
{
    std::cout << "benchmark_B9_Lex: \n";
    Timer timer;

    using FnCmp = alg::CmpLex;
    using Poly = alg::Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner<FnCmp>;
    constexpr auto GenExp = Poly::GenExp;

    int n_max = 9;
    alg::array gen_degs;
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
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
                rel += GenExp(index1, 1) * GenExp(index2, 1);
            }
            rels.push_back(std::move(rel));
        }
    }

    Gb gb(alg::DEG_MAX, alg::CRI_ON);
    std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) { return p1.GetDeg(gen_degs) < p2.GetDeg(gen_degs); });
    alg::AddRels(gb, std::move(rels), alg::DEG_MAX, gen_degs);
    size_t answer = 402;
    std::cout << "new: " << gb.size() << "==" << answer << '\n';
}

void benchmark_B9_Revlex()
{
    std::cout << "benchmark_B9_Revlex: \n";
    Timer timer;

    using FnCmp = alg::CmpRevlex;
    using Poly = alg::Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner<FnCmp>;
    constexpr auto GenExp = Poly::GenExp;

    int n_max = 9;
    alg::array gen_degs;
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
                rel += GenExp(index1, 1) * GenExp(index2, 1);
            }
            rels.push_back(std::move(rel));
        }
    }

    Gb gb(alg::DEG_MAX, alg::CRI_ON);
    std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) { return p1.GetDeg(gen_degs) < p2.GetDeg(gen_degs); });
    alg::AddRels(gb, std::move(rels), alg::DEG_MAX, gen_degs);
    size_t answer = 462;
    std::cout << "new: " << gb.size() << "==" << answer << '\n';
}

int main()
{
    benchmark_B9_Lex();
    benchmark_B9_Revlex();
    return 0;
}