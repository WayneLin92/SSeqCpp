#include "B9.h"
#include "algebras/algebras.h"
#include "algebras/groebner.h"
#include "algebras/benchmark.h"
#include "algebras/database.h"
#include <iostream>

int main()
{
    Timer timer;
    int n_max = 8;

    std::vector<Deg> gen_degs;
    array gen_degs_t;
    for (int j = 1; j <= n_max; ++j) {
        for (int i = j - 1; i >= 0; --i) {
            gen_degs.push_back(Deg{1, (1 << j) - (1 << i), j - i});
            gen_degs_t.push_back((1 << j) - (1 << i));
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
                auto p1 = std::find(gen_degs_t.begin(), gen_degs_t.end(), a);
                auto p2 = std::find(gen_degs_t.begin(), gen_degs_t.end(), b);
                int index1 = int(p1 - gen_degs_t.begin());
                int index2 = int(p2 - gen_degs_t.begin());
                rel = add(rel, Poly{ {{index1, 1}} } * Poly{ {{index2, 1}} });
            }
            rels.push_back(std::move(rel));
        }
    }

    grbn::GbWithCache gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) {
        return get_deg_t(p1, gen_degs) < get_deg_t(p2, gen_degs); });
    grbn::AddRels(gb, std::move(rels), gen_degs_t, -1);
    std::cout << "t_max=" << get_deg_t(gb.gb.back(), gen_degs) << '\n';
    size_t gb_size = gb.size();
    size_t answer = 163;
    std::cout << "gb_size=" << gb_size << '\n';
    return 0;
}