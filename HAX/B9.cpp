#include "algebras/algebras.h"
#include "algebras/groebner.h"
#include "algebras/benchmark.h"
#include "algebras/database.h"
#include <iostream>

int main()
{
    int n_max = 9;
    alg::array gen_degs_t;
    std::vector<alg::Deg> gen_degs;
    std::vector<std::string> gen_names;
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
            gen_degs_t.push_back((1 << j) - (1 << i));
            gen_degs.push_back(alg::Deg{1, (1 << j) - (1 << i), j - i});
            gen_names.push_back("R_{" + std::to_string(i) + std::to_string(j) + "}");
        }
    }
    alg::Poly1d rels;
    for (int d = 2; d <= n_max; d++) {
        for (int i = 0; i <= n_max - d; i++) {
            int j = i + d;
            alg::Poly rel;
            for (int k = i + 1; k < j; k++) {
                int a = (1 << k) - (1 << i);
                int b = (1 << j) - (1 << k);
                auto p1 = std::find(gen_degs_t.begin(), gen_degs_t.end(), a);
                auto p2 = std::find(gen_degs_t.begin(), gen_degs_t.end(), b);
                int index1 = int(p1 - gen_degs_t.begin());
                int index2 = int(p2 - gen_degs_t.begin());
                rel = add(rel, alg::Poly{ {{index1, 1}} } * alg::Poly{ {{index2, 1}} });
            }
            rels.push_back(std::move(rel));
        }
    }

    alg::GroebnerLex gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs_t](const alg::Poly& p1, const alg::Poly& p2) {
        return alg::get_deg(p1, gen_degs_t) < get_deg(p2, gen_degs_t); });
    alg::AddRelsV2(gb, std::move(rels), gen_degs_t, -1);

    Database db("/Users/weinanlin/MyData/Math_AlgTop/databases/B9.db");
    db.begin_transaction();
    db.execute_cmd("delete from B9_generators");
    db.execute_cmd("delete from B9_relations");
    db.save_generators("B9_generators", gen_names, gen_degs);
    db.save_gb("B9_relations", gb.gb, gen_degs);
    db.end_transaction();
    
    return 0;
}
