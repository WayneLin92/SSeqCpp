#include "algebras/algebras.h"
#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/groebner.h"
#include "algebras/myexception.h"
#include <iostream>

int main()
{
    using Poly = alg::PolyLex;
    int n_max = 9;
    alg::array gen_degs_t;
    std::vector<alg::MayDeg> gen_degs;
    std::vector<std::string> gen_names;
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
            gen_degs_t.push_back((1 << j) - (1 << i));
            gen_degs.push_back(alg::MayDeg{1, (1 << j) - (1 << i), j - i});
            gen_names.push_back("R_{" + std::to_string(i) + std::to_string(j) + "}");
        }
    }
    alg::GroebnerLex gb(alg::DEG_MAX);
    myio::DbAlg db(std::string(myio::dir_db) + std::string("B9.db"));

    for (int n = 1; n <= n_max; ++n) {
        for (int m = 1; m <= n; ++m) {
            int j = n;
            int i = n - m;
            Poly rel;
            for (int k = i + 1; k < j; k++) {
                int a = (1 << k) - (1 << i);
                int b = (1 << j) - (1 << k);
                auto p1 = std::find(gen_degs_t.begin(), gen_degs_t.end(), a);
                auto p2 = std::find(gen_degs_t.begin(), gen_degs_t.end(), b);
                int index1 = int(p1 - gen_degs_t.begin());
                int index2 = int(p2 - gen_degs_t.begin());
                rel += Poly::Gen(index1) * Poly::Gen(index2);
            }

            alg::AddRels(gb, {rel}, alg::DEG_MAX, gen_degs_t);
            std::string table_prefix = "B" + std::to_string(n) + std::to_string(m);

            std::cout << table_prefix << '\n';

            db.begin_transaction();
            db.create_generators_and_delete(table_prefix);
            db.create_relations_and_delete(table_prefix);
            db.save_gen_names(table_prefix, gen_names);
            db.save_gen_maydegs(table_prefix, gen_degs);
            db.save_gb(table_prefix, gb, gen_degs);
            db.end_transaction();
        }
    }

    return 0;
}
