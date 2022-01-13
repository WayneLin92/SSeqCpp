#include "algebras/algebras.h"
#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/groebner.h"
#include "algebras/myexception.h"
#include <iostream>

int main()
{
    int n_max = 9;
    alg::array gen_degs_t;
    std::vector<alg::MayDeg> gen_degs;
    std::vector<std::string> gen_names;
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
            gen_degs_t.push_back((1 << j) - (1 << i));
            gen_degs.push_back(alg::MayDeg{ 1, (1 << j) - (1 << i), j - i });
            gen_names.push_back("R_{" + std::to_string(i) + std::to_string(j) + "}");
        }
    }
    alg::GroebnerLex gb;
    alg::Database db("/Users/weinanlin/MyData/Math_AlgTop/databases/B9.db");

    for (int n = 1; n <= n_max; ++n) {
        for (int m = 1; m <= n; ++m) {
            int j = n;
            int i = n - m;
            alg::Poly rel;
            for (int k = i + 1; k < j; k++) {
                int a = (1 << k) - (1 << i);
                int b = (1 << j) - (1 << k);
                auto p1 = std::find(gen_degs_t.begin(), gen_degs_t.end(), a);
                auto p2 = std::find(gen_degs_t.begin(), gen_degs_t.end(), b);
                int index1 = int(p1 - gen_degs_t.begin());
                int index2 = int(p2 - gen_degs_t.begin());
                rel = add(rel, alg::Poly{ { { index1, 1 } } } * alg::Poly{ { { index2, 1 } } });
            }

            alg::AddRelsV2(gb, { rel }, gen_degs_t, -1);
            std::string table_prefix = "B" + std::to_string(n) + std::to_string(m);

            std::cout << table_prefix << '\n';

            db.begin_transaction();
            try {
                db.execute_cmd("CREATE TABLE " + table_prefix + "_generators (gen_id INTEGER PRIMARY KEY, gen_name TEXT UNIQUE, s SMALLINT, t SMALLINT, v SMALLINT);");
                db.execute_cmd("CREATE TABLE " + table_prefix + "_relations (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
            }
            catch (MyException&) {
            }
            db.execute_cmd("delete from " + table_prefix + "_generators");
            db.execute_cmd("delete from " + table_prefix + "_relations");
            db.save_generators(table_prefix + "_generators", gen_names, gen_degs);
            db.save_gb(table_prefix + "_relations", gb.gb, gen_degs);
            db.end_transaction();
        }
    }

    return 0;
}
