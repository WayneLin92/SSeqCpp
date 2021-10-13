#include "MayE2.h"
#include "algebras/benchmark.h"

bool is_bij(const alg::Deg& deg)
{
    if (deg.s == 2 && deg.t % 2 == 0 && deg.v % 2 == 0) {
        int t = deg.t / 2;
        int v = deg.v / 2;
        int j = bit_length(t);
        int i = j - (v + 1);
        if ((1 << j) - (1 << i) == t)
            return true;
    }
    return false;
}

int odd_part(int t)
{
    while (t % 2 == 0)
        t /= 2;
    return t;
}

bool compare(const alg::Deg& d1, const alg::Deg& d2)
{
    if (is_bij(d1) < is_bij(d2))
        return true;
    else if (is_bij(d1) == is_bij(d2)) {
        if (d1.s < d2.s) {
            return true;
        }
        else if (d1.s == d2.s) {
            if (odd_part(d1.t) < odd_part(d2.t))
                return true;
            else if (odd_part(d1.t) == odd_part(d2.t)) {
                if (d1.t < d2.t)
                    return true;
                else
                    return false;
            }
        }
    }
    return false;
}

/* This function reorders the generator by t and reproduce the Groebner basis */
std::pair<alg::array, alg::Poly1d> ReorderGens(const std::vector<alg::Deg>& gen_degs, const alg::Groebner& gb, int s_max)
{
    alg::array map_gen_id_inv = alg::range((int)gen_degs.size()); /* the i`th new generator is the old map_gen_id_inv[i]`th generator */
    std::sort(map_gen_id_inv.begin(), map_gen_id_inv.end(), [&gen_degs](int i, int j) { return compare(gen_degs[i], gen_degs[j]); });
    alg::array map_gen_id;
    map_gen_id.resize(gen_degs.size()); /* Generator id change: i -> map_gen_id[i] */
    alg::array gen_degs_new;
    for (int i = 0; i < (int)gen_degs.size(); ++i) {
        map_gen_id[map_gen_id_inv[i]] = i;
        gen_degs_new.push_back(gen_degs[map_gen_id_inv[i]].s);
    }

    alg::GbBuffer buffer;
    for (const alg::Poly& g : gb)
        buffer[get_deg(g, gen_degs).s].push_back(alg::subs(g, map_gen_id));
    alg::Groebner gb_new;
    alg::AddRelsB(gb_new, buffer, gen_degs_new, s_max, s_max);
    return std::make_pair(std::move(map_gen_id_inv), std::move(gb_new.gb));
}

/* Reorder and remove decomposables */
void ReorderHX(int n, int s_max)
{
    alg::Database db("/Users/weinanlin/MyData/Math_AlgTop/databases/HX9.db");
    std::string table_prefix = "HX" + std::to_string(n) + std::to_string(n);
    std::string table1_prefix = "HX" + std::to_string(n);

    try {
        db.execute_cmd("CREATE TABLE " + table1_prefix + "_generators (gen_id INTEGER PRIMARY KEY, gen_name TEXT UNIQUE, gen_diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
        db.execute_cmd("CREATE TABLE " + table1_prefix + "_relations (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
    }
    catch (MyException&) {
    }

    db.execute_cmd("DELETE FROM " + table1_prefix + "_generators");
    db.execute_cmd("DELETE FROM " + table1_prefix + "_relations");

    std::cout << "A=" << table_prefix << '\n';
    std::cout << "B=" << table1_prefix << '\n';

    /* Reorder and produce gen_degs1, gen_names1, gen_reprs1, gb1 */
    std::vector<alg::Deg> gen_degs = db.load_gen_degs(table_prefix + "_generators");
    std::vector<std::string> gen_names = db.load_gen_names(table_prefix + "_generators");
    alg::Poly1d gen_reprs = db.load_gen_reprs(table_prefix + "_generators");
    alg::Groebner gb = db.load_gb_s(table_prefix + "_relations", s_max);
    auto [map_gen_id_inv, gb1] = ReorderGens(gen_degs, gb, s_max);

    std::cout << "gb1.size()=" << gb1.size() << '\n';

    std::vector<alg::Deg> gen_degs1;
    std::vector<std::string> gen_names1;
    alg::Poly1d gen_reprs1;
    for (int i = 0; i < (int)gen_degs.size(); ++i) {
        gen_degs1.push_back(gen_degs[map_gen_id_inv[i]]);
        gen_names1.push_back(gen_names[map_gen_id_inv[i]]);
        gen_reprs1.push_back(gen_reprs[map_gen_id_inv[i]]);
    }

    /* Delete decomposables and produce gen_degs2, gen_names2, gen_reprs2, gb2 */
    std::vector<alg::Deg> gen_degs2;
    std::vector<std::string> gen_names2;
    alg::Poly1d gen_reprs2;
    alg::Poly1d gb2;

    alg::array indices_decomposables;
    alg::array map_gen_id1(gen_degs1.size());
    int count = 0;
    for (int i = 0; i < (int)gen_degs1.size(); ++i) {
        if (alg::Reduce({ { { i, 1 } } }, gb1) != alg::Poly{ { { i, 1 } } }) {
            indices_decomposables.push_back(i);
            ++count;
            map_gen_id1[i] = -1;
        }
        else {
            map_gen_id1[i] = i - count;
            gen_degs2.push_back(gen_degs1[i]);
            gen_names2.push_back(gen_names1[i]);
            gen_reprs2.push_back(gen_reprs1[i]);
        }
    }
    std::cout << "num of decomposables = " << indices_decomposables.size() << '\n';

    for (const alg::Poly& g : gb1)
        if (!std::binary_search(indices_decomposables.begin(), indices_decomposables.end(), g[0][0].gen))
            gb2.push_back(alg::subs(g, map_gen_id1));

    db.begin_transaction();
    db.save_generators(table1_prefix + "_generators", gen_names2, gen_degs2, gen_reprs2);
    db.save_gb(table1_prefix + "_relations", gb2, gen_degs2);
    db.end_transaction();
}
