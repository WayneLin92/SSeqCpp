#include "MayE2.h"

bool is_bij(const alg::MayDeg& deg)
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

bool compare(const alg::MayDeg& d1, const alg::MayDeg& d2)
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
std::pair<alg::array, alg::GroebnerRevlex> ReorderGens(const std::vector<alg::MayDeg>& gen_degs, const alg::GroebnerRevlex& gb, int t_max)
{
    alg::array map_gen_id_inv = ut::range((int)gen_degs.size()); /* the i`th new generator is the old map_gen_id_inv[i]`th generator */
    std::sort(map_gen_id_inv.begin(), map_gen_id_inv.end(), [&gen_degs](int i, int j) { return compare(gen_degs[i], gen_degs[j]); });
    alg::array map_gen_id;
    map_gen_id.resize(gen_degs.size()); /* Generator id change: i -> map_gen_id[i] */
    alg::array gen_degs_new;
    for (int i = 0; i < (int)gen_degs.size(); ++i) {
        map_gen_id[map_gen_id_inv[i]] = i;
        gen_degs_new.push_back(gen_degs[map_gen_id_inv[i]].t);
    }

    alg::GbBuffer buffer;
    alg::PolyRevlex1d rels;
    for (const alg::PolyRevlex& g : gb.data)
        rels.push_back(alg::subs<alg::CmpRevlex>(g.data, map_gen_id));
    alg::GroebnerRevlex gb_new;
    alg::AddRels(gb_new, rels, gen_degs_new, t_max);
    return std::make_pair(std::move(map_gen_id_inv), std::move(gb_new));
}

/* Reorder and remove decomposables */
void ReorderHX(int n, int t_max)
{
    myio::DbAlg db("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/HX9.db");
    std::string table_prefix = "HX" + std::to_string(n) + std::to_string(n);
    std::string table1_prefix = "HX" + std::to_string(n);

    std::cout << "A=" << table_prefix << '\n';
    std::cout << "B=" << table1_prefix << '\n';

    db.create_generators_and_delete(table1_prefix);
    db.create_relations_and_delete(table1_prefix);

    /* Reorder and produce gen_degs1, gen_names1, gen_reprs1, gb1 */
    std::vector<alg::MayDeg> gen_degs = db.load_gen_maydegs(table_prefix);
    std::vector<std::string> gen_names = db.load_gen_names(table_prefix);
    alg::PolyLex1d gen_reprs = db.load_gen_reprs<alg::CmpLex>(table_prefix);
    alg::GroebnerRevlex gb = db.load_gb<alg::CmpRevlex>(table_prefix, t_max);
    auto [map_gen_id_inv, gb1] = ReorderGens(gen_degs, gb, t_max);

    std::cout << "gb1.size()=" << gb1.size() << '\n';

    std::vector<alg::MayDeg> gen_degs1;
    std::vector<std::string> gen_names1;
    alg::PolyLex1d gen_reprs1;
    for (int i = 0; i < (int)gen_degs.size(); ++i) {
        gen_degs1.push_back(gen_degs[map_gen_id_inv[i]]);
        gen_names1.push_back(gen_names[map_gen_id_inv[i]]);
        gen_reprs1.push_back(gen_reprs[map_gen_id_inv[i]]);
    }

    /* Delete decomposables and produce gen_degs2, gen_names2, gen_reprs2, gb2 */
    std::vector<alg::MayDeg> gen_degs2;
    std::vector<std::string> gen_names2;
    alg::PolyLex1d gen_reprs2;
    alg::GroebnerRevlex gb2;

    alg::array indices_decomposables;
    alg::array map_gen_id1(gen_degs1.size());
    int count = 0;
    for (int i = 0; i < (int)gen_degs1.size(); ++i) {
        if (gb1.Reduce(alg::PolyRevlex::Gen(i)) != alg::PolyRevlex::Gen(i)) {
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

    for (const alg::PolyRevlex& g : gb1.data)
        if (!std::binary_search(indices_decomposables.begin(), indices_decomposables.end(), g.GetLead()[0].gen))
            gb2.push_back(alg::subs<alg::CmpRevlex>(g.data, map_gen_id1));

    db.begin_transaction();

    db.save_gen_maydegs(table1_prefix, gen_degs2);
    db.save_gen_names(table1_prefix, gen_names2);
    db.save_gen_reprs(table1_prefix, gen_reprs);
    db.save_gb(table1_prefix, gb2, gen_degs2);

    db.end_transaction();
}
