#include "algebras/groebner.h"
#include "algebras/myio.h"
#include "algebras/groebner_steenrod.h"
#include "algebras/dbAdamsSS.h"
#include <iostream>

void test()
{
    using namespace alg2;
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

void test1()
{
    using namespace algZ;
    auto x1 = Poly::Gen(1, 1, 3, true);
    auto x2 = Poly::Gen(2, 1, 3, true);
    auto x3 = Poly::Gen(3, 1, 3, true);
    auto x4 = Poly::Gen(4, 1, 3, true);
    auto x5 = Poly::Gen(5, 1, 3, true);
    auto y6 = Poly::Gen(6, 1, 1, false);
    auto y7 = Poly::Gen(7, 1, 1, false);
    auto two = Poly::twoTo(1);

    Mon m1 = (x1 * x3 * x5).data[0];
    Mon m2 = (x1 * x4 * x5).data[0];
    Mon m3 = (x1 * y6).data[0];

    Poly p1 = y6 * y6 + x1 + Mon::O(10);
    Poly p2 = y7 * y7 + x2 + Mon::O(4);
    Poly p3 = y7 * y7 * y6 * y6 + Mon::O(20);

    auto gb = Groebner(100, {{1, 1}, {3, 5}, {3, 5}, {3, 5}, {3, 5}, {3, 5}, {1, 2}, {1, 2}});
    ut::map_seq2d<int, 0> possEinf;
    gb.AddRels({Poly(Mon::two_x_square(6, 1)), Poly(Mon::two_x_square(7, 1))}, 100, possEinf);
    gb.AddRels({p1, p2, p3}, 100, possEinf);

    for (auto& rel : gb.data())
        std::cout << rel.GetLead().fil() << " " << rel << '\n';
}

void test2()
{
    using namespace alg2;
    auto x2 = Poly::Gen(2);
    auto x3 = Poly::Gen(3);
    auto x4 = Poly::Gen(4);

    auto gb = Groebner(100, {1, 1, 1, 1, 2, 2});
    gb.AddRels({x2 * x2 * x2}, 100);
    gb.AddRels({x3 * x4}, 100);
    gb.AddRels({x2 * x2 * x4}, 100);

    for (auto& g : gb.data())
        std::cout << g.Str() << '\n';
}

void test3()
{
    using namespace steenrod;
    auto a = Milnor::P(0, 1);
    auto b = Milnor::P(1, 8);
    uint32_t xi[] = {0, 0, 0, 0, 0, 0, 0, 1};
    std::cout << MMilnor::Xi(xi).StrXi() << '\n';

    /*for (size_t i = 0; i < MMILNOR_E_BITS; ++i)
        std::cout << int(MMILNOR_GEN_I[i]) << ' ' << int(MMILNOR_GEN_J[i]) << '\n';*/
}

void test4()
{
    using namespace algZ;
    Mod a = MMod({}, 0, 1) + MMod::O(1);
    Mod b = MMod({}, 0, 1) + MMod::O(2);

    std::cout << a.GetLead().fil() << '\n';
    std::cout << b.GetLead().fil() << '\n';
    std::cout << (a.GetLead() < b.GetLead()) << '\n';
    std::cout << a << '\n';
    std::cout << b << '\n';
}

void test5()
{
    using namespace algZ;
    auto x1 = Poly::Gen(1, 1, 3, true);
    auto x2 = Poly::Gen(2, 1, 3, true);
    auto x3 = Poly::Gen(3, 1, 3, true);
    auto x4 = Poly::Gen(4, 1, 3, true);
    auto x5 = Poly::Gen(5, 1, 3, true);
    auto y6 = Poly::Gen(6, 1, 1, false);
    auto y7 = Poly::Gen(7, 1, 1, false);
    auto y8 = Poly::Gen(8, 1, 1, false);
    auto x9 = Poly::Gen(9, 1, 3, true);
    auto x10 = Poly::Gen(10, 1, 3, true);
    auto two = Poly::twoTo(1);

    Poly p1 = y7 * y8 * y7 + Mon::O(20);
    Poly p2 = x1 * y6 * y8 * x9 + Mon::O(20);

    std::cout << p1 << '\n';

    std::cout << p1 * p2 << '\n';
    std::cout << p2 * p1 << '\n';
}

void test_for_each_pair_par()
{
    std::vector<std::vector<size_t>> p(19);
    ut::for_each_pair_par(p.size(), [&p](size_t i, size_t j) {
        p[i].push_back(j);
        p[j].push_back(i);
    });
    std::cout << "Done\n";
}

void test_groebner_steenrod()
{
    using namespace steenrod;
    int t_trunc = 50;
    int1d v_degs;
    std::vector<uint64_t> bcs;
    std::unordered_map<uint64_t, uint64_t> bc2v;
    for (int t = 1; t <= t_trunc; ++t) {
        for (int c = 0; c <= t / 3; ++c) {
            int b = t - 3 * c;
            v_degs.push_back(t);
            auto bc = ut::Bind((uint64_t)b, (uint64_t)c);
            bcs.push_back(bc);
            bc2v[bc] = bcs.size() - 1;
        }
    }
    Mod1d rels;
    Mod tmp;

    for (size_t i = 0; (1 << i) <= t_trunc; ++i) {
        size_t a = (size_t)1 << i;
        for (size_t j = 0; j < bcs.size(); ++j) {
            int d = v_degs[j];
            if (a + d > t_trunc)
                continue;
            uint64_t b, c;
            ut::UnBind(bcs[j], b, c);
            Mod rel = MMilnor::Sq((uint32_t)a) * MMod(MMilnor(), j);
            for (size_t n = 0; n <= (a + d) / 3; ++n) {
                size_t m = a + d - 3 * n; 
                if (!((a + 2 * c - 2 * n) & (b + c - n)) && c <= n && !(c & (n - c))) {
                    auto mn = ut::Bind((uint64_t)m, (uint64_t)n);
                    auto v_mn = bc2v.at(mn);
                    rel.iaddP(MMod(MMilnor(), v_mn), tmp);
                }
                //std::cout << "rel=" << rel << '\n';
            }
            rels.push_back(std::move(rel));
        }
    }

    Groebner gb(t_trunc, {}, v_degs);
    gb.AddRels(rels, t_trunc);
    gb.MinimizeOrderedGens();

    /*for (auto& rel : gb.data())
        std::cout << rel << '\n';*/
    std::cout << "v_degs=" << myio::StrCont("{", ",", "}", "{}", gb.v_degs(), [](int i) { return std::to_string(i); }) << '\n';
}

void minimize_rels()
{
    using namespace alg2;
    myio::DbAdamsSS db("C:/Users/lwnpk/Documents/Projects/algtop_cpp_build/bin/Release/main/S0_AdamsSS_t261.db");
    auto gb = Groebner(261, {}, db.load_gb("S0_AdamsE2", 261));
    auto gen_degs = db.load_gen_adamsdegs("S0_AdamsE2");
    int1d gen_degs_t;
    for (auto& d : gen_degs)
        gen_degs_t.push_back(d.t);
    auto gb1 = Groebner(261, gen_degs_t, {}, true);
    int1d min_rels;
    gb1.AddRels(gb.data(), 261, min_rels);

    std::cout << "gb.data().size()=" << gb.data().size() << '\n';
    std::cout << "gb1.data().size()=" << gb1.data().size() << '\n';
    std::cout << "min_rels.size()=" << min_rels.size() << '\n';
}

void reduce_gb()
{
    using namespace alg2;
    myio::DbAdamsSS db("C:/Users/lwnpk/Documents/Projects/algtop_cpp_build/bin/Release/main/S0_AdamsSS_t261.db");
    auto gb = Groebner(261, {}, db.load_gb("S0_AdamsE2", 261));
    auto gen_degs = db.load_gen_adamsdegs("S0_AdamsE2");
    gb.ReducedGb();
    std::map<AdamsDeg, Poly1d> gb_group_by_deg;
    for (auto& rel : gb.data()) {
        auto deg = GetDeg(rel.GetLead(), gen_degs);
        gb_group_by_deg[deg].push_back(rel);
    }
    db.begin_transaction();
    db.drop_and_create_relations("S0_AdamsE2");
    db.save_gb("S0_AdamsE2", gb_group_by_deg);
    db.end_transaction();
}

int main()
{
    minimize_rels();

    return 0;
}
