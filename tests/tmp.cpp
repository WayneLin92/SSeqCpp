#include "algebras/benchmark.h"
#include "algebras/database.h"
#include "algebras/groebner.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include <execution>

void compare_computations()
{
    alg::Database db("/Users/weinanlin/MyData/Math_AlgTop/databases/HX9.db");
    alg::Database db_py("/Users/weinanlin/MyData/Math_AlgTop/databases/HX9_py_gen.db");

    alg::Groebner gb = db.load_gb_s("HX9_relations", 12);
    alg::Groebner gb_py = db_py.load_gb_s("HX9_relations", 12);

    std::vector<std::string> gen_names = db.load_gen_names("HX9_generators");
    std::vector<std::string> gen_names_py = db_py.load_gen_names("HX9_generators");

    alg::array map_gen_id(gen_names.size());
    alg::array map_gen_id_inv(gen_names.size());
    for (size_t i = 0; i < gen_names.size(); ++i) {
        size_t new_index = std::find(gen_names_py.begin(), gen_names_py.end(), gen_names[i]) - gen_names_py.begin();
        map_gen_id[i] = (int)new_index;
        map_gen_id_inv[new_index] = (int)i;
    }

    int count = 0;
    for (auto& rel : gb) {
        auto rel1 = alg::subs(rel, map_gen_id);
        if (!alg::Reduce(rel1, gb_py).empty()) {
            dump_PolyV2(std::cout, rel1, gen_names_py);
            std::cout << '\n';
            dump_PolyV2(std::cout, alg::Reduce(rel1, gb_py), gen_names_py);
            std::cout << "\n\n";

            ++count;
        }
    }
    std::cout << "count=" << count << '\n';
}

void compute_ann()
{
    int s_max = 12;
    alg::Database db("/Users/weinanlin/MyData/Math_AlgTop/databases/HX9.db");
    alg::Groebner gb = db.load_gb_s("HX9_relations", s_max);
    std::vector<std::string> gen_names = db.load_gen_names("HX9_generators");
    std::vector<alg::Deg> gen_degs = db.load_gen_degs("HX9_generators");
    alg::array gen_degs_s;
    for (auto p = gen_degs.begin(); p < gen_degs.end(); ++p)
        gen_degs_s.push_back(p->s);

    alg::Poly x = GetPolyByName(gen_names, "h_0(1,3)");
    alg::Poly2d a2d = alg::ann_seq(gb, {x}, gen_degs_s, s_max);
    alg::Poly1d a;
    for (alg::Poly1d& v : a2d)
        a.push_back(std::move(v[0]));

    for (auto& ai : a) {
        dump_PolyV2(std::cout, ai, gen_names);
        std::cout << '\n';
    }
}

void compute()
{
    int s_max = 12;
    alg::Database db("/Users/weinanlin/MyData/Math_AlgTop/databases/HX9.db");
    alg::Groebner gb = db.load_gb_s("HX9_relations", s_max);
    std::vector<std::string> gen_names = db.load_gen_names("HX9_generators");
    std::vector<alg::Deg> gen_degs = db.load_gen_degs("HX9_generators");
    alg::array gen_degs_s;
    for (auto p = gen_degs.begin(); p < gen_degs.end(); ++p)
        gen_degs_s.push_back(p->s);

    alg::Poly h0_1 = GetPolyByName(gen_names, "h_0(1)");
    alg::Poly h5_1 = GetPolyByName(gen_names, "h_5(1)");
    alg::Poly h1 = GetPolyByName(gen_names, "h_1");
    alg::Poly h5 = GetPolyByName(gen_names, "h_5");
    alg::Poly h7 = GetPolyByName(gen_names, "h_7");
    alg::Poly b26 = GetPolyByName(gen_names, "b_{26}");
    alg::Poly b27 = GetPolyByName(gen_names, "b_{27}");

    alg::Poly result = h1 * (h5_1 * b26 + h5 * h7 * b27);
    dump_PolyV2(std::cout, result, gen_names);
    std::cout << '\n';

    result = alg::Reduce(result, gb);
    dump_PolyV2(std::cout, result, gen_names);
    std::cout << '\n';
}

void benchmark_old()
{
    int n_max = 8;
    alg::array gen_degs;
    for (int d = 1; d <= n_max; d++) {
        for (int i = 0; i <= n_max - d; i++) {
            int j = i + d;
            gen_degs.push_back((1 << j) - (1 << i));
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
                auto p1 = std::find(gen_degs.begin(), gen_degs.end(), a);
                auto p2 = std::find(gen_degs.begin(), gen_degs.end(), b);
                int index1 = int(p1 - gen_degs.begin());
                int index2 = int(p2 - gen_degs.begin());
                rel = add(rel, alg::Poly{{{index1, 1}}} * alg::Poly{{{index2, 1}}});
            }
            rels.push_back(std::move(rel));
        }
    }

    alg::Groebner gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs](const alg::Poly& p1, const alg::Poly& p2) { return alg::get_deg(p1, gen_degs) < get_deg(p2, gen_degs); });
    alg::AddRelsV2(gb, std::move(rels), gen_degs, -1);
    size_t answer = 163;
    std::cout << "old: " << gb.size() << "==" << answer << '\n';
}

void benchmark_new()
{
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

    Gb gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) { return p1.Deg(gen_degs) < p2.Deg(gen_degs); });
    alg::AddRels(gb, std::move(rels), gen_degs, -1);
    size_t answer = 163;
    std::cout << "new: " << gb.size() << "==" << answer << '\n';
}

int main()
{
    Timer timer;
    int n = 1;
    for (int i = 0; i < n; ++i)
        benchmark_new();
    std::cout << timer.Elapsed() / n << '\n';

    return 0;
}