#include "algebras/database.h"
#include "algebras/groebner.h"
#include "algebras/myio.h"

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
        map_gen_id[i] = new_index;
        map_gen_id_inv[new_index] = i;
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
    alg::Poly2d a2d = alg::ann_seq(gb, { x }, gen_degs_s, s_max);
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

int main()
{
    compute_ann();
    return 0;
}