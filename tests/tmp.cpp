#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/steenrod.h"
#include "algebras/thread_pool.h"
#include "algebras/utility.h"

void compare_computations()
{
    myio::DbAlg db("/Users/weinanlin/MyData/Math_AlgTop/databases/HX9.db");
    myio::DbAlg db_py("/Users/weinanlin/MyData/Math_AlgTop/databases/HX9_py_gen.db");

    alg::GroebnerRevlex gb = db.load_gb<alg::CmpRevlex>("HX9_relations", 12);
    alg::GroebnerRevlex gb_py = db_py.load_gb<alg::CmpRevlex>("HX9_relations", 12);

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
    for (auto& rel : gb.data()) {
        auto rel1 = alg::subs<alg::CmpRevlex>(rel.data, map_gen_id);
        if (gb_py.Reduce(rel1)) {
            /*std::cout << StrPoly(rel1.data, gen_names_py) << '\n';
            std::cout << StrPoly(gb_py.Reduce(rel1).data, gen_names_py) << "\n\n";*/
            ++count;
        }
    }
    std::cout << "count=" << count << '\n';
}

void test_homology()
{
    bench::Timer timer;
    using FnCmp = alg::CmpLex;
    using Poly = alg::Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner<FnCmp>;

    int n_max = 5;
    alg::MayDeg1d gen_degs;
    alg::array gen_degs_t;
    std::vector<std::string> gen_names;
    Poly1d gen_diffs;
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
            gen_degs.push_back({1, (1 << j) - (1 << i), j - i - 1});
            gen_degs_t.push_back(gen_degs.back().t);
            gen_names.push_back("R_{" + std::to_string(i) + std::to_string(j) + "}");
        }
    }
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
            Poly diff;
            for (int k = i + 1; k < j; k++) {
                int a = (1 << k) - (1 << i);
                int b = (1 << j) - (1 << k);
                auto p1 = std::find(gen_degs_t.begin(), gen_degs_t.end(), a);
                auto p2 = std::find(gen_degs_t.begin(), gen_degs_t.end(), b);
                int index1 = int(p1 - gen_degs_t.begin());
                int index2 = int(p2 - gen_degs_t.begin());
                diff += Poly::Gen(index1) * Poly::Gen(index2);
                // std::cout << i << j << k << " " << gen_names.back() << " " << gen_names[index1] << " " << gen_names[index2] << '\n';
            }
            gen_diffs.push_back(std::move(diff));
        }
    }
    Gb gb(alg::DEG_MAX);

    alg::GroebnerRevlex gb_h(alg::DEG_MAX);
    alg::MayDeg1d gen_degs_h;
    Poly1d gen_repr_h;
    std::vector<std::string> gen_names_h;
    alg::Homology(gb, gen_degs, gen_diffs, alg::MayDeg{1, 0, -1}, gb_h, gen_degs_h, gen_names_h, gen_repr_h, 80);

    /*for (size_t i = 0; i < gen_names_h.size(); ++i)
        std::cout << "t=" << gen_repr_h[i].GetMayDeg(gen_degs) << " " << gen_names_h[i] << "=" << StrPoly(gen_repr_h[i].data, gen_names) << '\n';*/
    std::cout << '\n';
    /*for (auto& g : gb_h.data())
        std::cout << StrPoly(g.data, gen_names_h) << "=0\n";*/
}

void TestMilnorProduct()
{
    using namespace steenrod;
    uint32_t xi1[XI_MAX] = {0, 0, 1, 2, 2};
    uint32_t xi2[XI_MAX] = {2, 0, 0, 0, 1};
    auto a1 = MMilnor::Xi(xi1);
    auto a2 = MMilnor::Xi(xi2);
    auto prod = a1 * a2;
    std::cout << prod << '\n';
}

constexpr int N = 10000000;
constexpr int M = N/24;

void Bench1()
{
    std::vector<double> arr_cos(N);

    std::cout << "bench1:\n";
    bench::Timer timer;
    for (int i = 0; i < N; ++i) {
        arr_cos[i] = std::cos(std::cos(i));
        arr_cos[i] += std::cos(std::cos(i));
        arr_cos[i] += std::cos(std::cos(i));
        arr_cos[i] += std::cos(std::cos(i));
        arr_cos[i] += std::cos(std::cos(i));
        arr_cos[i] += std::cos(std::cos(i));
    }
    double sum = 0;
    for (int i = 0; i < N; ++i)
        sum += arr_cos[i];
    std::cout << "sum=" << sum << std::endl;
}

void Bench2()
{
    ut::ThreadPool pool(24);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    std::vector<double> arr_cos(N);

    std::cout << "bench2:\n";
    bench::Timer timer;
    for (int i = 0; i < N; i += M) {
        int end = std::min(i + M, N);
        pool.push_task([&arr_cos, i, end]() {
            for (int j = i; j < end; ++j) {
                arr_cos[j] = std::cos(std::cos(j));
                arr_cos[j] += std::cos(std::cos(j));
                arr_cos[j] += std::cos(std::cos(j));
                arr_cos[j] += std::cos(std::cos(j));
                arr_cos[j] += std::cos(std::cos(j));
                arr_cos[j] += std::cos(std::cos(j));
            }
        });
    }
    pool.wait_for_tasks();
    double sum = 0;
    for (int i = 0; i < N; ++i)
        sum += arr_cos[i];
    std::cout << "sum=" << sum << std::endl;
}

void TestThreadPool()
{
    //Bench1();
    Bench2();
    //Bench1();
    //Bench2();
}

int main()
{
    TestThreadPool();

    return 0;
}