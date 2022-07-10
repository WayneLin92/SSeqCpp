#define CATCH_CONFIG_MAIN
#include "algebras/groebner.h"
#include "algebras/myio.h"
#include <catch2/catch.hpp>

namespace Catch {
template <>
struct StringMaker<alg::Poly>
{
    static std::string convert(alg::Poly const& p)
    {
        std::stringstream ss;
        ss << p;
        return ss.str();
    }
};
}  // namespace Catch

TEST_CASE("Compute the Groebner basis of B7 wrt revlex ordering", "[alg::AddRels]")
{
    using namespace alg;
    constexpr auto Gen = Poly::Gen;
#ifndef NDEBUG
    int n_max = 7;
#else
    int n_max = 8;
#endif
    int1d gen_degs;
    for (int d = 1; d <= n_max; d++) {
        for (int i = 0; i <= n_max - d; i++) {
            int j = i + d;
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
                rel += Gen(index1, 1) * Gen(index2, 1);
            }
            rels.push_back(std::move(rel));
        }
    }

    Groebner gb(DEG_MAX, gen_degs);
    gb.AddRels(rels, DEG_MAX);
#ifndef NDEBUG
    size_t answer = 65; /* n_max = 7 */
#else
    size_t answer = 163; /* n_max = 8 */
    // size_t answer = 462; /* n_max = 9 */
#endif
    REQUIRE(gb.size() == answer);
}

//TEST_CASE("Compute the Groebner basis of B7 wrt lex ordering", "[alg::AddRels]")
//{
//    using FnCmp = alg::CmpLex;
//    using Poly = alg::Polynomial<FnCmp>;
//    using Poly1d = std::vector<Poly>;
//    using Gb = alg::Groebner<FnCmp>;
//    constexpr auto GenExp = Poly::GenExp;
//
//#ifndef NDEBUG
//    int n_max = 7;
//#else
//    int n_max = 8;
//#endif
//    alg::array gen_degs;
//    for (int i = 0; i < n_max; ++i) {
//        for (int j = i + 1; j <= n_max; ++j) {
//            gen_degs.push_back((1 << j) - (1 << i));
//        }
//    }
//    Poly1d rels;
//    for (int d = 2; d <= n_max; d++) {
//        for (int i = 0; i <= n_max - d; i++) {
//            int j = i + d;
//            Poly rel;
//            for (int k = i + 1; k < j; k++) {
//                int a = (1 << k) - (1 << i);
//                int b = (1 << j) - (1 << k);
//                auto p1 = std::find(gen_degs.begin(), gen_degs.end(), a);
//                auto p2 = std::find(gen_degs.begin(), gen_degs.end(), b);
//                int index1 = int(p1 - gen_degs.begin());
//                int index2 = int(p2 - gen_degs.begin());
//                rel += GenExp(index1, 1) * GenExp(index2, 1);
//            }
//            rels.push_back(std::move(rel));
//        }
//    }
//
//    Gb gb(alg::DEG_MAX);
//    std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) { return p1.GetDeg(gen_degs) < p2.GetDeg(gen_degs); });
//    alg::AddRels(gb, std::move(rels), alg::DEG_MAX, gen_degs);
//
//#ifndef NDEBUG
//    size_t answer = 78; /* n_max = 7 */
//#else
//    size_t answer = 181; /* n_max = 8 */
//    // size_t answer = 402; /* n_max = 9 */
//#endif
//    REQUIRE(gb.size() == answer);
//}
//
//TEST_CASE("Compute annihilator with FnCmp=CmpLex", "[alg::AnnSeq]")
//{
//    using FnCmp = alg::CmpLex;
//    using Poly = alg::Polynomial<FnCmp>;
//    using Poly1d = std::vector<Poly>;
//    using Gb = alg::Groebner<FnCmp>;
//    constexpr auto GenExp = Poly::GenExp;
//
//    alg::array gen_degs = {1, 1, 1};
//    Poly p1 = GenExp(0, 2) + GenExp(1, 2); /* x_0^2 + x_1^2 */
//    Poly p2 = GenExp(0, 1) * GenExp(2, 1) + GenExp(1, 1) * GenExp(2, 1);                /* x_0x_2 + x_1x_2 */
//    Poly1d polys = {p1, p2};
//    Gb gb(alg::DEG_MAX);
//    alg::AddRels(gb, polys, alg::DEG_MAX, gen_degs);
//
//    auto x = GenExp(0, 1) + GenExp(1, 1);
//
//    auto ann = alg::AnnSeq(gb, {x}, gen_degs, alg::DEG_MAX);
//    /*for (auto& v : ann)
//        std::cout << v[0].data << '\n';*/
//
//    REQUIRE(ann.size() == 2);
//}
//
//TEST_CASE("Compute annihilator with FnCmp=CmpLex V2", "[alg::AnnSeq]")
//{
//    using FnCmp = alg::CmpLex;
//    using Poly = alg::Polynomial<FnCmp>;
//    using Poly1d = std::vector<Poly>;
//    using Gb = alg::Groebner<FnCmp>;
//    constexpr auto GenExp = Poly::GenExp;
//
//    alg::array gen_degs = {1, 1, 1, 1};
//    Poly p1 = GenExp(0, 1) * GenExp(3, 1) + GenExp(1, 2); /* x_0x_3 + x_1^2 */
//    Poly p2 = GenExp(1, 2) + GenExp(2, 1) * GenExp(3, 1); /* x_1^2 + x_2x_3 */
//    Poly1d polys = {p1, p2};
//    Gb gb(alg::DEG_MAX);
//    alg::AddRels(gb, polys, alg::DEG_MAX, gen_degs);
//
//    auto x = GenExp(3, 1);
//
//    auto ann = alg::AnnSeq(gb, {x}, gen_degs, alg::DEG_MAX);
//    /*for (auto& v : ann)
//        std::cout << v[0].data << '\n';*/
//
//    REQUIRE(ann.size() == 1);
//}
