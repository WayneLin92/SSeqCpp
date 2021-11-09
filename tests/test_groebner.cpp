#define CATCH_CONFIG_MAIN
#include "algebras/groebner.h"
#include "algebras/myio.h"
#include <catch2/catch.hpp>

namespace Catch {
template <>
struct StringMaker<alg::PolyLex>
{
    static std::string convert(alg::PolyLex const& p)
    {
        std::stringstream ss;
        ss << p.data;
        return ss.str();
    }
};
}  // namespace Catch

TEST_CASE("Compute the Reduction of f by g", "[alg::detail::Reduce]")
{
    constexpr auto GenExp = alg::PolyLex::GenExp;
    alg::PolyLex f1 = GenExp(0, 1) + GenExp(1, 2) + GenExp(2, 1);                                 /* x_0 + x_1^2 + x_2 */
    alg::PolyLex g1 = GenExp(1, 1) + GenExp(3, 1) + GenExp(4, 1);                                 /* x_1 + x_3 + x_4 */
    alg::PolyLex h1 = GenExp(0, 1) + GenExp(2, 1) + GenExp(1, 1) * (GenExp(3, 1) + GenExp(4, 1)); /* x_0 + x_2 + x_1x_3 + x_1x_4 */
    alg::detail::Reduce(f1, g1, 1);
    REQUIRE(f1 == h1);

    alg::PolyLex f2 = GenExp(0, 1) + GenExp(1, 2) + GenExp(1, 1) * GenExp(3, 1) + GenExp(2, 1); /* x_0 + x_1^2 + x_1x_3 + x_2 */
    alg::PolyLex g2 = GenExp(1, 1) + GenExp(3, 1);                                              /* x_1 + x_3 */
    alg::PolyLex h2 = GenExp(0, 1) + GenExp(2, 1);                                              /* x_0 + x_2 */
    alg::detail::Reduce(f2, g2, 1);
    REQUIRE(f2 == h2);
}

TEST_CASE("Reduce wrt lex ordering", "[alg::Reduce]")
{
    constexpr auto GenExp = alg::PolyLex::GenExp;
    alg::PolyLex p1 = GenExp(0, 2) + GenExp(1, 2); /* x_0^2 + x_1^2 */
    alg::PolyLex p2 = GenExp(1, 4) + GenExp(3, 1); /* x_1^4 + x_3 */
    alg::PolyLex1d polys = {p1, p2};
    alg::Groebner<alg::CmpLex> gb(polys);
    alg::PolyLex q1 = GenExp(0, 9);                /* x_0^9 */
    alg::PolyLex q2 = GenExp(0, 1) * GenExp(3, 2); /* x_0x_3^2 */
    REQUIRE(gb.Reduce(q1) == q2);
}

TEST_CASE("Compute a Groebner basis wrt lex ordering", "[alg::AddRelsV2Lex]")
{
    alg::array gen_degs = {1, 1, 1};
    alg::Poly p1 = {{{1, 2}}, {{0, 2}}}; /* x_0^2 + x_1^2 */
    alg::Poly p2 = {{{0, 3}}};           /* x_0^3 */
    alg::Poly1d polys = {p1, p2};
    alg::GroebnerLex gb;
    alg::AddRels(gb, polys, gen_degs, -1);
    alg::Poly q1 = {{{1, 4}}}; /* x_1^4 */
    alg::Poly q2 = {};         /* 0 */
    REQUIRE(alg::Reduce(q1, gb) == q2);
}

TEST_CASE("Compute a Groebner basis wrt lex ordering", "[alg::AddRels]")
{
    constexpr auto GenExp = alg::PolyLex::GenExp;
    alg::array gen_degs = {1, 1, 1};
    alg::PolyLex p1 = GenExp(0, 2) + GenExp(1, 2); /* x_0^2 + x_1^2 */
    alg::PolyLex p2 = GenExp(0, 3);                /* x_0^3 */
    alg::PolyLex1d polys = {p1, p2};
    alg::Groebner<alg::CmpLex> gb;
    alg::GbBuffer buffer;
    alg::AddRels(gb, polys, gen_degs, -1);
    alg::PolyLex q1 = GenExp(1, 4); /* x_1^4 */
    alg::PolyLex q2 = {};           /* 0 */
    // REQUIRE(gb.Reduce(q1) == q2);
}

TEST_CASE("Compute the Groebner basis of B7 wrt revlex ordering", "[alg::AddRels]")
{
    int n_max = 7;
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
    size_t answer = 65;
    // size_t answer = 163;
    REQUIRE(gb.size() == answer);
}

TEST_CASE("Compute the Groebner basis of B7 wrt revlex ordering", "[alg::AddRels]")
{
    using FnCmp = alg::CmpRevlex;
    using Poly = alg::Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner<FnCmp>;
    constexpr auto GenExp = Poly::GenExp;
    int n_max = 8;
    alg::array gen_degs;
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
                rel += GenExp(index1, 1) * GenExp(index2, 1);
            }
            rels.push_back(std::move(rel));
        }
    }

    Gb gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) { return p1.Deg(gen_degs) < p2.Deg(gen_degs); });
    alg::AddRels(gb, std::move(rels), gen_degs, -1);
    // size_t answer = 65; /* n_max = 7 */
    size_t answer = 163; /* n_max = 8 */
    // size_t answer = 462; /* n_max = 9 */
    REQUIRE(gb.size() == answer);
}

TEST_CASE("Computes the Groebner basis of B7 wrt lex ordering", "[alg::AddRelsV2Lex]")
{
    int n_max = 8;
    alg::array gen_degs;
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
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

    alg::GroebnerLex gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs](const alg::Poly& p1, const alg::Poly& p2) { return alg::get_deg(p1, gen_degs) < get_deg(p2, gen_degs); });
    alg::AddRelsV2(gb, std::move(rels), gen_degs, -1);
    // size_t answer = 17; /* n_max = 5 */
    // size_t answer = 78; /* n_max = 7 */
    size_t answer = 181; /* n_max = 8 */
    REQUIRE(gb.size() == answer);
}

TEST_CASE("Compute the Groebner basis of B7 wrt lex ordering", "[alg::AddRels]")
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
    // size_t answer = 78; /* n_max = 7 */
    size_t answer = 181; /* n_max = 8 */
    // size_t answer = 402; /* n_max = 9 */
    REQUIRE(gb.size() == answer);
}
