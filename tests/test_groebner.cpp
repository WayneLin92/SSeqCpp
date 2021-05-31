#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <algebras/groebner.h>

TEST_CASE( "Reduce wrt lex ordering", "[alg::Reduce]" ) {
    alg::Poly p1 = {{{1, 2}}, {{0, 2}}}; /* x_0^2 + x_1^2 */
    alg::Poly p2 = {{{3, 1}}, {{1, 4}}}; /* x_1^4 + x_3 */
    alg::Poly1d polys = {p1, p2};
    alg::GroebnerLex gb(polys);
    alg::Poly q1 = {{{0, 9}}}; /* x_0^9 */
    alg::Poly q2 = {{{0, 1}, {3, 2}}}; /* x_0x_3^2 */
    REQUIRE( alg::Reduce(q1, gb) == q2 );
}

TEST_CASE( "Compute a Groebner basis wrt lex ordering", "[alg::AddRelsV2Lex]" ) {
    alg::array gen_degs = {1, 1, 1};
    alg::Poly p1 = {{{1, 2}}, {{0, 2}}}; /* x_0^2 + x_1^2 */
    alg::Poly p2 = {{{0, 3}}}; /* x_0^3 */
    alg::Poly1d polys = {p1, p2};
    alg::GroebnerLex gb;
    alg::AddRels(gb, polys, gen_degs, -1);
    alg::Poly q1 = {{{1, 4}}}; /* x_1^4 */
    alg::Poly q2 = {}; /* 0 */
    REQUIRE( alg::Reduce(q1, gb) == q2 );
}

TEST_CASE( "Compute the Groebner basis of B7 wrt revlex ordering", "[alg::AddRels]" ) {
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
                rel = add(rel, alg::Poly{ {{index1, 1}} } * alg::Poly{ {{index2, 1}} });
            }
            rels.push_back(std::move(rel));
        }
    }

    alg::Groebner gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs](const alg::Poly& p1, const alg::Poly& p2) {
        return alg::get_deg(p1, gen_degs) < get_deg(p2, gen_degs); });
    alg::AddRelsV2(gb, std::move(rels), gen_degs, -1);
    size_t answer = 65;
    //size_t answer = 163;
    REQUIRE( gb.size() == answer );
}

TEST_CASE( "Computes the Groebner basis of B7 wrt lex ordering", "[alg::AddRelsV2Lex]" ) {
    int n_max = 7;
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
                rel = add(rel, alg::Poly{ {{index1, 1}} } * alg::Poly{ {{index2, 1}} });
            }
            rels.push_back(std::move(rel));
        }
    }

    alg::GroebnerLex gb;
    std::sort(rels.begin(), rels.end(), [&gen_degs](const alg::Poly& p1, const alg::Poly& p2) {
        return alg::get_deg(p1, gen_degs) < get_deg(p2, gen_degs); });
    alg::AddRelsV2(gb, std::move(rels), gen_degs, -1);
    size_t answer = 78;
    //size_t answer = 181;
    REQUIRE( gb.size() == answer );
}
