#include "groebner_res.h"

/****************************************************
 *                   S0
 ***************************************************/
void Coh_S0(int1d& v_degs, Mod1d& rels, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; (1 << i) <= t_max; ++i)
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
}

/****************************************************
 *                   tmf
 ***************************************************/
void Coh_tmf(int1d& v_degs, Mod1d& rels, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; i < 3 && ((1 << i) <= t_max); ++i)
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
}

/****************************************************
 *                   j
 ***************************************************/
void Coh_j(int1d& v_degs, Mod1d& rels, int t_max)
{
    v_degs = {0, 7};
    rels.clear();
    for (int i = 0; i < 3; ++i)
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
    Mod rel;
    rel = MMilnor::Sq(8) * MMod(MMilnor(), 0) + MMilnor::Sq(1) * MMod(MMilnor(), 1);
    rels.push_back(std::move(rel));

    rel = MMilnor::Sq(7) * MMod(MMilnor(), 1);
    rels.push_back(std::move(rel));

    rel = MMilnor::Sq(4) * (MMilnor::Sq(6) * MMod(MMilnor(), 1)) + MMilnor::Sq(6) * (MMilnor::Sq(4) * MMod(MMilnor(), 1));
    rels.push_back(std::move(rel));
}

/****************************************************
 *                   RPn
 ***************************************************/
void Coh_RP(int1d& v_degs, Mod1d& rels, int n1, int n2, int t_max)
{
    int1d v_degs2;
    for (int n = n1; n <= n2 && n <= t_max; ++n)
        v_degs2.push_back(n);

    Mod1d rels2;
    Mod tmp;
    /* Sq^{k} * x^n = (k, n - k) x^{n + k} */
    for (size_t i = 0; (1 << i) <= t_max; ++i) {
        size_t k = (size_t)1 << i;
        for (int n = n1; n <= n2 && n + k <= t_max; ++n) {
            Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), n - n1);
            if (k <= n && !(k & (n - k)) && n + k <= n2)
                rel.iaddP(MMod(MMilnor(), n + k - n1), tmp);
            rels2.push_back(std::move(rel));
        }
    }

    Groebner gb(t_max, {}, v_degs2);
    gb.AddRels(rels2, t_max);
    gb.MinimizeOrderedGens();

    v_degs = gb.v_degs();
    rels = gb.data();
}

/****************************************************
 *                   X2
 ***************************************************/
void Coh_X2(int1d& v_degs, Mod1d& rels, int t_max)
{
    int1d v_degs2;
    std::vector<uint64_t> pairs_bc;
    std::unordered_map<uint64_t, uint64_t> bc2v;
    for (int t = 0; t <= t_max; ++t) {
        for (int c = 0; c <= t / 3; ++c) {
            int b = t - 3 * c;
            v_degs2.push_back(t);
            auto bc = ut::Bind((uint64_t)b, (uint64_t)c);
            pairs_bc.push_back(bc);
            bc2v[bc] = pairs_bc.size() - 1;
        }
    }

    Mod1d rels2;
    Mod tmp;
    for (size_t i = 0; (1 << i) <= t_max; ++i) {
        size_t a = (size_t)1 << i;
        for (size_t j = 0; j < pairs_bc.size(); ++j) {
            int d = v_degs2[j];
            if (a + d > t_max)
                continue;
            uint64_t b, c;
            ut::UnBind(pairs_bc[j], b, c);
            Mod rel = MMilnor::Sq((uint32_t)a) * MMod(MMilnor(), j);
            for (size_t n = 0; n <= (a + d) / 3; ++n) {
                size_t m = a + d - 3 * n;
                if (a + 2 * c >= 2 * n && b + c >= n && !((a + 2 * c - 2 * n) & (b + c - n)) && c <= n && !(c & (n - c))) {
                    auto mn = ut::Bind((uint64_t)m, (uint64_t)n);
                    auto v_mn = bc2v.at(mn);
                    rel.iaddP(MMod(MMilnor(), v_mn), tmp);
                }
            }
            rels2.push_back(std::move(rel));
        }
    }

    Groebner gb(t_max, {}, v_degs2);
    gb.AddRels(rels2, t_max);
    gb.MinimizeOrderedGens();

    v_degs = gb.v_degs();
    rels = gb.data();
}
