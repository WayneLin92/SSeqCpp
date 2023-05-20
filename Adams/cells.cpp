#include "groebner_res.h"
#include <regex>
#include <fmt/core.h>

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
 *                   ko
 ***************************************************/
void Coh_ko(int1d& v_degs, Mod1d& rels, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; i < 2 && ((1 << i) <= t_max); ++i)
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
}

/****************************************************
 *                   X2
 ***************************************************/
void Coh_X2(int1d& v_degs, Mod1d& rels, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int j = 3; j <= 5; ++j)
        for (int i = 0; (((1 << j) - 1) << i) <= t_max; ++i)
            rels.push_back(MMod(MMilnor::P(i, i + j), 0));
}

/****************************************************
 *                    Ch
 *              Cofiber of h_n, n=0,1,2,3
 ***************************************************/
void Coh_Chn(int1d& v_degs, Mod1d& rels, int n, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; (1 << i) <= t_max; ++i)
        if (i != n)
            rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
    for (int i = 0; (1 << i) + (1 << n) <= t_max; ++i)
        rels.push_back(MMilnor::P(i, i + 1) * MMod(MMilnor::P(n, n + 1), 0));
}

/****************************************************
 *               C(2, h_n), n=4,5,6
 ***************************************************/
void Coh_C2hn(int1d& v_degs, Mod1d& rels, int n, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; (1 << i) <= t_max; ++i)
        if (i != 0 && i != n)
            rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
    for (int j = 0; j <= n; j += n)
        for (int i = 0; (1 << i) + (1 << j) <= t_max; ++i)
            rels.push_back(MMilnor::P(i, i + 1) * MMod(MMilnor::P(j, j + 1), 0));
}

void Coh_C2h5(int1d& v_degs, Mod1d& rels, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; (1 << i) <= t_max; ++i)
        if (i != 0 && i != 5)
            rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
    for (int n = 0; n <= 5; n += 5)
        for (int i = 0; (1 << i) + (1 << n) <= t_max; ++i)
            rels.push_back(MMilnor::P(i, i + 1) * MMod(MMilnor::P(n, n + 1), 0));
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
int Period_RP(int n)
{
    MyException::Assert(n >= 0, "n >= 0 in phi(n)");
    int phi = (n / 8) * 4, residue = n % 8;
    if (residue <= 1)
        phi += residue;
    else if (residue <= 3)
        phi += 2;
    else if (residue <= 7)
        phi += 3;
    return 1 << phi;
}

void normalize_RP(int n1, int n2, int& n1_, int& n2_)
{
    int T = Period_RP(n2 - n1);
    n1_ = 1 + ((n1 - 1) % T);
    n2_ = n2 - (n1 - n1_);
}

void Coh_RP(int1d& v_degs, Mod1d& rels, Mod1d& convert_v, int n1, int n2, int t_max)
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
    gb.MinimizeOrderedGens(convert_v);

    v_degs = gb.v_degs();
    rels = gb.data();
}

/****************************************************
                     # Maps
 ***************************************************/
void SetCohMapImages(std::string& cw1, std::string& cw2, Mod1d& images, int& sus)
{
    std::regex words_regex("^RP(\\d+)_(\\d+)$"); /* match example: RP1_4 */
    std::smatch match_RP;
    std::regex_search(cw1, match_RP, words_regex);
    std::ssub_match sub1 = match_RP[1], sub2 = match_RP[2];
    if (cw1 == "S0") {
        if (cw2 == "tmf" || cw2 == "ko" || cw2 == "X2") {
            images = {MMod(MMilnor(), 0)};
            sus = 0;
            return;
        }
    }
    else if (cw2 == "S0") {
        if (cw1 == "C2") {
            images = {MMod(MMilnor::P(0, 1), 0)};
            sus = 1;
            return;
        }
        else if (cw1 == "Ceta") {
            images = {MMod(MMilnor::P(1, 2), 0)};
            sus = 2;
            return;
        }
        else if (cw1 == "Cnu") {
            images = {MMod(MMilnor::P(2, 3), 0)};
            sus = 4;
            return;
        }
        else if (cw1 == "Csigma") {
            images = {MMod(MMilnor::P(3, 4), 0)};
            sus = 8;
            return;
        }
        else if (cw1 == "C2h4") {
            images = {MMod(MMilnor::P(4, 5), 0)};
            sus = 16;
            return;
        }
        else if (cw1 == "C2h5") {
            images = {MMod(MMilnor::P(5, 6), 0)};
            sus = 32;
            return;
        }
        else if (cw1 == "C2h6") {
            images = {MMod(MMilnor::P(6, 7), 0)};
            sus = 64;
            return;
        }
    }
    else if (cw1 == "C2") {
        if (cw2 == "C2h4") {
            images = {MMod(MMilnor(), 0)};
            sus = 0;
            return;
        }
        else if (cw2 == "C2h5") {
            images = {MMod(MMilnor(), 0)};
            sus = 0;
            return;
        }
        else if (cw2 == "C2h6") {
            images = {MMod(MMilnor(), 0)};
            sus = 0;
            return;
        }
    }
    else if (match_RP[0].matched) {
        int n1 = std::stoi(match_RP[1].str()), n2 = std::stoi(match_RP[2].str());
        MyException::Assert(n1 % 2 == 1 && n1 >= 1 && n1 < n2, "We need n1 % 2 == 1 && n1 >= 1 && n1 < n2");
        int n1_ = -1, n2_ = -1;
        normalize_RP(n1, n2, n1_, n2_);
        if (n1 != n1_) {
            std::string cw1_old = cw1;
            if (n2 - n1 == 1) {
                cw1 = "C2";
                n1_ = 0;
                n2_ = 1;
            }
            else 
                cw1 = fmt::format("RP{}_{}", n1_, n2_);
            fmt::print("We use {} instead of {} because of James periodicity\n", cw1, cw1_old);
        }

        std::regex_search(cw2, match_RP, words_regex);
        if (match_RP[0].matched) {
            int m1 = std::stoi(match_RP[1].str()), m2 = std::stoi(match_RP[2].str());
            MyException::Assert(n1 == m1 || n2 == m2, "n1 == m1 || n2 == m2");
            int m1_ = -1, m2_ = -1;
            normalize_RP(m1, m2, m1_, m2_);
            if (m1 != m1_) {
                std::string cw2_old = cw2;
                if (m2 - m1 == 1) {
                    cw2 = "C2";
                    m1_ = 0;
                    m2_ = 1;
                }
                else if (m2 - m1 == 0) {
                    cw2 = "S0";
                    m1_ = 0;
                    m2_ = 0;
                }
                else
                    cw2 = fmt::format("RP{}_{}", m1_, m2_);
                fmt::print("We use {} instead of {} because of James periodicity\n", cw2, cw2_old);
            }
            if (n1 == m1) {
                MyException::Assert(n2 < m2, "n2 < m2");

                int1d v_degs_n;
                Mod1d tmp, tmp1;
                Coh_RP(v_degs_n, tmp, tmp1, n1, n2, n2);
                int1d v_degs_m;
                Coh_RP(v_degs_m, tmp, tmp1, m1, m2, m2);
                images.clear();

                for (size_t i = 0; i < v_degs_n.size(); ++i)
                    images.push_back(MMod(MMilnor(), i));
                for (size_t i = v_degs_n.size(); i < v_degs_m.size(); ++i)
                    images.push_back({});
                sus = n1_ - m1_;
                return;
            }
            else if (n2 == m2) {
                MyException::Assert(n1 < m1, "n1 < m1");

                int1d v_degs_n;
                Mod1d tmp, tmp1, convert_v_n;
                Coh_RP(v_degs_n, tmp, convert_v_n, n1, n2, n2);
                int1d v_degs_m;
                Coh_RP(v_degs_m, tmp, tmp1, m1, m2, m2);
                images.clear();

                for (size_t i = 0; i < v_degs_m.size(); ++i)
                    images.push_back(convert_v_n[size_t(v_degs_m[i] - n1)]);
                sus = n2_ - m2_;
                return;
            }
        }
    }
    
    throw MyException(0x8636b4b2, "map not supported.");
}
