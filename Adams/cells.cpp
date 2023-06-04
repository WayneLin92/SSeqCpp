#include "groebner_res.h"
#include <fmt/core.h>
#include <regex>

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
 *                    tmf smash Ch
 *              Cofiber of h_n, n=0,1,2
 ***************************************************/
void Coh_tmf_Chn(int1d& v_degs, Mod1d& rels, int n, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; (1 << i) <= t_max && i <= 2; ++i)
        if (i != n)
            rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
    for (int i = 0; (1 << i) + (1 << n) <= t_max && i <= 2; ++i)
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
constexpr inline int N_MAX_RP = 261;

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

void Coh_RP(int1d& v_degs, Mod1d& rels, Mod1d& cell_reduced, int1d& min_rels, int n1, int n2, int t_max, const std::string& over)
{
    int i_max = 0;
    if (over == "S0")
        i_max = INT_MAX;
    else if (over == "tmf")
        i_max = 2;
    else {
        fmt::print("over={} not supported\n", over);
        throw MyException(0x87fe2a2b, "over not supported.");
    }

    int1d v_degs2;
    int offset = std::max(0, -n1);
    for (int n = n1; n <= n2 && n <= t_max; ++n)
        v_degs2.push_back(n + offset);

    Mod1d rels2;
    Mod tmp;
    /* Sq^k * x^n = (k, n - k) x^{n + k} when n >= 0
     * Sq^k * x^n = (k, -n - 1) x^{n + k} when n < 0
     */
    for (size_t i = 0; (1 << i) <= t_max && i <= i_max; ++i) {
        int k = 1 << i;
        for (int n = n1; n <= n2 && n + k <= t_max; ++n) {
            Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), n - n1);
            if (n + k <= n2)
                if (n >= 0 ? (k <= n && !(k & (n - k))) : !(k & (-n - 1)))
                    rel.iaddP(MMod(MMilnor(), n + k - n1), tmp);
            rels2.push_back(std::move(rel));
        }
    }

    Groebner gb(t_max, {}, v_degs2);
    min_rels.clear();
    gb.AddRels(rels2, t_max, min_rels);
    gb.MinimizeOrderedGensRels(cell_reduced, min_rels);

    v_degs = gb.v_degs();
    rels.clear();
    for (int i : min_rels)
        rels.push_back(gb.data()[i]);
}

void Coh_X2_RP(int1d& v_degs, Mod1d& rels, Mod1d& cell_reduced, int1d& min_rels, int n1, int n2, int t_max)
{
    int1d v_degs2;
    int offset = std::max(0, -n1);
    for (int n = n1; n <= n2 && n <= t_max; ++n)
        v_degs2.push_back(n + offset);

    Mod1d rels2;
    Mod tmp;
    /* Sq^k * x^n = (k, n - k) x^{n + k} when n >= 0
     * Sq^k * x^n = (k, -n - 1) x^{n + k} when n < 0
     */
    for (size_t i = 0; (1 << i) <= t_max; ++i) {
        int k = 1 << i;
        for (int n = n1; n <= n2 && n + k <= t_max; ++n) {
            Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), n - n1);
            if (n + k <= n2)
                if (n >= 0 ? (k <= n && !(k & (n - k))) : !(k & (-n - 1)))
                    rel.iaddP(MMod(MMilnor(), n + k - n1), tmp);
            rels2.push_back(std::move(rel));
        }
    }
    Groebner gb2(t_max, {}, v_degs2);
    gb2.AddRels(rels2, t_max);

    Mod1d rels3;
    for (int j = 3; j <= 5; ++j) {
        for (int i = 0; (((1 << j) - 1) << i) <= t_max; ++i) {
            int k = ((1 << j) - 1) << i;
            for (int n = n1; n <= n2 && n + k <= t_max; ++n) {
                Mod rel = MMilnor::P(i, i + j) * MMod(MMilnor(), n - n1);
                if (n + k <= n2 && !gb2.Reduce(rel + MMod(MMilnor(), n + k - n1)))
                    rel.iaddP(MMod(MMilnor(), n + k - n1), tmp);
                rels3.push_back(std::move(rel));
            }
        }
    }
    Groebner gb3(t_max, {}, v_degs2);
    min_rels.clear();
    gb3.AddRels(rels3, t_max, min_rels);
    gb3.MinimizeOrderedGensRels(cell_reduced, min_rels);

    v_degs = gb3.v_degs();
    rels.clear();
    for (int i : min_rels)
        rels.push_back(gb3.data()[i]);
}

void Coh_Fphi(int1d& v_degs, Mod1d& rels, Mod1d& cell_reduced, int1d& min_rels, int t_max)
{
    int1d v_degs2 = {0};
    for (int n = 1; n <= N_MAX_RP && n <= t_max; ++n)
        v_degs2.push_back(n + 1);

    Mod1d rels2;
    Mod tmp;
    /* Sq^k * x^n = (k, n - k) x^{n + k} */
    for (size_t i = 0; (1 << i) <= t_max; ++i) {
        size_t k = (size_t)1 << i;
        for (int n = 1; n <= N_MAX_RP && n + k <= t_max; ++n) {
            Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), n);
            if (k <= n && !(k & (n - k)) && n + k <= N_MAX_RP)
                rel.iaddP(MMod(MMilnor(), n + k), tmp);
            rels2.push_back(std::move(rel));
        }
    }
    /* Sq^k * x^{-1} = (k, 1) x^{k-1} */
    for (size_t i = 0; (1 << i) <= t_max; ++i) {
        size_t k = (size_t)1 << i;

        Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), 0);
        if (k > 1)
            rel.iaddP(MMod(MMilnor(), k - 1), tmp);
        rels2.push_back(std::move(rel));
    }

    Groebner gb(t_max, {}, v_degs2);
    min_rels.clear();
    gb.AddRels(rels2, t_max, min_rels);
    gb.MinimizeOrderedGensRels(cell_reduced, min_rels);

    v_degs = gb.v_degs();
    rels.clear();
    for (int i : min_rels)
        rels.push_back(gb.data()[i]);
}

/****************************************************
                     # Maps
 ***************************************************/
bool IsRP(const std::string& cw)
{
    std::regex is_RP_regex("^RPm{0,1}\\d+_\\d+$"); /* match example: RP1_4, RPm10_4 */
    std::smatch match_RP;
    std::regex_search(cw, match_RP, is_RP_regex);
    return match_RP[0].matched;
}

void GetRPNums(const std::string& cw, int& n1, int& n2)
{
    std::regex words_regex("^RP(m{0,1})(\\d+)_(\\d+)$"); /* match example: RP1_4, RPm10_4 */
    std::smatch match_RP;
    std::regex_search(cw, match_RP, words_regex);

    n1 = std::stoi(match_RP[2].str());
    n2 = std::stoi(match_RP[3].str());
    if (!match_RP[1].str().empty())
        n1 = -n1;
}

void SetCohMapRP2RP(const std::string& cw1, const std::string& cw2, std::string& from, std::string& to, Mod1d& images, int& sus, int& fil)
{
    int n1, n2, m1, m2;
    GetRPNums(cw1, n1, n2);
    GetRPNums(cw2, m1, m2);
    MyException::Assert(n1 <= n2 && m1 <= m2, "n1 < n2 && m1 < m2");

    /*# Normalize cw1 */

    int n1_ = 0, n2_ = 0;
    normalize_RP(n1, n2, n1_, n2_);
    if (n2 - n1 == 1) {
        from = "C2";
        n1_ = 0;
        n2_ = 1;
    }
    else if (n2 - n1 == 0) {
        from = "S0";
        n1_ = 0;
        n2_ = 0;
    }
    else
        from = fmt::format("RP{}_{}", n1_, n2_);
    if (from != cw1)
        fmt::print("We use {} instead of {} because of James periodicity\n", from, cw1);

    /*# Normalize cw2 */

    int m1_ = 0, m2_ = 0;
    normalize_RP(m1, m2, m1_, m2_);
    if (m2 - m1 == 1) {
        to = "C2";
        m1_ = 0;
        m2_ = 1;
    }
    else if (m2 - m1 == 0) {
        to = "S0";
        m1_ = 0;
        m2_ = 0;
    }
    else
        to = fmt::format("RP{}_{}", m1_, m2_);

    if (to != cw2)
        fmt::print("We use {} instead of {} because of James periodicity\n", to, cw2);

    /*# Compute map */

    if (n1 == m1) { /*## subcomplex */
        MyException::Assert(n2 < m2, "n2 < m2");

        int1d v_degs_n, tmp_ind1d;
        Mod1d tmp, tmp1;
        Coh_RP(v_degs_n, tmp, tmp1, tmp_ind1d, n1, n2, n2, "S0");
        int1d v_degs_m;
        Coh_RP(v_degs_m, tmp, tmp1, tmp_ind1d, m1, m2, m2, "S0");

        images.clear();
        for (size_t i = 0; i < v_degs_n.size(); ++i)
            images.push_back(MMod(MMilnor(), i));
        for (size_t i = v_degs_n.size(); i < v_degs_m.size(); ++i)
            images.push_back({});
        sus = n1_ - m1_;
        return;
    }
    else if (n2 == m2) { /*## quotient complex */
        MyException::Assert(n1 < m1, "n1 < m1");

        int1d v_degs_n, tmp_ind1d;
        Mod1d cell_reduced_n, tmp, tmp1;
        Coh_RP(v_degs_n, tmp, cell_reduced_n, tmp_ind1d, n1, n2, n2, "S0");
        int1d v_degs_m;
        Coh_RP(v_degs_m, tmp, tmp1, tmp_ind1d, m1, m2, m2, "S0");

        images.clear();
        for (size_t i = 0; i < v_degs_m.size(); ++i)
            images.push_back(cell_reduced_n[size_t(v_degs_m[i] - n1)]);
        sus = n2_ - m2_;
        return;
    }
    else if (n1 == m2 + 1) { /*## connecting map */
        int1d v_degs_n, tmp_ind1d;
        Mod1d cell_reduced_n, tmp, tmp1;
        Coh_RP(v_degs_n, tmp, cell_reduced_n, tmp_ind1d, n1, n2, n2, "S0");
        int1d v_degs_m;
        Mod1d rels_m;
        Coh_RP(v_degs_m, rels_m, tmp1, tmp_ind1d, m1, m2, n2, "S0");
        int1d v_degs_l;
        Mod1d rels_l;
        Coh_RP(v_degs_l, rels_l, tmp1, tmp_ind1d, m1, n2, n2, "S0");
        Groebner gb(n2, rels_l, v_degs_l);

        images.clear();
        for (size_t i = 0; i < rels_m.size(); ++i) {
            if (gb.Reduce(rels_m[i])) {
                int cell = rels_m[i].GetLead().deg_m() + v_degs_m[rels_m[i].GetLead().v()];
                images.push_back(cell_reduced_n[size_t(cell - n1)]);
            }
            else {
                images.push_back({});
            }
        }
        sus = n1_ - m2_;
        fil = 1;
        return;
    }
    else {
        fmt::print("The map between RP is not supported\n");
        throw MyException(0x31dd0ad6, "The map between RP is not supported");
    }
}

void SetCohMap(const std::string& cw1, const std::string& cw2, std::string& from, std::string& to, Mod1d& images, int& sus, int& fil)
{
    sus = 0;
    fil = 0;
    if (cw1 == "S0") {
        from = "S0";
        if (cw2 == "tmf" || cw2 == "ko" || cw2 == "X2") {
            to = cw2;
            images = {MMod(MMilnor(), 0)};
            return;
        }
    }
    if (cw2 == "S0") {
        to = "S0";
        {
            from = cw1;
            constexpr std::array<std::string_view, 7> Cs = {"C2", "Ceta", "Cnu", "Csigma", "C2h4", "C2h5", "C2h6"};
            int i = 0;
            while (i < (int)Cs.size()) {
                if (cw1 == Cs[i])
                    break;
                ++i;
            }
            if (i < (int)Cs.size()) {
                images = {MMod(MMilnor::P(i, i + 1), 0)};
                sus = 1 << i;
                return;
            }
        }
        if (cw1 == "RP1_383") {
            from = cw1;
            images = {{}};
            for (int i = 1; i <= 8; ++i) {
                images.push_back(MMod(MMilnor(), i - 1));
            }
            sus = 0;
            fil = 1;
            return;
        }
    }
    if (cw1 == "C2") {
        if (cw2 == "C2h4" || cw2 == "C2h5" || cw2 == "C2h6") {
            images = {MMod(MMilnor(), 0)};
            return;
        }
    }
    if (IsRP(cw1) && IsRP(cw2)) {
        SetCohMapRP2RP(cw1, cw2, from, to, images, sus, fil);
        return;
    }

    fmt::print("map not supported.\n");
    throw MyException(0x8636b4b2, "map not supported.");
}
