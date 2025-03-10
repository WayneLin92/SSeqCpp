#include "algebras/linalg.h"
#include "algebras/myio.h"
#include "groebner_res.h"
#include <fmt/core.h>
#include <regex>

namespace ut {
inline int get(const nlohmann::json& js, std::string key, int default_)
{
    return js.contains(key) ? js.at(key).get<int>() : default_;
}
}  // namespace ut

using json = nlohmann::json;

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
void Coh_A_over_An(int1d& v_degs, Mod1d& rels, int n, int t_max)
{
    v_degs = {0};
    rels.clear();
    for (int i = 0; i <= n && ((1 << i) <= t_max); ++i)
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
 *                     M
 ***************************************************/
void Coh_M(int1d& v_degs, Mod1d& rels, int t_max)
{
    v_degs = {0};
    rels.clear();
    Mod cell = MMod({}, 0);
    int dim = 0;
    for (int j = 0; (1 << j) - 1 <= t_max; ++j) {
        for (int i = 0; (1 << i) + dim <= t_max; ++i)
            if (i != j)
                rels.push_back(MMilnor::P(i, i + 1) * cell);
        cell = MMilnor::P(j, j + 1) * cell;
        dim += 1 << j;
    }
}

/****************************************************
 *                   j
 ***************************************************/
void Coh_j(int1d& v_degs, Mod1d& rels, int)
{
    auto& Sq = MMilnor::Sq;

    v_degs = {0, 7};
    rels.clear();
    for (int i = 0; i < 3; ++i)
        rels.push_back(MMod(MMilnor::P(i, i + 1), 0));
    Mod rel;
    rel = Sq(8) * MMod(MMilnor(), 0) + Sq(1) * MMod(MMilnor(), 1);
    rels.push_back(std::move(rel));

    rel = Sq(7) * MMod(MMilnor(), 1);
    rels.push_back(std::move(rel));

    rel = Sq(4) * (Sq(6) * MMod(MMilnor(), 1)) + Sq(6) * (Sq(4) * MMod(MMilnor(), 1));
    rels.push_back(std::move(rel));
}

/****************************************************
 *                   j/2
 ***************************************************/
void Coh_j_C2(int1d& v_degs, Mod1d& rels, int)
{
    auto& Sq = MMilnor::Sq;

    v_degs = {0, 7};
    rels.clear();

    Mod i7Sq1 = Sq(8) * MMod(MMilnor(), 0) + Sq(1) * MMod(MMilnor(), 1);

    Mod rel;
    rel = Sq(2) * MMod(MMilnor(), 0);
    rels.push_back(std::move(rel));

    rel = Sq(2) * (Sq(1) * MMod(MMilnor(), 0));
    rels.push_back(std::move(rel));

    rel = Sq(4) * MMod(MMilnor(), 0);
    rels.push_back(std::move(rel));

    rel = Sq(7) * MMod(MMilnor(), 1) + Sq(6) * i7Sq1;
    rels.push_back(std::move(rel));

    rel = Sq(4) * (Sq(6) * MMod(MMilnor(), 1)) + Sq(6) * (Sq(4) * MMod(MMilnor(), 1)) + Sq(3) * (Sq(6) * i7Sq1) + Sq(6) * (Sq(3) * i7Sq1) + Sq(4) * (Sq(5) * i7Sq1) + Sq(5) * (Sq(4) * i7Sq1);
    rels.push_back(std::move(rel));
}

/****************************************************
 *                   RPn
 ***************************************************/
constexpr inline int N_MAX_RP = 261;

int Period_RP(int n)
{
    ErrorIdMsg::Assert(n >= 0, "n >= 0 in phi(n)");
    int phi = (n / 8) * 4, residue = n % 8;
    if (residue <= 1)
        phi += residue;
    else if (residue <= 3)
        phi += 2;
    else if (residue <= 7)
        phi += 3;
    return phi <= 31 ? 1 << phi : -1;
}

void normalize_RP(int n1, int n2, int& n1_, int& n2_)  //// TODO: normalize CP and HP
{
    n1_ = n1;
    n2_ = n2;
    int T = Period_RP(n2 - n1);
    if (T >= 0) {
        if (n1 < 0 && T < -n1 * 2) {
            int shift = ((-n1) / T + 1) * T;
            n1 += shift;
            n2 += shift;
        }
        if (n1 > 0) {
            n1_ = 1 + ((n1 - 1) % T);
            n2_ = n2 - (n1 - n1_);
        }
    }
}

void Coh_P(int1d& v_degs, Mod1d& rels, Mod1d& cells, int1d& min_rels, int n1, int n2, int t_max, const std::string& over, int hopf)
{
    size_t i_max = 0;
    if (over == "S0")
        i_max = 1ULL << 30;
    else if (over == "tmf")
        i_max = 2;
    else {
        fmt::print("over={} not supported\n", over);
        throw ErrorIdMsg(0x87fe2a2b, "over not supported.");
    }

    int1d v_degs2;
    int offset = std::max(0, -n1);
    for (int n = n1; n <= n2 && n <= t_max; ++n)
        v_degs2.push_back((n + offset) << hopf);

    Mod1d rels2;
    Mod tmp;
    /* Sq^{k} * x^n = (k1, n - k1) x^{n + k1} when n >= 0 and k = (k1 << hopf)
     * Sq^{k} * x^n = (k1, -n - 1) x^{n + k1} when n < 0 and k = (k1 << hopf)
     */
    for (size_t i = 0; (1 << i) <= t_max && i <= i_max; ++i) {
        int k = 1 << i;
        int k1 = k >> hopf;
        for (int n = n1; n <= n2 && ((n << hopf) + k) <= t_max; ++n) {
            Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), n - n1);
            if (k1 && n + k1 <= n2)
                if (n >= 0 ? (k1 <= n && !(k1 & (n - k1))) : !(k1 & (-n - 1)))
                    rel.iaddP(MMod(MMilnor(), n + k1 - n1), tmp);
            rels2.push_back(std::move(rel));
        }
    }

    Groebner gb(t_max, {}, v_degs2);
    min_rels.clear();
    gb.AddRels(rels2, t_max, min_rels);
    gb.MinimizeOrderedGensRels(cells, min_rels);

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

void Coh_Fphi_n(int1d& v_degs, Mod1d& rels, Mod1d& cell_reduced, int1d& min_rels, int n_max, int t_max)
{
    int1d v_degs2 = {0};
    for (int n = 1; n <= n_max && n <= t_max; ++n)
        v_degs2.push_back(n + 1);

    Mod1d rels2;
    Mod tmp;
    /* Sq^k * x^n = (k, n - k) x^{n + k} */
    for (size_t i = 0; (1 << i) <= t_max; ++i) {
        size_t k = (size_t)1 << i;
        for (int n = 1; n <= n_max && n + k <= t_max; ++n) {
            Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), n);
            if (k <= n && !(k & (n - k)) && n + k <= n_max)
                rel.iaddP(MMod(MMilnor(), n + k), tmp);
            rels2.push_back(std::move(rel));
        }
    }
    /* Sq^k * x^{-1} = (k, 1) x^{k-1} */
    for (size_t i = 0; (1 << i) <= t_max; ++i) {
        size_t k = (size_t)1 << i;

        Mod rel = MMilnor::Sq((uint32_t)k) * MMod(MMilnor(), 0);
        if (k > 1 && k - 1 <= n_max)
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

struct Cohomology
{
    int1d v_degs;
    Mod1d rels;
    Mod1d cells; /* minimal presentation of the i'th cell */
    int1d min_rels;
    std::map<int, int> num_cells;     /* num_cells[d] is the number of d-cells */
    std::map<int, int> indices_cells; /* indices_cells[d] is the first index of d-cell */
};

inline std::string Serialize(const int1d& arr)
{
    return myio::StrCont("", ",", "", "", arr, [](int i) { return std::to_string(i); });
}

int GetCohFromJson(const json& js, const std::string& name_, int t_max, Cohomology& coh)
{
    try {
        auto& cws = js.at("CW_complexes");
        std::smatch match;
        std::string ring = "S0";
        std::string name = name_;
        if (std::regex_search(name_, match, std::regex("^tmf_(\\w+)$")); match[0].matched) { /* tmf_C2 */
            name = match[1].str();
            ring = "tmf";
        }
        if (!cws.contains(name)) {
            return -1;
        }
        auto& cw_json = cws.at(name);
        if (!cw_json.contains("operations")) {
            fmt::print("json missing key: operations\n");
            return -2;
        }
        int1d gen_degs;
        int index = 0;

        /* Set coh.num_cells and coh.indices_cells */
        for (auto& c : cw_json.at("cells")) {
            int d = c.get<int>();
            gen_degs.push_back(d);
            ++coh.num_cells[d];
            if (!ut::has(coh.indices_cells, d))
                coh.indices_cells[d] = index;
            ++index;
        }
        int1d cells_gen;
        for (auto& c : cw_json.at("cells_gen"))
            cells_gen.push_back(c.get<int>());
        std::map<int, std::map<int, int1d>> ops;
        for (auto& op : cw_json.at("operations")) {
            int c0, i0;
            if (op[0].is_number()) {
                c0 = op[0].get<int>();
                i0 = coh.indices_cells.at(c0);
            }
            else {
                int1d arr = op[0].get<std::vector<int>>();
                c0 = arr[0];
                ErrorIdMsg::Assert(arr[1] < coh.num_cells.at(c0), "arr[1] < cells.at(c0)");
                i0 = coh.indices_cells.at(c0) + arr[1];
            }
            int c1;
            int1d i1s;
            if (op[1].is_number()) {
                c1 = op[1].get<int>();
                i1s.push_back(coh.indices_cells.at(c1));
            }
            else {
                int1d arr = op[1].get<std::vector<int>>();
                c1 = arr[0];
                for (size_t i = 1; i < arr.size(); ++i) {
                    ErrorIdMsg::Assert(arr[i] < coh.num_cells.at(c1), "arr[i] < cells.at(c1)");
                    i1s.push_back(coh.indices_cells.at(c1) + arr[i]);
                }
            }

            int n = c1 - c0;
            ErrorIdMsg::Assert(!(n & (n - 1)), "n is a power of 2");
            for (int i1 : i1s)
                ops[i0][n].push_back(i1);
        }

        Mod1d rels2;
        Mod tmp;
        if (ring == "S0") {
            for (size_t i = 0; i < gen_degs.size(); ++i) {
                for (int j = 0; (1 << j) + gen_degs[i] <= t_max; ++j) {
                    Mod rel = MMilnor::P(j, j + 1) * MMod(MMilnor(), i);
                    if (ut::has(ops, (int)i) && ut::has(ops.at((int)i), (int)(1 << j)))
                        for (int v1 : ops.at((int)i).at(int(1 << j)))
                            rel.iaddP(MMod(MMilnor(), v1), tmp);
                    rels2.push_back(std::move(rel));
                }
            }
        }
        else if (ring == "tmf") {
            for (size_t i = 0; i < gen_degs.size(); ++i) {
                for (int j = 0; j <= 2 && (1 << j) + gen_degs[i] <= t_max; ++j) {
                    Mod rel = MMilnor::P(j, j + 1) * MMod(MMilnor(), i);
                    if (ut::has(ops, (int)i) && ut::has(ops.at((int)i), (int)(1 << j)))
                        for (int v1 : ops.at((int)i).at(int(1 << j)))
                            rel.iaddP(MMod(MMilnor(), v1), tmp);
                    rels2.push_back(std::move(rel));
                }
            }
        }

        Groebner gb(t_max, {}, gen_degs);
        gb.AddRels(rels2, t_max, coh.min_rels);
        gb.MinimizeOrderedGensRels(coh.cells, coh.min_rels);

        int1d cells_gen_v2;
        for (size_t i = 0; i < coh.cells.size(); ++i)
            if (coh.cells[i].data.size() == 1 && coh.cells[i].GetLead().deg_m() == 0)
                cells_gen_v2.push_back(gen_degs[i]);
        ut::RemoveIf(cells_gen, [t_max](int n) { return n > t_max; });
        ut::RemoveIf(cells_gen_v2, [t_max](int n) { return n > t_max; });
        ErrorIdMsg::Assert(cells_gen == cells_gen_v2, fmt::format("cells_gen == cells_gen_v2 = {}", Serialize(cells_gen_v2)));

        coh.v_degs = gb.v_degs();
        ut::RemoveIf(coh.v_degs, [t_max](int i) { return i > t_max; });
        for (int i : coh.min_rels)
            coh.rels.push_back(gb.data()[i]);
    }
    catch (nlohmann::detail::exception& e) {
        fmt::print("Error({:#x}) - {}\n", e.id, e.what());
        throw e;
    }
    return 0;
}

int GetCohFromJson(const std::string& name, int t_max, int1d& v_degs, Mod1d& rels)
{
    Cohomology coh;
    auto js = myio::load_json("Adams.json");
    if (int error = GetCohFromJson(js, name, t_max, coh))
        return error;
    v_degs = std::move(coh.v_degs);
    rels = std::move(coh.rels);
    return 0;
}

Mod Combination(Mod1d& xs, uint32_t n)
{
    Mod x;
    for (int i : ut::two_exp(n))
        if (i < xs.size())
            x += xs[i];
    return x;
}

/****************************************************
                     # d2
 ***************************************************/

int GetD2FromJson(const std::string& name, int t_max, Mod1d& h_d2_images, const int1d& /*nTry*/)
{
    if (name == "j") {
        auto nTry1 = int1d{1, 1, 0, 2, 0, 5};
        ut::get(h_d2_images, 3) = MMod({}, 1); /* stem=7 */

        Mod1d xs8 = {Mod(MMilnor::Sqm({8}), 0)};
        Mod1d xs14 = {Mod(MMilnor::Sqm({0, 0, 1}), 1), Mod(MMilnor::Sqm({0, 0, 2}), 0)};
        Mod1d xs16 = {Mod(MMilnor::Sqm({0, 3}), 1), Mod(MMilnor::Sqm({2, 0, 1}), 1), Mod(MMilnor::Sqm({16}), 0)};
        Mod1d xs17 = {Mod(MMilnor::Sqm({10}), 1), Mod(MMilnor::Sqm({0, 1, 1}), 1)};
        Mod1d xs19 = {Mod(MMilnor::Sqm({12}), 1), Mod(MMilnor::Sqm({0, 4}), 1), Mod(MMilnor::Sqm({2, 1, 1}), 1)};
        Mod1d xs20 = {Mod(MMilnor::Sqm({10, 1}), 1), Mod(MMilnor::Sqm({0, 2, 1}), 1), Mod(MMilnor::Sqm({8, 4}), 0)};
        Mod1d xs22 = {Mod(MMilnor::Sqm({2, 2, 1}), 1), Mod(MMilnor::Sqm({0, 5}), 1), Mod(MMilnor::Sqm({8, 0, 1}), 1), Mod(MMilnor::Sqm({0, 0, 0, 1}), 1), Mod(MMilnor::Sqm({8, 0, 2}), 0)};

        ut::get(h_d2_images, 4) = Combination(xs8, ut::get(nTry1, 0, 0));
        ut::get(h_d2_images, 5) = Combination(xs14, ut::get(nTry1, 1, 0));
        ut::get(h_d2_images, 6) = Combination(xs16, ut::get(nTry1, 2, 0));
        ut::get(h_d2_images, 7) = Combination(xs17, ut::get(nTry1, 3, 0));
        ut::get(h_d2_images, 8) = Combination(xs19, ut::get(nTry1, 4, 0));
        ut::get(h_d2_images, 9) = Combination(xs20, ut::get(nTry1, 5, 0));
        ut::get(h_d2_images, 10) = Combination(xs22, ut::get(nTry1, 6, 0));
        return 0;
    }
    Cohomology coh;
    if (!myio::FileExists("Adams.json"))
        return 0;
    auto js = myio::load_json("Adams.json");
    if (int error = GetCohFromJson(js, name, t_max, coh))
        return error;
    if (!js.at("CW_complexes").contains(name))
        return 0;
    if (!js.at("CW_complexes").at(name).contains("secondary"))
        return 0;
    auto& secondary = js.at("CW_complexes").at(name).at("secondary");
    for (auto it = secondary.begin(); it != secondary.end(); ++it) {
        int i = std::stoi(it.key());
        int c = it.value().get<int>();
        int iCell = coh.indices_cells.at(c);
        ut::get(h_d2_images, i) = coh.cells[iCell];
    }
    return 0;
}

/****************************************************
                     # Maps
 ***************************************************/
bool IsFP(const std::string& cw, const std::string& field)
{
    std::regex is_RP_regex("^" + field + "P(?:m|)\\d+_(?:m|)\\d+$"); /* match example: RP1_4, RPm10_m4 */
    std::smatch match_RP;
    std::regex_search(cw, match_RP, is_RP_regex);
    return match_RP[0].matched;
}

bool IsFP(const std::string& cw)
{
    std::regex is_RP_regex("^(?:R|C|H)P(?:m|)\\d+_(?:m|)\\d+$"); /* match example: RP1_4, CPm10_m4 */
    std::smatch match_RP;
    std::regex_search(cw, match_RP, is_RP_regex);
    return match_RP[0].matched;
}

int IsFphi_n(const std::string& cw)
{
    std::regex regex("^Fphi(\\d+)$"); /* match example: Fphi4 */
    std::smatch match;
    std::regex_search(cw, match, regex);
    return match[0].matched ? std::stoi(match[1].str()) : 0;
}

// TODO: write myio::parse to parse a string to assign value to variables
void ParseFP(const std::string& cw, int& hopf, int& n1, int& n2)
{
    std::regex words_regex("^(R|C|H)P(m|)(\\d+)_(m|)(\\d+)$"); /* match example: CP1_4 */
    std::smatch match_P;
    std::regex_search(cw, match_P, words_regex);

    std::string field = match_P[1].str();
    hopf = field == "R" ? 0 : (field == "C" ? 1 : 2);
    n1 = std::stoi(match_P[3].str());
    n2 = std::stoi(match_P[5].str());
    if (!match_P[2].str().empty())
        n1 = -n1;
    if (!match_P[4].str().empty())
        n2 = -n2;
}

void GetCohMapFP2FP(const std::string& cw1, const std::string& cw2, std::string& from, std::string& to, Mod1d& images, int& sus, int& fil)
{
    int hopf1, hopf2, n1, n2, m1, m2;
    ParseFP(cw1, hopf1, n1, n2);
    ParseFP(cw2, hopf2, m1, m2);
    ErrorIdMsg::Assert(n1 <= n2 && m1 <= m2, "n1 < n2 && m1 < m2");
    ErrorIdMsg::Assert(hopf1 == hopf2 || hopf1 + 1 == hopf2, "hopf1 == hopf2 || hopf1 + 1 == hopf2");
    int n1_ = n1, n2_ = n2;
    int m1_ = m1, m2_ = m2;
    std::string field1 = hopf1 == 0 ? "R" : (hopf1 == 1 ? "C" : "H");
    std::string field2 = hopf2 == 0 ? "R" : (hopf2 == 1 ? "C" : "H");

    /*# Normalize cw1 */
    if (hopf1 == 0)
        normalize_RP(n1, n2, n1_, n2_);
    if (n2 - n1 == 1) {
        from = hopf1 == 0 ? "C2" : (hopf1 == 1 ? "Ceta" : "Cnu");
        n1_ = 0;
        n2_ = 1;
    }
    else if (n2 - n1 == 0) {
        from = "S0";
        n1_ = 0;
        n2_ = 0;
    }
    else if (n1 == 3 && n2 == 5) {
        from = hopf1 == 0 ? "CW_2_V_eta" : (hopf1 == 1 ? "CW_eta_V_nu" : "CW_nu_V_sigma");
        n1_ = 0;
        n2_ = 2;
    }
    else if (n1 == 2 && n2 == 4) {
        from = hopf1 == 0 ? "CW_2_A_eta" : (hopf2 == 1 ? "CW_eta_A_nu" : "CW_nu_A_sigma");
        n1_ = 0;
        n2_ = 2;
    }
    else
        from = fmt::format("{}P{}_{}", field1, n1_, n2_);
    if (from != cw1)
        fmt::print("We use {} instead of {}\n", from, cw1);

    /*# Normalize cw2 */
    if (hopf2 == 0)
        normalize_RP(m1, m2, m1_, m2_);
    if (m2 - m1 == 1) {
        to = hopf2 == 0 ? "C2" : (hopf2 == 1 ? "Ceta" : "Cnu");
        m1_ = 0;
        m2_ = 1;
    }
    else if (m2 - m1 == 0) {
        to = "S0";
        m1_ = 0;
        m2_ = 0;
    }
    else if (m1 == 3 && m2 == 5) {
        to = hopf2 == 0 ? "CW_2_V_eta" : (hopf2 == 1 ? "CW_eta_V_nu" : "CW_nu_V_sigma");
        m1_ = 0;
        m2_ = 2;
    }
    else if (m1 == 2 && m2 == 4) {
        to = hopf2 == 0 ? "CW_2_A_eta" : (hopf2 == 1 ? "CW_eta_A_nu" : "CW_nu_A_sigma");
        m1_ = 0;
        m2_ = 2;
    }
    else
        to = fmt::format("{}P{}_{}", field2, m1_, m2_);
    if (to != cw2)
        fmt::print("We use {} instead of {}\n", to, cw2);

    if (hopf1 == hopf2) {
        int hopf = hopf1;
        /*# Compute map */

        if (n1 == m1) { /*## subcomplex */
            ErrorIdMsg::Assert(n2 < m2, "n2 < m2");

            int1d v_degs_n, tmp_ind1d;
            Mod1d tmp, tmp1;
            Coh_P(v_degs_n, tmp, tmp1, tmp_ind1d, n1, n2, std::min(n2 << hopf, 256), "S0", hopf);
            int1d v_degs_m;
            Coh_P(v_degs_m, tmp, tmp1, tmp_ind1d, m1, m2, std::min(m2 << hopf, 256), "S0", hopf);

            images.clear();
            for (size_t i = 0; i < v_degs_n.size(); ++i)
                images.push_back(MMod(MMilnor(), i));
            for (size_t i = v_degs_n.size(); i < v_degs_m.size(); ++i)
                images.push_back({});
            sus = (n1_ - m1_) << hopf;
            return;
        }
        else if (n2 == m2) { /*## quotient complex */
            ErrorIdMsg::Assert(n1 < m1, "n1 < m1");

            int1d v_degs_n, tmp_ind1d;
            Mod1d cell_reduced_n, tmp, tmp1;
            Coh_P(v_degs_n, tmp, cell_reduced_n, tmp_ind1d, n1, n2, std::min(n2 << hopf, 256), "S0", hopf);
            int1d v_degs_m;
            Coh_P(v_degs_m, tmp, tmp1, tmp_ind1d, m1, m2, std::min(m2 << hopf, 256), "S0", hopf);

            images.clear();
            for (size_t i = 0; i < v_degs_m.size(); ++i)
                images.push_back(cell_reduced_n[size_t((v_degs_m[i] << hopf) - n1)]);
            sus = (n2_ - m2_) << hopf;
            return;
        }
        else if (n1 == m2 + 1) { /*## attaching map */
            int1d v_degs_n, tmp_ind1d;
            Mod1d cell_reduced_n, tmp, tmp1;
            Coh_P(v_degs_n, tmp, cell_reduced_n, tmp_ind1d, n1, n2, std::min(n2 << hopf, 256), "S0", hopf);
            int1d v_degs_m;
            Mod1d rels_m;
            Coh_P(v_degs_m, rels_m, tmp1, tmp_ind1d, m1, m2, std::min(n2 << hopf, 256), "S0", hopf);
            int1d v_degs_l;
            Mod1d rels_l;
            Coh_P(v_degs_l, rels_l, tmp1, tmp_ind1d, m1, n2, std::min(n2 << hopf, 256), "S0", hopf);
            Groebner gb(n2, rels_l, v_degs_l);

            images.clear();
            for (size_t i = 0; i < rels_m.size(); ++i) {
                if (gb.Reduce(rels_m[i])) {
                    int cell = ((rels_m[i].GetLead().deg_m() + v_degs_m[rels_m[i].GetLead().v()]) << hopf) + std::min(m1, 0);
                    images.push_back(cell_reduced_n[size_t(cell - n1)]);
                }
                else {
                    images.push_back({});
                }
            }
            sus = (n1_ - m2_) << hopf;
            fil = 1;
            return;
        }
        else {
            fmt::print("The map between FP is not supported\n");
            throw ErrorIdMsg(0x31dd0ad6, "The map between FP is not supported");
        }
    }
    else if (hopf2 == hopf1 + 1) {
        ErrorIdMsg::Assert(n1 == 2 * m1 - 1 && n2 == 2 * m2, "n1 == 2 * m1 - 1 && n2 == 2 * m2");

        int1d v_degs_n, tmp_ind1d;
        Mod1d cell_reduced_n, tmp, tmp1;
        Coh_P(v_degs_n, tmp, cell_reduced_n, tmp_ind1d, n1, n2, std::min(n2 << hopf1, 256), "S0", hopf1);
        int1d v_degs_m;
        Coh_P(v_degs_m, tmp, tmp1, tmp_ind1d, m1, m2, std::min(m2 << hopf2, 256), "S0", hopf2);

        images.clear();
        for (size_t i = 0; i < v_degs_m.size(); ++i)
            images.push_back(cell_reduced_n[size_t((v_degs_m[i] - n1) >> hopf1)]);
        sus = (n2_ - 2 * m2_) << hopf1;
        return;
    }
}

int GetCoh(int1d& v_degs, Mod1d& rels, int d_max, const std::string& cw)
{
    std::regex is_P_regex("^(tmf_|X2_|)((R|C|H)P(?:m|)\\d+_(?:m|)\\d+)$"); /* match example: RP1_4, CPm1_10 */
    std::smatch match;
    if (cw == "S0")
        Coh_S0(v_degs, rels, d_max);
    else if (cw == "X2")
        Coh_X2(v_degs, rels, d_max);
    else if (cw == "tmf")
        Coh_A_over_An(v_degs, rels, 2, d_max);
    else if (cw == "A_A3")
        Coh_A_over_An(v_degs, rels, 3, d_max);
    else if (cw == "A_A4")
        Coh_A_over_An(v_degs, rels, 4, d_max);
    else if (cw == "A_A5")
        Coh_A_over_An(v_degs, rels, 5, d_max);
    else if (cw == "ko")
        Coh_ko(v_degs, rels, d_max);
    else if (cw == "j")
        Coh_j(v_degs, rels, d_max);
    else if (cw == "j_C2")
        Coh_j_C2(v_degs, rels, d_max);
    else if (cw == "Fphi") {
        Mod1d tmp;
        int1d tmp_int1d;
        Coh_Fphi(v_degs, rels, tmp, tmp_int1d, d_max);
    }
    else if (std::regex_search(cw, match, std::regex("^Fphi(\\d+)$")); match[0].matched) { /* Fphin */
        Mod1d tmp;
        int1d tmp_int1d;
        int n = std::stoi(match[1].str());
        Coh_Fphi_n(v_degs, rels, tmp, tmp_int1d, n, d_max);
    }
    else if (std::regex_search(cw, match, is_P_regex); match[0].matched) {
        std::string ring = match[1].str();
        int hopf, n1, n2;
        ParseFP(match[2].str(), hopf, n1, n2);
        std::string field = match[3].str();
        ErrorIdMsg::Assert(n1 + 1 < n2, "n1 + 1 < n2");
        if (field == "R") {
            int n1_ = n1, n2_ = n2;
            normalize_RP(n1, n2, n1_, n2_);
            if (n1 != n1_) {
                fmt::print("Please use n1={} n2={} because of James periodicity\n", n1_, n2_);
                return -1;
            }
        }
        if (n1 == 2 && n2 == 4) {
            fmt::print("Please use CW_2_A_eta instead\n");
            return -2;
        }
        if (n1 == 3 && n2 == 5) {
            fmt::print("Please use CW_2_V_eta instead\n");
            return -3;
        }

        Mod1d tmp_Mod1d;
        int1d tmp_ind1d;
        if (ring == "")
            Coh_P(v_degs, rels, tmp_Mod1d, tmp_ind1d, n1, n2, d_max, "S0", hopf);
        else if (ring == "tmf_")
            Coh_P(v_degs, rels, tmp_Mod1d, tmp_ind1d, n1, n2, d_max, "tmf", hopf);
        else if (ring == "X2_")
            Coh_X2_RP(v_degs, rels, tmp_Mod1d, tmp_ind1d, n1, n2, d_max);  ////
        else {
            fmt::print("ring={} not supported\n", ring);
            return -4;
        }
    }
    else if (cw == "M") {
        Coh_M(v_degs, rels, d_max);
    }
    else if (int error = GetCohFromJson(cw, d_max, v_degs, rels)) {
        fmt::print("Error({}) - Unsupported arugment cw={}\n", error, cw);
        return -5;
    }
    return 0;
}

Mod cell2Mod(const Cohomology& coh, const json& c)
{
    Mod result;
    if (c.is_number()) {
        int i = coh.indices_cells.at(c.get<int>());
        result = coh.cells[i];
    }
    else {
        int1d arr = c.get<std::vector<int>>();
        if (!arr.empty()) {
            int d = arr[0];
            for (size_t i = 1; i < arr.size(); ++i)
                result += coh.cells[size_t(coh.indices_cells.at(d) + arr[i])];
        }
    }
    return result;
}

constexpr int MAP_T_MAX = 265;

int GetCohMapFromJson(const std::string& name_, std::string& from_, std::string& to_, Mod1d& images, int& sus, int& fil)
{
    auto js = myio::load_json("Adams.json");
    const int t_max = MAP_T_MAX;
    try {
        auto& cws = js.at("CW_complexes");
        auto& maps = js.at("maps");

        std::smatch match;
        std::string ring = "";
        std::string name = name_;
        if (std::regex_search(name_, match, std::regex("^tmf(?:_(\\w+)|)__tmf(?:_(\\w+)|)$")); match[0].matched) { /* tmf_C2__tmf */
            auto cw1 = match[1].str();
            auto cw2 = match[2].str();
            if (cw1.empty())
                cw1 = "S0";
            if (cw2.empty())
                cw2 = "S0";
            ring = "tmf_";
            name = fmt::format("{}__{}", cw1, cw2);
        }

        if (!maps.contains(name))
            return -1;
        auto& map_json = maps.at(name);
        auto from = map_json.at("from").get<std::string>();
        auto to = map_json.at("to").get<std::string>();
        from_ = from == "S0" ? ring : fmt::format("{}{}", ring, from);
        to_ = to == "S0" ? ring : fmt::format("{}{}", ring, to);
        if (to_.empty())
            to_ = "S0";
        sus = myio::get(map_json, "sus", 0);
        images = {};
        if (map_json.contains("images")) {
            fil = 0;
            int1d cells_gen_to;
            for (auto& c : cws.at(to).at("cells_gen"))
                cells_gen_to.push_back(c.get<int>());
            Cohomology coh_from;
            GetCohFromJson(js, from, t_max, coh_from);

            auto& images_json = map_json.at("images");
            if (images_json.size() != cells_gen_to.size())
                return -2;
            for (size_t i = 0; i < images_json.size(); ++i)
                images.push_back(cell2Mod(coh_from, images_json[i]));

            /* Verify the map */
            Cohomology coh_to;
            GetCohFromJson(js, to, t_max, coh_to);
            Groebner gb_from(t_max, {}, coh_from.v_degs);
            gb_from.AddRels(coh_from.rels, t_max);
            for (const auto& rel : coh_to.rels) {
                int deg_f_rel = rel.GetLead().deg_m() + coh_to.v_degs[rel.GetLead().v()] + sus;
                if (deg_f_rel <= MAP_T_MAX && gb_from.Reduce(subs(rel, images))) {
                    fmt::print("the map is not an A-module homomorphism");
                    return 2;
                }
            }
        }
        else if (map_json.contains("total")) {
            fil = 1;
            std::string total = map_json.at("total");

            Cohomology coh_to;
            GetCohFromJson(js, to, t_max, coh_to);
            Groebner gb_to(t_max, {}, coh_to.v_degs);
            gb_to.AddRels(coh_to.rels, t_max);

            Cohomology coh_from;
            GetCohFromJson(js, from, t_max, coh_from);

            Cohomology coh_total;
            GetCohFromJson(js, total, t_max, coh_total);
            Groebner gb(t_max, {}, coh_total.v_degs);
            gb.AddRels(coh_total.rels, t_max);

            std::string tmp1, tmp2;
            int sus_q, fil_q;
            Mod1d images_q;
            GetCohMapFromJson(fmt::format("{}__{}", total, from), tmp1, tmp2, images_q, sus_q, fil_q);

            int1d cells_from;
            for (auto& c : cws.at(from).at("cells"))
                cells_from.push_back(c.get<int>());

            Mod1d sub_lift;  // TODO: Compute this without specifying it in json
            if (map_json.contains("sub_lift")) {
                for (auto& i : map_json.at("sub_lift")) {
                    if (i.is_array()) {
                        if (i.size())
                            return -9;
                        sub_lift.push_back({});
                    }
                    else
                        sub_lift.push_back(steenrod::MMod({}, i.get<int>()));
                }
            }

            for (int i : coh_to.min_rels) {
                steenrod::Mod rel0 = gb_to.data()[i];
                Mod rel = sub_lift.empty() ? rel0 : subs(rel0, sub_lift);
                rel = gb.Reduce(std::move(rel));
                if (rel) {
                    int c = rel.GetLead().deg_m() + coh_total.v_degs[rel.GetLead().v()] - sus_q;
                    int j_max = 1 << coh_from.num_cells[c];
                    Mod x;
                    int j = 1;
                    for (; j < j_max; ++j) {  // TODO: Use linear algebra
                        x.data.clear();
                        for (int je : ut::two_exp((uint32_t)j))
                            x += coh_from.cells[coh_from.indices_cells[c] + je];
                        auto y = subs(x, images_q);
                        y = gb.Reduce(std::move(y));
                        if (y == rel) {
                            images.push_back(x);
                            break;
                        }
                    }
                    if (j == j_max) {
                        return -6;
                    }
                }
                else
                    images.push_back({});
            }
        }
        else
            return -8;
    }
    catch (nlohmann::detail::exception& e) {
        fmt::print("json error: {:#x} - {}\n", e.id, e.what());
        throw e;
    }
    /*catch (MyException&) {
    }*/
    return 0;
}

void SetCohMap(const std::string& cw1, const std::string& cw2, std::string& from, std::string& to, Mod1d& images, int& sus, int& fil)
{
    from = cw1;
    to = cw2;
    sus = 0;
    fil = 0;
    if (cw1 == "S0") {
        if (cw2 == "tmf" || cw2 == "ko" || cw2 == "X2" || cw2 == "A_A3" || cw2 == "A_A4" || cw2 == "A_A5") {
            images = {MMod(MMilnor(), 0)};
            return;
        }
    }
    if (cw2 == "S0") {
        if (cw1 == "RP1_256") {
            images = {{}};
            for (int i = 1; i <= 8; ++i)
                images.push_back(MMod(MMilnor(), i - 1));
            fil = 1;
            return;
        }
        if (cw1 == "CP1_128") {
            images = {{}, {}};
            for (int i = 2; i <= 8; ++i)
                images.push_back(MMod(MMilnor(), i - 2));
            fil = 1;
            sus = -1;
            return;
        }
        if (cw1 == "HP1_64") {
            images = {{}, {}, {}};
            for (int i = 3; i <= 8; ++i)
                images.push_back(MMod(MMilnor(), i - 3));
            fil = 1;
            sus = -3;
            return;
        }
    }
    if (cw1 == "RP1_256" && cw2 == "tmf_RP1_256") {
        images = {};
        int1d v_degs;
        int1d v_degs_n, tmp_ind1d;
        const int t_max = 256;
        Mod1d cell_reduced, tmp;
        Coh_P(v_degs, tmp, cell_reduced, tmp_ind1d, 1, t_max, t_max, "S0", 0);
        for (int i = 0; i < 2; ++i)
            images.push_back(MMod({}, i));
        for (int i = 2; 7 + 8 * (i - 2) <= t_max; ++i)
            images.push_back(cell_reduced[size_t(6 + 8 * (i - 2))]);
        return;
    }
    if (cw2 == fmt::format("tmf_{}", cw1)) {
        images = {MMod(MMilnor(), 0)};
        return;
    }
    if (cw1 == "CW_sigma_nu_eta_2" && cw2 == "tmf") {
        images = {MMod(MMilnor(), 0)};
        return;
    }
    if (cw1 == "C2") {
        if (cw2 == "j_C2") {
            images = {MMod(MMilnor(), 0), {}};
            return;
        }
    }
    if (cw1 == "j" && cw2 == "j_C2") {
        images = {MMod(MMilnor(), 0), MMod(MMilnor(), 1)};
        return;
    }
    if (cw1 == "j_C2" && cw2 == "j") {
        Mod i7Sq1 = MMilnor::Sq(8) * MMod(MMilnor(), 0) + MMilnor::Sq(1) * MMod(MMilnor(), 1);
        images = {MMod(MMilnor::P(0, 1), 0), i7Sq1};
        sus = 1;
        return;
    }

    if (cw1 == "CW_eta_2" && cw2 == "Fphi4") {
        images = {MMod(MMilnor(), 0)};
        return;
    }
    if (IsFP(cw1) && IsFP(cw2)) {
        GetCohMapFP2FP(cw1, cw2, from, to, images, sus, fil);
        return;
    }
    if (int n1 = IsFphi_n(cw1), n2 = IsFphi_n(cw2); 0 < n1 && n1 < n2) {
        images = {MMod(MMilnor(), 0)};
        return;
    }
    if (cw1 == "Fphi" && cw2 == "RP1_256") {
        images = {};
        for (int i = 1; i <= 8; ++i)
            images.push_back(MMod(MMilnor::P(i, i + 1), 0));
        sus = 1;
        return;
    }
    if (int n = IsFphi_n(cw1)) {
        if (cw2 == fmt::format("RP1_{}", n)) {
            images = {};
            for (int i = 1; (1 << i) - 1 <= n; ++i)
                images.push_back(MMod(MMilnor::P(i, i + 1), 0));
            sus = 1;
            return;
        }
    }
    std::string name = fmt::format("{}__{}", cw1, cw2);
    if (int error = GetCohMapFromJson(name, from, to, images, sus, fil)) {
        fmt::print("Error({}) - map not supported.\n", error);
        throw ErrorIdMsg(0x8636b4b2, "map not supported.");
    }
}

int1d SteenrodMod2Indices(Mod x, const MMod1d& basis)
{
    int1d result;
    Mod tmp;
    size_t j = 0;
    for (int i = 0; i < (int)basis.size() && j < x.data.size(); ++i) {
        if (x.data[j] == basis[i]) {
            result.push_back(i);
            ++j;
        }
    }
#ifndef NDEBUG
    if (j < (int)x.data.size()) {
        fmt::print("MyException(0xfd3d2814): Index not found\n");
        throw ErrorIdMsg(0xfd3d2814, "Index not found");
    }
#endif
    return result;
}

int1d SteenrodMod2Indices(Mod x, const MMod1d& basis, const int2d& inv)
{
    int1d result;
    Mod tmp;
    size_t j = 0;
    for (int i = 0; i < (int)basis.size() && j < x.data.size(); ++i) {
        if (x.data[j] == basis[i]) {
            result = lina::add(result, inv[i]);
            ++j;
        }
    }
#ifndef NDEBUG
    if (j < (int)x.data.size()) {
        fmt::print("MyException(0xfd3d2814): Index not found\n");
        throw ErrorIdMsg(0xfd3d2814, "Index not found");
    }
#endif
    return result;
}

Mod SteenrodIndices2Mod(const int1d& indices, const MMod1d& basis)
{
    Mod result;
    for (int i : indices)
        result += basis[i];
    return result;
}

void GenerateCells(const int1d& v_degs, const Groebner& gb, int t_max, MMod2d& mbasis, int3d& basis, int3d& inv)
{
    Mod tmp;
    for (size_t i = 0; i < v_degs.size(); ++i) {
        if (v_degs[i] > t_max)
            break;
        ut::get(basis, v_degs[i]).push_back({(int)ut::get(mbasis, v_degs[i]).size()});
        ut::get(mbasis, v_degs[i]).push_back(MMod({}, i));
    }
    Mod1d basis_t_tri;
    Mod1d basis_t;
    for (int t = 0; t <= t_max; ++t) {
        basis_t_tri.clear();
        basis_t.clear();
        for (int n = 1; n <= t; n <<= 1) {
            auto m = MMilnor::Sq(n);
            int t1 = t - n;
            if (t1 >= (int)mbasis.size())
                continue;
            for (const auto& b : mbasis[t1]) {
                auto b_new = gb.Reduce(m * b);
                auto b_new1 = b_new;
                for (auto& b_t : basis_t_tri)
                    if (ut::has(b_new1.data, b_t.GetLead()))
                        b_new1.iaddP(b_t, tmp);
                if (b_new1) {
                    // fmt::print("m={}, b={}, b_new={}\n", m, b, b_new);
                    ut::get(mbasis, t).push_back(b_new1.GetLead());
                    basis_t_tri.push_back(b_new1);
                    basis_t.push_back(b_new);
                }
            }
        }
        if (t < (int)mbasis.size()) {
            std::sort(mbasis[t].begin(), mbasis[t].end());
            for (auto& b : basis_t)
                ut::get(basis, (size_t)t).push_back(SteenrodMod2Indices(b, mbasis[t]));
            int2d image, kernel, g;
            lina::SetLinearMap(basis[t], image, kernel, g);
            ut::get(inv, (size_t)t) = {};
            for (int i = 0; i < (int)image.size(); ++i)
                inv[t].push_back(lina::GetImage(image, g, {i}));
        }
    }
}

/* Generate cell structure from generators and relations */
void GenerateCellStructure(const int1d& v_degs, const Mod1d& rels, int t_max)
{
    Groebner gb(t_max, {}, v_degs);
    gb.AddRels(rels, t_max);
    MMod2d mbasis;
    int3d basis, inv;
    GenerateCells(v_degs, gb, t_max, mbasis, basis, inv);

    json js = "{\"cells_gen\": [], \"cells\": [], \"operations\": []}"_json;
    for (size_t t = 0; t < basis.size(); ++t) {
        for (size_t i = 0; i < basis[t].size(); ++i) {
            js.at("cells").push_back(t);
            // fmt::print("t={}, {}\n", t, SteenrodIndices2Mod(basis[t][i], mbasis[t]));
        }
        // if (basis[t].size())
        // fmt::print("\n");
    }
    for (int d : v_degs) {
        if (d > t_max)
            break;
        js.at("cells_gen").push_back(d);
    }
    for (size_t t = 0; t < basis.size(); ++t) {
        for (size_t i = 0; i < basis[t].size(); ++i) {
            Mod basis_t_i;
            for (int j : basis[t][i])
                basis_t_i.data.push_back(mbasis[t][j]);
            for (size_t j = 1; j + t < basis.size(); j <<= 1) {
                auto y = MMilnor::Sq((uint32_t)j) * basis_t_i;
                y = gb.Reduce(std::move(y));
                int t_y = int(j + t);
                if (y) {
                    if (basis[t_y].size() == 1) {
                        if (basis[t].size() == 1) {
                            js.at("operations").push_back({t, t_y});
                        }
                        else {
                            js.at("operations").push_back({{t, i}, t_y});
                        }
                    }
                    else {
                        int1d iy = SteenrodMod2Indices(y, mbasis[t_y], inv[t_y]);

                        iy.insert(iy.begin(), t_y);
                        if (basis[t].size() == 1) {
                            js.at("operations").push_back({t, iy});
                        }
                        else {
                            js.at("operations").push_back({{t, i}, iy});
                        }
                    }
                }
            }
        }
    }
    fmt::print("{}\n", js.dump(2));
}

struct SmashBasis
{
    size_t t1, t2, i1, i2;
    bool operator<(const SmashBasis& rhs) const
    {
        if (t1 != rhs.t1)
            return t1 < rhs.t1;
        if (i1 != rhs.i1)
            return i1 < rhs.i1;
        return i2 < rhs.i2;
    }
    bool operator==(const SmashBasis& rhs) const
    {
        return t1 == rhs.t1 && i1 == rhs.i1 && i2 == rhs.i2;
    }
};
using SmashBasis1d = std::vector<SmashBasis>;
using SmashBasis2d = std::vector<SmashBasis1d>;

struct SmashGenDeg
{
    size_t t, index;
    bool operator<(const SmashGenDeg& rhs) const
    {
        return t < rhs.t || (t == rhs.t && index < rhs.index);
    }
    bool operator==(const SmashGenDeg& rhs) const
    {
        return t == rhs.t && index == rhs.index;
    }
};
using SmashGenDeg1d = std::vector<SmashGenDeg>;

/**
 * Sort the sequence and each time remove a pair of identical elements
 */
void SortMod2(SmashBasis1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i].t1 = ~0ull;
        }
    ut::RemoveIf(data, [](const SmashBasis& m) { return m.t1 == ~0ull; });
}

/* Generate cell structure from generators and relations */
void GenerateSmash(const int1d& v_degs1, const Mod1d& rels1, const int1d& v_degs2, const Mod1d& rels2, int t_max)
{
    MMod2d mbasis1, mbasis2;
    int3d basis1, basis2, inv1, inv2;
    Groebner gb1(t_max, {}, v_degs1);
    gb1.AddRels(rels1, t_max);
    GenerateCells(v_degs1, gb1, t_max, mbasis1, basis1, inv1);
    Groebner gb2(t_max, {}, v_degs2);
    gb2.AddRels(rels2, t_max);
    GenerateCells(v_degs2, gb2, t_max, mbasis2, basis2, inv2);

    SmashBasis2d basis_smash;
    SmashGenDeg1d gen_degs_smash;
    std::vector<std::pair<size_t, size_t>> indices;
    for (size_t t = 0; t < (size_t)t_max; ++t) {
        for (size_t t1 = 0; t1 < basis1.size(); ++t1) {
            if (t < t1)
                break;
            size_t t2 = t - t1;
            for (size_t i1 = 0; i1 < ut::get(basis1, t1).size(); ++i1) {
                for (size_t i2 = 0; i2 < ut::get(basis2, t2).size(); ++i2) {
                    gen_degs_smash.push_back({t, ut::get(basis_smash, t).size()});
                    basis_smash[t].push_back({t1, t2, i1, i2});
                }
            }
        }
        if (t < basis_smash.size())
            std::sort(basis_smash[t].begin(), basis_smash[t].end());
    }

    Cohomology coh_smash;
    Mod tmp;
    SmashBasis1d op;
    size1d op_indices;
    for (size_t i = 0; i < gen_degs_smash.size(); ++i) {
        const auto& smash_m = basis_smash[gen_degs_smash[i].t][gen_degs_smash[i].index];
        const auto& m1 = mbasis1[smash_m.t1][smash_m.i1];
        const auto& m2 = mbasis2[smash_m.t2][smash_m.i2];

        for (int j = 0; (1 << j) + (int)gen_degs_smash[i].t <= t_max; ++j) {
            Mod rel = MMilnor::P(j, j + 1) * MMod(MMilnor(), i);

            op.clear();
            op_indices.clear();
            for (int k1 = 0; k1 <= (1 << j); ++k1) {
                int k2 = (1 << j) - k1;
                auto sqk1m1 = MMilnor::Sq((uint32_t)k1) * m1;
                sqk1m1 = gb1.Reduce(std::move(sqk1m1));
                auto sqk2m2 = MMilnor::Sq((uint32_t)k2) * m2;
                sqk2m2 = gb2.Reduce(std::move(sqk2m2));

                int t1 = (int)smash_m.t1 + k1;
                int t2 = (int)smash_m.t2 + k2;
                if (t1 < mbasis1.size() && t2 < mbasis2.size()) {
                    int1d i1s = SteenrodMod2Indices(sqk1m1, mbasis1[t1]);
                    int1d i2s = SteenrodMod2Indices(sqk2m2, mbasis2[t2]);
                    for (int i1 : i1s)
                        for (int i2 : i2s)
                            op.push_back({(size_t)t1, (size_t)t2, (size_t)i1, (size_t)i2});
                }
            }
            SortMod2(op);
            for (auto& m : op) {
                int index = ut::IndexOfInSorted(basis_smash[size_t(m.t1 + m.t2)], m);
                ErrorIdMsg::Assert(index != -1, "index != -1");
                int index1 = ut::IndexOf(gen_degs_smash, SmashGenDeg{m.t1 + m.t2, (size_t)index});
                ErrorIdMsg::Assert(index1 != -1, "index1 != -1");
                op_indices.push_back(size_t(index1));
            }

            for (size_t index : op_indices)
                rel.iaddP(MMod(MMilnor(), index), tmp);
            coh_smash.rels.push_back(std::move(rel));
        }
    }
    for (auto& d : gen_degs_smash)
        coh_smash.v_degs.push_back((int)d.t);
    Groebner gb_smash(t_max, {}, coh_smash.v_degs);
    gb_smash.AddRels(coh_smash.rels, t_max, coh_smash.min_rels);
    gb_smash.MinimizeOrderedGensRels(coh_smash.cells, coh_smash.min_rels);

    coh_smash.v_degs = gb_smash.v_degs();
    coh_smash.rels.clear();
    for (int i : coh_smash.min_rels)
        coh_smash.rels.push_back(gb_smash.data()[i]);

    GenerateCellStructure(coh_smash.v_degs, coh_smash.rels, t_max);
}

extern const char* PROGRAM;
extern const char* VERSION;

int main_cellstructure(int argc, char** argv, int& index, const char* desc)
{
    std::string cw;
    int t_max = 100;

    myio::CmdArg1d args = {{"cw", &cw}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    int1d v_degs;
    Mod1d rels;
    GetCoh(v_degs, rels, t_max, cw);
    GenerateCellStructure(v_degs, rels, t_max);
    return 0;
}

int main_smash(int argc, char** argv, int& index, const char* desc)
{
    std::string cw1, cw2;
    int t_max = 100;

    myio::CmdArg1d args = {{"cw1", &cw1}, {"cw2", &cw2}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    int1d v_degs1, v_degs2;
    Mod1d rels1, rels2;
    GetCoh(v_degs1, rels1, t_max, cw1);
    GetCoh(v_degs2, rels2, t_max, cw2);

    GenerateSmash(v_degs1, rels1, v_degs2, rels2, t_max);
    return 0;
}