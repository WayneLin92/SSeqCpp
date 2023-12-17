#include "algebras/linalg.h"
#include "json.h"
#include "main.h"
#include "mylog.h"
#include <filesystem>
#include <fmt/os.h>
#include <fstream>
#include <regex>

struct GenCell
{
    int cell = -1;
    Poly poly;
};
using GenCell1d = std::vector<GenCell>;
using GenCell2d = std::vector<GenCell1d>;

class MyDB : public DBSS
{
    using Statement = myio::Statement;

public:
    MyDB() = default;
    explicit MyDB(const std::string& filename) : DBSS(filename) {}

    void create_pi_basis_products(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_basis_products (id1 INTEGER, id2 INTEGER, prod TEXT, O SMALLINT, PRIMARY KEY(id1, id2));");
    }

    void create_pi_basis_maps_S0() const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS S0_pi_basis_maps (id INTEGER PRIMARY KEY, to_C2 TEXT, to_Ceta TEXT, to_Cnu TEXT, to_Csigma TEXT, O_C2 SMALLINT, O_Ceta SMALLINT, O_Cnu SMALLINT, O_Csigma SMALLINT);");
    }

    GenCell1d load_gen_cells(const std::string& table_prefix) const
    {
        GenCell1d result;
        Statement stmt(*this, "SELECT coalesce(cell, -1), coalesce(cell_coeff, \"\") FROM " + table_prefix + "_generators ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW)
            result.push_back(GenCell{stmt.column_int(0), myio::Deserialize<Poly>(stmt.column_str(1))});
        return result;
    }

    void save_gen_cells(const std::string& table_prefix, const std::vector<std::string>& gen_names, const GenCell1d& gen_cells) const
    {
        Statement stmt(*this, "UPDATE " + table_prefix + "_generators SET name=?1, cell=?2, cell_coeff=?3 WHERE id=?4;");
        for (size_t i = 0; i < gen_names.size(); ++i)
            if (!gen_names[i].empty())
                stmt.bind_and_step(gen_names[i], gen_cells[i].cell, myio::Serialize(gen_cells[i].poly), (int)i);
    }
};

constexpr double BULLET_RADIUS = 0.08f;

double GetRadius(double n)
{
    double length = BULLET_RADIUS * 3 * (n - 1);
    if (length > 1 - BULLET_RADIUS * 6)
        length = 1 - BULLET_RADIUS * 6;
    double sep = n > 1.01 ? length / (n - 1) : BULLET_RADIUS * 3;
    return sep / 3;
}

double radii_get_r(const std::map<AdamsDeg, double>& radii, AdamsDeg deg)
{
    return ut::has(radii, deg) ? radii.at(deg) : BULLET_RADIUS;
}

void SmoothenRadii(std::map<AdamsDeg, double>& radii)
{
    constexpr std::array<AdamsDeg, 8> nbhd = {AdamsDeg(1, 1), AdamsDeg(-1, -1), AdamsDeg(0, 1), AdamsDeg(0, -1), AdamsDeg(1, 2), AdamsDeg(-1, 0), AdamsDeg(1, 0), AdamsDeg(-1, -2)};
    constexpr double BULLET_RADIUS_RATIO = 1.07;
    const double BULLET_RADIUS_RATIO_V2 = std::pow(BULLET_RADIUS_RATIO, 1.414);
    std::map<AdamsDeg, double> radii_ub;
    for (int _ = 0; _ < 4; ++_) {
        radii_ub = radii;
        for (auto& [d, r] : radii) {
            for (size_t i = 0; i < 4; ++i)
                radii_ub[d + nbhd[i]] = std::min(radii_ub[d + nbhd[i]], r * BULLET_RADIUS_RATIO);
            for (size_t i = 4; i < 8; ++i)
                radii_ub[d + nbhd[i]] = std::min(radii_ub[d + nbhd[i]], r * BULLET_RADIUS_RATIO_V2);
        }
        for (auto& [d, r] : radii)
            if (r > radii_ub.at(d) * 1.005)
                r = std::sqrt(r * radii_ub.at(d));
    }
}

void LoadJson(const std::string& diagram_name, nlohmann::json& root_json, nlohmann::json& diag_json);

/* -20 degree */
constexpr double COS_BULLET_ANGLE = 0.9396926228224405;
constexpr double SIN_BULLET_ANGLE = -0.3420201377303428;

void plotBullets(const Staircases1d& nodes_ss, std::variant<const RingSp*, const ModSp*> pCw, const std::map<AdamsDeg, int>& deg2id, nlohmann::json& js)
{
    std::map<AdamsDeg, double> radii;
    for (auto& [d, sc] : nodes_ss.front())
        radii[d] = GetRadius((int)sc.levels.size());
    SmoothenRadii(radii);

    for (auto& [deg, sc] : nodes_ss.front()) {
        int n = (int)sc.levels.size();
        double bottom_right_x = (double)deg.stem() + radii.at(deg) * 1.5 * COS_BULLET_ANGLE * (n - 1);
        double bottom_right_y = (double)deg.s + radii.at(deg) * 1.5 * SIN_BULLET_ANGLE * (n - 1);
        int stable_level = Diagram::GetFirstFixedLevelForPlot(nodes_ss, deg);
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            js["bullets"].push_back(nlohmann::json::object());
            auto& bullet = js["bullets"].back();
            bullet["x"] = bottom_right_x - radii.at(deg) * 3 * COS_BULLET_ANGLE * double((int)i);
            bullet["y"] = bottom_right_y - radii.at(deg) * 3 * SIN_BULLET_ANGLE * double((int)i);
            bullet["r"] = radii.at(deg);
            bullet["b"] = sc.basis[i];

            if (!sc.basis[i].empty()) {
                size_t index = (size_t)sc.basis[i].front();
                bool isGen = std::holds_alternative<const RingSp*>(pCw) ? std::get<const RingSp*>(pCw)->basis.at(deg)[index].IsGen() : std::get<const ModSp*>(pCw)->basis.at(deg)[index].IsGen();
                bullet["c"] = isGen ? "blue" : "black";
            }
            else
                bullet["c"] = "black";

            if (sc.diffs[i] == NULL_DIFF)
                bullet["d"] = nullptr;
            else
                bullet["d"] = sc.diffs[i];

            if (sc.levels[i] >= stable_level) {
                if (sc.diffs[i] == NULL_DIFF)
                    bullet["p"] = R_PERM;
                else
                    bullet["p"] = 10000 - sc.levels[i];
            }
            else if (sc.levels[i] > 5000 || sc.diffs[i] == NULL_DIFF)
                bullet["p"] = R_PERM;
            else
                bullet["p"] = sc.levels[i];

            bullet["l"] = sc.levels[i];
            bullet["i0"] = deg2id.at(deg);
        }
    }
}

void plotBullets(const CofSeq& cofseq, size_t iCs, const Diagram& diagram, const std::map<AdamsDeg, int>& deg2id, nlohmann::json& js)
{
    auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
    auto& nodes_ss = *cofseq.nodes_ss[iCs];

    std::map<AdamsDeg, double> radii;
    for (auto& [d, sc_d] : nodes_cofseq.front())
        radii[d] = (double)sc_d.levels.size();
    for (auto& [d, _] : nodes_ss.front())
        radii[d] += 1.0;
    for (auto& [d, _] : radii)
        radii[d] = GetRadius(radii[d]);

    SmoothenRadii(radii);
    for (auto& [deg, _] : radii) {
        int n = 0;
        if (ut::has(nodes_cofseq.front(), deg))
            n = (int)nodes_cofseq.front().at(deg).levels.size();
        int extra_b = 0;

        if (Diagram::PossMoreEinf(nodes_ss, deg)) {
            extra_b = 1;
            js["bullets_p"].push_back(nlohmann::json::object());
            auto& bullet = js["bullets_p"].back();

            bullet["x"] = (double)deg.stem() + radii.at(deg) * 1.5 * COS_BULLET_ANGLE * (n - 1 + extra_b);
            bullet["y"] = (double)deg.s - radii.at(deg) * 1.5 * SIN_BULLET_ANGLE * (n - 1 + extra_b);
            bullet["r"] = radii.at(deg);
            bullet["c"] = "grey";
        }

        if (!ut::has(nodes_cofseq.front(), deg))
            continue;
        auto& sc = nodes_cofseq.front().at(deg);
        double bottom_right_x = (double)deg.stem() + radii.at(deg) * 1.5 * COS_BULLET_ANGLE * (n - 1 + extra_b);
        double bottom_right_y = (double)deg.s - radii.at(deg) * 1.5 * SIN_BULLET_ANGLE * (n - 1 + extra_b);
        int stable_level = Diagram::GetFirstFixedLevelForPlotCofseq(cofseq, iCs, deg);
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            js["bullets"].push_back(nlohmann::json::object());
            auto& bullet = js["bullets"].back();
            bullet["x"] = bottom_right_x - radii.at(deg) * 3 * COS_BULLET_ANGLE * double((int)i + extra_b);
            bullet["y"] = bottom_right_y + radii.at(deg) * 3 * SIN_BULLET_ANGLE * double((int)i + extra_b);
            bullet["r"] = radii.at(deg);
            bullet["b"] = sc.basis[i];

            if (!sc.basis[i].empty()) {
                size_t index = (size_t)sc.basis[i].front();
                bool isGen = cofseq.isRing[iCs] ? diagram.GetRings()[cofseq.indexCw[iCs]].basis.at(deg)[index].IsGen() : diagram.GetModules()[cofseq.indexCw[iCs]].basis.at(deg)[index].IsGen();
                bullet["c"] = isGen ? "blue" : "black";
            }
            else
                bullet["c"] = "black";

            if (sc.diffs[i] == NULL_DIFF)
                bullet["d"] = nullptr;
            else
                bullet["d"] = sc.diffs[i];

            if (sc.levels[i] >= stable_level) {
                if (sc.diffs[i] == NULL_DIFF)
                    bullet["p"] = R_PERM;
                else
                    bullet["p"] = 10000 - sc.levels[i];
            }
            else if (sc.levels[i] > 5000 || sc.diffs[i] == NULL_DIFF)
                bullet["p"] = R_PERM;
            else
                bullet["p"] = sc.levels[i];

            bullet["l"] = sc.levels[i];
            bullet["i0"] = deg2id.at(deg);
        }
    }
}

void plotRingStrLines(const Staircases1d& nodes_ss, const RingSp& ring, const std::map<AdamsDeg, int>& deg2id, const AdamsDeg1d& degs_factors, const Poly1d& bjs, const int1d& indices_factors, nlohmann::json& js, int forStrl, bool forCofseq = false)
{
    using json = nlohmann::json;
    size_t first_PC = 0;
    for (auto& [deg, sc] : nodes_ss.front()) {
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            Poly bi = Indices2Poly(sc.basis[i], ring.basis.at(deg));
            auto key = std::to_string(deg2id.at(deg) + (int)i);
            for (size_t j = 0; j < degs_factors.size(); ++j) {
                const AdamsDeg deg_prod = degs_factors[j] + deg;
                if (deg_prod.t > ring.t_max)
                    break;
                auto alg_prod = ring.gb.Reduce(bjs[j] * bi);
                if (!alg_prod)
                    continue;
                int1d prod = Poly2Indices(alg_prod, ring.basis.at(deg_prod));
                if (forCofseq) {
                    prod = Residue(std::move(prod), ring.nodes_ss, deg_prod, LEVEL_PERM);
                    if (prod.empty())
                        continue;
                }

                prod = lina::GetInvImage(nodes_ss.front().at(deg_prod).basis, prod);

                auto& prod_json = json::object();
                prod_json["p"] = json::array();
                for (size_t k = 0; k < prod.size(); ++k)
                    prod_json["p"].push_back(deg2id.at(deg_prod) + prod[k]);
                prod_json["l"] = forStrl;
                js["prods"][key].push_back(prod_json);
            }
        }
    }
}

void plotModuleStrLines(const Staircases1d& nodes_ss, const ModSp& mod, const std::map<AdamsDeg, int>& deg2id, const AdamsDeg1d& degs_factors, const Poly1d& bjs, const int1d& indices_factors, nlohmann::json& js, int forStrl, bool forCofseq = false)
{
    using json = nlohmann::json;
    size_t first_PC = 0;
    for (auto& [deg, sc] : nodes_ss.front()) {
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            Mod bi = Indices2Mod(sc.basis[i], mod.basis.at(deg));
            auto key = std::to_string(deg2id.at(deg) + (int)i);
            for (size_t j = 0; j < degs_factors.size(); ++j) {
                const AdamsDeg deg_prod = degs_factors[j] + deg;
                if (deg_prod.t > mod.t_max)
                    break;
                auto alg_prod = mod.gb.Reduce(bjs[j] * bi);
                if (!alg_prod)
                    continue;

                int1d prod = Mod2Indices(alg_prod, mod.basis.at(deg_prod));
                if (forCofseq) {
                    prod = Residue(std::move(prod), mod.nodes_ss, deg_prod, LEVEL_PERM);
                    if (prod.empty())
                        continue;
                }
                prod = lina::GetInvImage(nodes_ss.front().at(deg_prod).basis, prod);

                auto& prod_json = json::object();
                prod_json["p"] = json::array();
                for (size_t k = 0; k < prod.size(); ++k)
                    prod_json["p"].push_back(deg2id.at(deg_prod) + prod[k]);
                prod_json["l"] = forStrl;
                js["prods"][key].push_back(prod_json);
            }
        }
    }
}

std::map<AdamsDeg, int> get_deg2id(const Staircases1d& nodes_ss)
{
    std::map<AdamsDeg, int> deg2id;
    int index = 0;
    for (auto& [d, sc_d] : nodes_ss.front()) {
        deg2id[d] = index;
        index += (int)sc_d.levels.size();
    }
    return deg2id;
}

int main_plot_ss(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;
    std::string plot_dir = root_json.at("dir_website_ss").get<std::string>() + "/" + diag_json.at("dir_plot").get<std::string>();
    myio::AssertFolderExists(plot_dir);

    Diagram diagram(diagram_name, DeduceFlag::no_op);
    const auto& rings = diagram.GetRings();
    const auto& mods = diagram.GetModules();
    const auto& maps = diagram.GetMaps();

    /*
    {
      "type": "ring",
      "gen_names": [],
      "basis": [[1, 2], [1, 2, 3, 4]],
      "bullets": [{"x": 0, "y": 0, "r"(radius): 1, "c"(color): "blue", "b"(basis): [0, 1], "d"(diff): [2], "l"(level): 2, "p"(page): 2, "i0"(index): 0}],
      "degs_factors": [[0, 1], [1, 1]],
      "prods": {"0": [{"p": [0, 1], "l"(is structure line): 0}]},
      "diffs": [{"i": 0, "j": [0, 1], "r": 2}],
      "nds": [{"i": 0, "r": 2}],
    }
    */

    constexpr double BULLET_ANGLE = -20.0f / 180 * 3.1415926f;
    const double COS_BULLET_ANGLE = std::cos(BULLET_ANGLE);
    const double SIN_BULLET_ANGLE = std::sin(BULLET_ANGLE);

    std::vector<std::map<AdamsDeg, int>> all_deg2id;
    size_t cw_size = rings.size() + mods.size();
    for (size_t iCw = 0; iCw < cw_size; ++iCw) {
        bool isRing = iCw < rings.size();
        size_t iMod = isRing ? -1 : iCw - rings.size();
        auto& name = isRing ? rings[iCw].name : mods[iMod].name;
        auto& nodes_ss = isRing ? rings[iCw].nodes_ss : mods[iMod].nodes_ss;
        auto& ring = isRing ? rings[iCw] : rings[mods[iMod].iRing];

        json js;
        js["type"] = isRing ? "ring" : "module";
        if (!isRing)
            js["over"] = ring.name;

        /* gen_names */
        {
            auto path = isRing ? diag_json.at("rings")[iCw].at("path").get<std::string>() : diag_json.at("modules")[iMod].at("path").get<std::string>();
            MyDB db(db_dir + "/" + path);
            auto gen_names = db.load_gen_names(fmt::format("{}_AdamsE2", name));
            auto gen_degs = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
            std::map<AdamsDeg, int> gen_index;
            char letter = isRing ? 'x' : 'v';
            std::string key = isRing ? "gen_names" : "v_names";
            for (size_t i = 0; i < gen_names.size(); ++i) {
                ++gen_index[gen_degs[i]];
                if (gen_names[i].empty())
                    js[key].push_back(fmt::format("{}_{{{},{}{}}}", letter, gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[gen_degs[i]])));
                else
                    js[key].push_back(gen_names[i]);
            }
        }

        /* basis */
        if (isRing) {
            for (auto& [d, basis_d] : ring.basis) {
                for (auto& m : basis_d) {
                    js["basis"].push_back(json::array());
                    for (auto& p : m) {
                        js["basis"].back().push_back(p.g());
                        js["basis"].back().push_back(p.e());
                    }
                }
            }
        }
        else {
            for (auto& [d, basis_d] : mods[iMod].basis) {
                for (auto& m : basis_d) {
                    js["basis"].push_back(json::array());
                    for (auto& p : m.m) {
                        js["basis"].back().push_back(p.g());
                        js["basis"].back().push_back(p.e());
                    }
                    js["basis"].back().push_back(m.v);
                }
            }
        }

        /* ss */
        all_deg2id.push_back(get_deg2id(nodes_ss));
        const auto& deg2id = all_deg2id.back();
        if (isRing)
            plotBullets(nodes_ss, &ring, deg2id, js);
        else
            plotBullets(nodes_ss, &mods[iMod], deg2id, js);

        /* struct lines */
        js["degs_factors"] = json::array();
        js["prods"] = json::object();
        const std::map<AdamsDeg, int>& deg2id_ring = isRing ? deg2id : all_deg2id[mods[iMod].iRing];
        AdamsDeg1d degs_factors;
        Poly1d bjs;
        int1d indices_factors;
        for (auto& strt_factor : diag_json.at("rings")[isRing ? iCw : mods[iMod].iRing].at("plot_factors")[0]) {
            int stem = strt_factor[0].get<int>(), s = strt_factor[1].get<int>(), i_factor = strt_factor[2].get<int>();
            degs_factors.push_back(AdamsDeg(s, stem + s));
            js["degs_factors"].push_back(json::array({stem, s}));
            bjs.push_back(ring.basis.at(degs_factors.back())[i_factor]);
            indices_factors.push_back(deg2id_ring.at(degs_factors.back()) + i_factor);
        }

        if (isRing)
            plotRingStrLines(nodes_ss, ring, deg2id, degs_factors, bjs, indices_factors, js, 1);
        else
            plotModuleStrLines(nodes_ss, mods[iMod], deg2id, degs_factors, bjs, indices_factors, js, 1);

        // Products
        degs_factors.clear();
        bjs.clear();
        indices_factors.clear();
        for (auto& strt_factor : diag_json.at("rings")[isRing ? iCw : mods[iMod].iRing].at("plot_factors")[1]) {
            int stem = strt_factor[0].get<int>(), s = strt_factor[1].get<int>(), i_factor = strt_factor[2].get<int>();
            degs_factors.push_back(AdamsDeg(s, stem + s));
            js["degs_factors"].push_back(json::array({stem, s}));
            bjs.push_back(ring.basis.at(degs_factors.back())[i_factor]);
            indices_factors.push_back(deg2id_ring.at(degs_factors.back()) + i_factor);
        }

        if (isRing)
            plotRingStrLines(nodes_ss, ring, deg2id, degs_factors, bjs, indices_factors, js, 0);
        else
            plotModuleStrLines(nodes_ss, mods[iMod], deg2id, degs_factors, bjs, indices_factors, js, 0);

        /* diff lines */
        js["diffs"] = json::array();
        for (auto& [deg, sc] : nodes_ss.front()) {
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                int src = deg2id.at(deg) + (int)i;
                if (sc.levels[i] > 9000 && sc.diffs[i] != NULL_DIFF) {
                    int r = LEVEL_MAX - sc.levels[i];
                    AdamsDeg deg_tgt = deg + AdamsDeg(r, r - 1);
                    int1d tgt = lina::GetInvImage(nodes_ss.front().at(deg_tgt).basis, sc.diffs[i]);

                    js["diffs"].push_back(json::object());
                    auto& diff_json = js["diffs"].back();
                    diff_json["i"] = src;
                    diff_json["j"] = json::array();
                    for (size_t l = 0; l < tgt.size(); ++l)
                        diff_json["j"].push_back(deg2id.at(deg_tgt) + tgt[l]);
                    diff_json["r"] = r;
                }
            }
        }

        /* unknown diff lines */
        js["nds"] = json::array();
        for (auto& [deg, sc] : nodes_ss.front()) {
            if (deg.stem() <= 126) {
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    int src = deg2id.at(deg) + (int)i;
                    if (sc.levels[i] > 9000 && sc.diffs[i] == NULL_DIFF) {
                        int r = LEVEL_MAX - sc.levels[i];
                        js["nds"].push_back(json::object());
                        auto& nd = js["nds"].back();
                        nd["i"] = src;
                        nd["r"] = r;
                    }
                }
            }
        }

        auto out = fmt::output_file(plot_dir + "/" + name + ".js");
        out.print("globalThis.DATA_JSON_{} = {};\n", name, js.dump());
    }

    /*
    {
      "type": "map",
      "from": "C2",
      "to": "S0",
      "maps": {"2": [0, 1]}
    }
    */
    for (auto& map : maps) {
        if (map->IsMul())
            continue;
        size_t from, to;
        json js;
        js["type"] = "map";
        js["maps"] = json::object();
        auto& map_json = js["maps"];

        if (map->IsFromRing(from)) {
            if (!map->IsToRing(to))
                throw MyException(0x189448f, "Incorrect map type");
            auto& nodes_ss = rings[from].nodes_ss;
            js["from"] = rings[from].name;
            js["to"] = rings[to].name;
            js["sus"] = 0;
            auto& deg2id1 = all_deg2id[from];
            auto& deg2id2 = all_deg2id[to];
            for (auto& [deg, sc] : nodes_ss.front()) {
                if (deg.t > map->t_max)
                    break;
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    int1d fx = map->map(sc.basis[i], deg, diagram);
                    if (!fx.empty()) {
                        auto key = std::to_string(deg2id1.at(deg) + (int)i);
                        map_json[key] = json::array();
                        int1d fx_ss = lina::GetInvImage(rings[to].nodes_ss.front().at(deg).basis, fx);
                        for (int j : fx_ss)
                            map_json[key].push_back(deg2id2.at(deg) + j);
                    }
                }
            }
        }
        else {
            if (map->IsToRing(to)) {
                auto& nodes_ss = mods[from].nodes_ss;
                js["from"] = mods[from].name;
                js["to"] = rings[to].name;
                js["sus"] = -map->deg.stem();
                auto& deg2id1 = all_deg2id[from + rings.size()];
                auto& deg2id2 = all_deg2id[to];
                for (auto& [deg, sc] : nodes_ss.front()) {
                    if (deg.t > map->t_max)
                        break;
                    AdamsDeg deg_fx = deg + map->deg;
                    for (size_t i = 0; i < sc.levels.size(); ++i) {
                        int1d fx = map->map(sc.basis[i], deg, diagram);
                        if (!fx.empty()) {
                            auto key = std::to_string(deg2id1.at(deg) + (int)i);
                            map_json[key] = json::array();
                            int1d fx_ss = lina::GetInvImage(rings[to].nodes_ss.front().at(deg_fx).basis, fx);
                            for (int j : fx_ss)
                                map_json[key].push_back(deg2id2.at(deg_fx) + j);
                        }
                    }
                }
            }
            else {
                auto& nodes_ss = mods[from].nodes_ss;
                js["from"] = mods[from].name;
                js["to"] = mods[to].name;
                js["sus"] = -map->deg.stem();
                auto& deg2id1 = all_deg2id[from + rings.size()];
                auto& deg2id2 = all_deg2id[to + rings.size()];
                for (auto& [deg, sc] : nodes_ss.front()) {
                    if (deg.t > map->t_max)
                        break;
                    AdamsDeg deg_fx = deg + map->deg;
                    for (size_t i = 0; i < sc.levels.size(); ++i) {
                        int1d fx = map->map(sc.basis[i], deg, diagram);
                        if (!fx.empty()) {
                            auto key = std::to_string(deg2id1.at(deg) + (int)i);
                            map_json[key] = json::array();
                            int1d fx_ss = lina::GetInvImage(mods[to].nodes_ss.front().at(deg_fx).basis, fx);
                            for (int j : fx_ss)
                                map_json[key].push_back(deg2id2.at(deg_fx) + j);
                        }
                    }
                }
            }
        }

        auto out = fmt::output_file(plot_dir + "/" + map->name + ".js");
        out.print("globalThis.DATA_JSON_{} = {};\n", map->name, js.dump());
    }

    return 0;
}

int main_plot_cofseq(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;
    std::string plot_dir = root_json.at("dir_website_ss").get<std::string>() + "/" + diag_json.at("dir_plot").get<std::string>();
    myio::AssertFolderExists(plot_dir);

    Diagram diagram(diagram_name, DeduceFlag::cofseq);
    const auto& rings = diagram.GetRings();
    const auto& mods = diagram.GetModules();
    const auto& maps = diagram.GetMaps();
    const auto& cofseqs = diagram.GetCofSeqs();

    /*
    {
      "type": "cofseq",
      "names": ["S0", "C2", "S0"],
      "degs_maps": [0, 1, 0],
      "ss_groups": [
        {
          "bullets": {"x": 0, "y": 0, "r"(radius): 1, "c"(color): "blue", "b"(basis): [0, 1], "d"(diff): [2], "l"(level): 2, "p"(page): 2, "i0"(index): 0},
          "prods": {"i": 0, "j"(factor): 1, "p": [0, 1], "l"(is structure line): 0},
          "diffs": {"i": 0, "j": [0, 1], "r": 2},
          "nds": {"i": 0, "r": 2}
        },
        {}, {}
      ]
    }
    */
    constexpr double BULLET_ANGLE = -20.0f / 180 * 3.1415926f;
    const double COS_BULLET_ANGLE = std::cos(BULLET_ANGLE);
    const double SIN_BULLET_ANGLE = std::sin(BULLET_ANGLE);

    for (auto& cofseq : cofseqs) {
        json js;
        js["type"] = "cofseq";
        js["names"] = cofseq.nameCw;
        js["degs_maps"] = {{cofseq.degMap[0].stem(), cofseq.degMap[0].s}, {cofseq.degMap[1].stem(), cofseq.degMap[1].s}, {cofseq.degMap[2].stem(), cofseq.degMap[2].s}};
        js["cofseq_groups"] = {json::object(), json::object(), json::object()};
        for (size_t iCs = 0; iCs < cofseq.degMap.size(); ++iCs) {
            /* ss */

            auto& jsi = js["cofseq_groups"][iCs];
            jsi["type"] = "cofseq_gp";

            const auto& nodes_ss = *cofseq.nodes_ss[iCs];
            const auto& nodes_cofseq = cofseq.nodes_cofseq[iCs];
            std::map<AdamsDeg, int> deg2id_ss = get_deg2id(nodes_ss);
            std::map<AdamsDeg, int> deg2id_cofseq = get_deg2id(nodes_cofseq);

            plotBullets(cofseq, iCs, diagram, deg2id_ss, jsi);

            /* struct lines */
            jsi["prods"] = json::array();
            AdamsDeg1d degs_factors;
            Poly1d bjs;
            int1d indices_factors;
            std::map<AdamsDeg, int> deg2id_ring;
            size_t iRing = cofseq.isRing[iCs] ? cofseq.indexCw[iCs] : mods[cofseq.indexCw[iCs]].iRing;
            if (cofseq.isRing[iCs])
                deg2id_ring = deg2id_ss;
            else
                deg2id_ring = get_deg2id(rings[iRing].nodes_ss);
            for (auto& strt_factor : diag_json.at("rings")[iRing].at("plot_factors")[0]) {
                int stem = strt_factor[0].get<int>();
                if (stem != 0)
                    continue;
                int s = strt_factor[1].get<int>(), i_factor = strt_factor[2].get<int>();
                degs_factors.push_back(AdamsDeg(s, stem + s));
                bjs.push_back(rings[iRing].basis.at(degs_factors.back())[i_factor]);
                indices_factors.push_back(deg2id_ring.at(degs_factors.back()) + i_factor);
            }

            if (cofseq.isRing[iCs])
                plotRingStrLines(nodes_cofseq, rings[cofseq.indexCw[iCs]], deg2id_cofseq, degs_factors, bjs, indices_factors, jsi, true);
            else
                plotModuleStrLines(nodes_cofseq, mods[cofseq.indexCw[iCs]], deg2id_cofseq, degs_factors, bjs, indices_factors, jsi, true);

            /* diff lines */
            jsi["diffs"] = json::array();
            int stem_map = cofseq.degMap[iCs].stem();
            const auto& nodes_cofseq_next = cofseq.nodes_cofseq[(iCs + 1) % 3];
            std::map<AdamsDeg, int> deg2id_cofseq_next = get_deg2id(nodes_cofseq_next);
            for (auto& [deg, sc] : nodes_cofseq.front()) {
                for (size_t i = 0; i < sc.levels.size(); ++i) {
                    int src = deg2id_cofseq.at(deg) + (int)i;
                    if (sc.levels[i] > 9000 && sc.diffs[i] != NULL_DIFF) {
                        int r = LEVEL_MAX - sc.levels[i];
                        AdamsDeg deg_tgt = deg + AdamsDeg(r, stem_map + r);
                        int1d tgt = lina::GetInvImage(nodes_cofseq_next.front().at(deg_tgt).basis, sc.diffs[i]);

                        jsi["diffs"].push_back(json::object());
                        auto& diff_json = jsi["diffs"].back();
                        diff_json["i"] = src;
                        diff_json["j"] = json::array();
                        for (size_t l = 0; l < tgt.size(); ++l)
                            diff_json["j"].push_back(deg2id_cofseq_next.at(deg_tgt) + tgt[l]);
                        diff_json["r"] = r;
                    }
                }
            }

            /* unknown diff lines */
            jsi["nds"] = json::array();
            for (auto& [deg, sc] : nodes_cofseq.front()) {
                if (deg.stem() <= 126) {
                    for (size_t i = 0; i < sc.levels.size(); ++i) {
                        int src = deg2id_cofseq.at(deg) + (int)i;
                        if (sc.levels[i] > 9000 && sc.diffs[i] == NULL_DIFF) {
                            int r = LEVEL_MAX - sc.levels[i];
                            jsi["nds"].push_back(json::object());
                            auto& nd = jsi["nds"].back();
                            nd["i"] = src;
                            nd["r"] = r;
                        }
                    }
                }
            }
        }

        auto out = fmt::output_file(plot_dir + "/" + cofseq.name + ".js");
        out.print("globalThis.DATA_JSON_{} = {};\n", cofseq.name, js.dump());
    }

    return 0;
}

void ToIndices(const algZ::Poly& x, const PiBasis& basis, const std::map<AdamsDeg, int>& pi_deg2id, int stem, int t_max, int1d& result, int& O)
{
    for (auto& m : x.data) {
        AdamsDeg d = AdamsDeg(m.fil(), stem + m.fil());
        if (m.IsUnKnown() || d.t > t_max) {
            O = m.fil();
            break;
        }
        else {
            const auto& basis_d = basis.at(d).nodes_pi_basis;
            const auto& p = std::lower_bound(basis_d.begin(), basis_d.end(), m);
            int index = (int)(p - basis_d.begin());
            result.push_back(pi_deg2id.at(d) + index);
        }
    }
}

void ToIndices(const algZ::Mod& x, const PiBasisMod& basis, const std::map<AdamsDeg, int>& pi_deg2id, int stem, int t_max, int1d& result, int& O)
{
    for (auto& m : x.data) {
        AdamsDeg d = AdamsDeg(m.fil(), stem + m.fil());
        if (m.IsUnKnown() || d.t > t_max) {
            O = m.fil();
            break;
        }
        else {
            const auto& basis_d = basis.at(d).nodes_pi_basis;
            const auto& p = std::lower_bound(basis_d.begin(), basis_d.end(), m);
            int index = (int)(p - basis_d.begin());
            result.push_back(pi_deg2id.at(d) + index);
        }
    }
}

int main_plot_pi(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Diagram diagram(diagram_name, DeduceFlag::pi);
    /* pi_basis_products */
    int count_ss = 0, count_homotopy = 0;
    diagram.SyncHomotopy(AdamsDeg(0, 0), count_ss, count_homotopy, 0);
    diagram.DeduceTrivialExtensions(0);

    auto& ssS0 = diagram.GetRings();
    auto& ssCofs = diagram.GetModules();
    // auto& all_basis_ss = diagram.GetAllBasisSs();

    // std::vector<MyDB> dbPlots;
    // dbPlots.reserve(dbnames.size()); /* To avoid reallocation is critical */
    // std::vector<std::string> tablesE2, complexNames;
    // std::vector<std::map<AdamsDeg, int>> pi_deg2ids(dbnames.size());
    // for (size_t i = 0; i < dbnames.size(); ++i) {
    //     std::string dbname = dbnames[i];
    //     dbname.insert(dbname.size() - 3, "_plot");
    //     dbPlots.emplace_back(dbname);
    //     tablesE2.push_back(GetE2TablePrefix(dbname));
    //     complexNames.push_back(GetComplexName(dbname));
    //     MyDB db(dbnames[i]);
    // }

    // for (auto& db : dbPlots)
    //     db.begin_transaction();

    //{
    //    algZ::Mon1d S0_pi_basis;
    //    AdamsDeg1d S0_pi_basis_deg;

    //    int pi_index = 0;
    //    AdamsDeg1d S0_degs = OrderDegsV2(diagram.GetRings().nodes_pi_basis.front());
    //    for (AdamsDeg d : S0_degs) {
    //        auto& basis_d = diagram.GetRings().nodes_pi_basis.front().at(d).nodes_pi_basis;
    //        pi_deg2ids[0][d] = pi_index;
    //        pi_index += (int)basis_d.size();
    //        for (auto& b : basis_d) {
    //            S0_pi_basis.push_back(b);
    //            S0_pi_basis_deg.push_back(d);
    //        }
    //    }

    //    algZ::MMod2d Cofs_pi_basis(diagram.GetModules().size());
    //    AdamsDeg2d Cofs_pi_basis_deg(diagram.GetModules().size());

    //    for (size_t iCof = 0; iCof < diagram.GetModules().size(); ++iCof) {
    //        auto& Cof = diagram.GetModules()[iCof];
    //        int pi_index = 0;
    //        AdamsDeg1d Cof_degs = OrderDegsV2(Cof.nodes_pi_basis.front());
    //        for (AdamsDeg d : Cof_degs) {
    //            auto& basis_d = Cof.nodes_pi_basis.front().at(d).nodes_pi_basis;
    //            pi_deg2ids[iCof + 1][d] = pi_index;
    //            pi_index += (int)basis_d.size();
    //            for (auto& b : basis_d) {
    //                Cofs_pi_basis[iCof].push_back(b);
    //                Cofs_pi_basis_deg[iCof].push_back(d);
    //            }
    //        }
    //    }

    //    const int1d arr_factors = {1, 3, 7, 15, 23, 29, 33, 39, 40, 42, 47};
    //    {
    //        auto& pi_gb = diagram.GetRings().pi_gb;
    //        int t_max = diagram.GetRings().t_max;
    //        dbPlots[0].drop_and_create_pi_basis_products(complexNames[0]);
    //        myio::Statement stmt(dbPlots[0], "INSERT INTO " + complexNames[0] + "_pi_basis_products (id1, id2, prod, O) VALUES (?1, ?2, ?3, ?4);");
    //        for (int i : arr_factors) {
    //            for (size_t j = 0; j < S0_pi_basis.size(); ++j) {
    //                const AdamsDeg deg_prod = S0_pi_basis_deg[i] + S0_pi_basis_deg[j];
    //                if (deg_prod.t <= t_max) {
    //                    algZ::Poly poly_prod = pi_gb.ReduceV2(S0_pi_basis[i] * S0_pi_basis[j]);
    //                    int1d prod;
    //                    int O = -1;
    //                    ToIndices(poly_prod, diagram.GetRings().nodes_pi_basis.front(), pi_deg2ids[0], deg_prod.stem(), t_max, prod, O);

    //                    if (O != -1 || !prod.empty()) {
    //                        stmt.bind_and_step((int)i, (int)j, myio::Serialize(prod), O);
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    dbPlots[0].drop_and_create_pi_basis_maps_S0();
    //    for (size_t iSS = 1; iSS < dbPlots.size(); ++iSS) {
    //        size_t iCof = iSS - 1;
    //        auto& pi_gb = ssCofs[iCof].pi_gb;
    //        int t_max = ssCofs[iCof].t_max;

    //        /* module structure */
    //        {
    //            dbPlots[iSS].drop_and_create_pi_basis_products(complexNames[iSS]);
    //            myio::Statement stmt(dbPlots[iSS], "INSERT INTO " + complexNames[iSS] + "_pi_basis_products (id1, id2, prod, O) VALUES (?1, ?2, ?3, ?4);");
    //            for (int i : arr_factors) {
    //                for (size_t j = 0; j < Cofs_pi_basis[iCof].size(); ++j) {
    //                    const AdamsDeg deg_prod = S0_pi_basis_deg[i] + Cofs_pi_basis_deg[iCof][j];
    //                    if (deg_prod.t <= t_max) {
    //                        algZ::Mod x_prod = pi_gb.ReduceV2(S0_pi_basis[i] * Cofs_pi_basis[iCof][j]);
    //                        int1d prod;
    //                        int O = -1;
    //                        ToIndices(x_prod, diagram.GetModules()[iCof].nodes_pi_basis.front(), pi_deg2ids[iSS], deg_prod.stem(), t_max, prod, O);

    //                        if (O != -1 || !prod.empty()) {
    //                            stmt.bind_and_step((int)i, (int)j, myio::Serialize(prod), O);
    //                        }
    //                    }
    //                }
    //            }
    //        }

    //        /* S0->Cof */
    //        {
    //            myio::Statement stmt(dbPlots[0], "INSERT INTO S0_pi_basis_maps (id, to_" + complexNames[iSS] + ", O_" + complexNames[iSS] + ") VALUES (?1, ?2, ?3) ON CONFLICT DO UPDATE SET to_" + complexNames[iSS] + "=excluded.to_"
    //                                                 + complexNames[iSS] + ", O_" + complexNames[iSS] + "=excluded.O_" + complexNames[iSS]);
    //            for (size_t i = 0; i < S0_pi_basis.size(); ++i) {
    //                const AdamsDeg deg = S0_pi_basis_deg[i];
    //                if (deg.t <= t_max) {
    //                    algZ::Mod x = pi_gb.ReduceV2(algZ::Mod(S0_pi_basis[i], 0, 0));
    //                    int1d image;
    //                    int O = -1;
    //                    ToIndices(x, diagram.GetModules()[iCof].nodes_pi_basis.front(), pi_deg2ids[iSS], deg.stem(), t_max, image, O);

    //                    if (O != -1 || !image.empty()) {
    //                        stmt.bind_and_step((int)i, myio::Serialize(image), O);
    //                    }
    //                }
    //            }
    //        }

    //        /* Cof->S0 */
    //        {
    //            dbPlots[iSS].drop_and_create_pi_basis_maps_Cof(complexNames[iSS]);
    //            myio::Statement stmt(dbPlots[iSS], "INSERT INTO " + complexNames[iSS] + "_pi_basis_maps (id, to_S0, O_S0) VALUES (?1, ?2, ?3)");
    //            for (size_t i = 0; i < Cofs_pi_basis[iCof].size(); ++i) {
    //                const AdamsDeg deg = Cofs_pi_basis_deg[iCof][i] - ssCofs[iCof].deg_qt;
    //                if (deg.t <= ssS0.t_max) {
    //                    algZ::Poly a = ssS0.pi_gb.ReduceV2(algZ::subs(Cofs_pi_basis[iCof][i], ssCofs[iCof].nodes_pi_qt.back(), ssCofs[iCof].pi_gb.v_degs()));
    //                    diagram.ExtendRelRing(deg.stem(), a);
    //                    int1d image;
    //                    int O = -1;
    //                    ToIndices(a, diagram.GetRings().nodes_pi_basis.front(), pi_deg2ids[0], deg.stem(), t_max, image, O);

    //                    if (O != -1 || !image.empty()) {
    //                        stmt.bind_and_step((int)i, myio::Serialize(image), O);
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}

    // for (auto& db : dbPlots)
    //     db.end_transaction();

    // diagram.save(dbnames, DeduceFlag::homotopy);
    return 0;
}

void AddCellColumn(MyDB& db, const std::string& name)
{
    std::string table = fmt::format("{}_AdamsE2_generators", name);
    if (!db.has_column(table, "cell")) {
        db.add_column(table, "cell SMALLINT");
        db.add_column(table, "cell_coeff TEXT");
    }
}

int main_rename_gen_reset(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;

    auto& rings = diag_json.at("rings");
    for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
        auto name = rings[iRing].at("name").get<std::string>();
        auto path = diag_json.at("rings")[iRing].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);

        db.begin_transaction();
        db.execute_cmd(fmt::format("UPDATE {}_AdamsE2_generators SET name=NULL", name));
        db.end_transaction();
    }
    auto& mods = diag_json.at("modules");
    for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
        auto name = mods[iMod].at("name").get<std::string>();
        auto path = mods[iMod].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);

        db.begin_transaction();
        AddCellColumn(db, name);
        db.execute_cmd(fmt::format("UPDATE {}_AdamsE2_generators SET name=NULL, cell=NULL, cell_coeff=NULL", name));
        db.end_transaction();
    }
    return 0;
}

int main_rename_gen_export(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";
    std::string gen_names_json;

    myio::CmdArg1d args = {{"gen_names_json", &gen_names_json}};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;

    Diagram diagram(diagram_name, DeduceFlag::no_op);
    const auto& rings = diagram.GetRings();
    const auto& mods = diagram.GetModules();

    /*
    {
        "S0": [
            [stem, s, index, "gen_name"]
        ]
    {
    */

    json js;
    size_t cw_size = rings.size() + mods.size();
    for (size_t iCw = 0; iCw < cw_size; ++iCw) {
        bool isRing = iCw < rings.size();
        size_t iMod = isRing ? -1 : iCw - rings.size();
        auto& name = isRing ? rings[iCw].name : mods[iMod].name;

        js[name] = json::array();

        auto path = isRing ? diag_json.at("rings")[iCw].at("path").get<std::string>() : diag_json.at("modules")[iMod].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        auto gen_names = db.load_gen_names(fmt::format("{}_AdamsE2", name));
        auto gen_degs = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));

        std::map<AdamsDeg, int> gen_index;
        for (size_t i = 0; i < gen_names.size(); ++i) {
            ++gen_index[gen_degs[i]];
            if (!gen_names[i].empty())
                js[name].push_back({gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]], gen_names[i]});
        }

        if (js[name].empty())
            js.erase(name);
    }

    auto out = fmt::output_file(gen_names_json);
    out.print("{}", js.dump());

    return 0;
}

int main_rename_gen_import(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";
    std::string gen_names_json;

    myio::CmdArg1d args = {{"gen_names_json", &gen_names_json}};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;

    json js;
    {
        std::ifstream ifs(gen_names_json);
        if (ifs.is_open())
            ifs >> js;
        else {
            fmt::print("File {} not found\n", gen_names_json);
            throw MyException(0x13889aa5, "File not found");
        }
    }

    auto& rings_json = diag_json.at("rings");
    for (auto& [name, json_gen_names] : js.items()) {
        size_t iRing = 0;
        for (; iRing < rings_json.size(); ++iRing) {
            if (rings_json[iRing].at("name").get<std::string>() == name)
                break;
        }
        if (iRing < rings_json.size()) {
            auto path = rings_json[iRing].at("path").get<std::string>();
            MyDB db(db_dir + "/" + path);
            auto gen_degs = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
            myio::string1d gen_names(gen_degs.size());

            std::map<AdamsDeg, int1d> gen_indices;
            for (size_t i = 0; i < gen_names.size(); ++i)
                gen_indices[gen_degs[i]].push_back((int)i);

            for (auto& json_gen_name : json_gen_names) {
                int stem = json_gen_name[0].get<int>();
                int s = json_gen_name[1].get<int>();
                int index = json_gen_name[2].get<int>();
                std::string gen_name = json_gen_name[3].get<std::string>();
                AdamsDeg deg(s, stem + s);
                gen_names[gen_indices.at(deg)[index - 1]] = std::move(gen_name);
            }

            db.begin_transaction();
            db.save_gen_names(fmt::format("{}_AdamsE2", name), gen_names);
            db.end_transaction();
        }
        else {
            size_t iMod = 0;
            auto& mods_json = diag_json.at("modules");
            for (; iMod < mods_json.size(); ++iMod) {
                if (mods_json[iMod].at("name").get<std::string>() == name)
                    break;
            }
            if (iMod >= mods_json.size())
                continue;
            auto path = mods_json[iMod].at("path").get<std::string>();
            MyDB db(db_dir + "/" + path);
            auto gen_degs = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
            myio::string1d gen_names(gen_degs.size());
            GenCell1d gen_cells(gen_degs.size());

            std::map<AdamsDeg, int1d> gen_indices;
            for (size_t i = 0; i < gen_names.size(); ++i)
                gen_indices[gen_degs[i]].push_back((int)i);

            for (auto& json_gen_name : json_gen_names) {
                int stem = json_gen_name[0].get<int>();
                int s = json_gen_name[1].get<int>();
                int index = json_gen_name[2].get<int>();
                std::string gen_name = json_gen_name[3].get<std::string>();
                int cell = json_gen_name[4].get<int>();
                Poly cell_coeff = myio::Deserialize<Poly>(json_gen_name[5].get<std::string>());
                AdamsDeg deg(s, stem + s);
                gen_names[gen_indices.at(deg)[index - 1]] = std::move(gen_name);
                gen_cells[gen_indices.at(deg)[index - 1]] = GenCell{cell, cell_coeff};
            }

            db.begin_transaction();
            AddCellColumn(db, name);
            db.save_gen_cells(fmt::format("{}_AdamsE2", name), gen_names, gen_cells);
            db.end_transaction();
        }
    }

    return 0;
}

int main_rename_gen_defaultS0(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;

    auto& rings_json = diag_json.at("rings");

    size_t iRing = 0;
    for (; iRing < rings_json.size(); ++iRing) {
        if (rings_json[iRing].at("name").get<std::string>() == "S0")
            break;
    }
    if (iRing < rings_json.size()) {
        auto path = rings_json[iRing].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        auto gen_degs = db.load_gen_adamsdegs("S0_AdamsE2");
        myio::string1d gen_names = db.load_gen_names("S0_AdamsE2");

        std::map<AdamsDeg, int> gen_index;
        for (size_t i = 0; i < gen_names.size(); ++i) {
            ++gen_index[gen_degs[i]];
            if (gen_names[i].empty())
                gen_names[i] = fmt::format("x_{{{},{}{}}}", gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[gen_degs[i]]));
        }

        db.begin_transaction();
        db.save_gen_names("S0_AdamsE2", gen_names);
        db.end_transaction();
    }

    return 0;
}

int main_rename_gen_cell(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;

    Diagram diagram(diagram_name, DeduceFlag::no_op);
    const auto& rings = diagram.GetRings();
    const auto& mods = diagram.GetModules();

    for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
        auto& name = mods[iMod].name;

        auto path = diag_json.at("modules")[iMod].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        auto gen_degs = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
        myio::string1d gen_names(gen_degs.size());
        GenCell1d gen_cells(gen_degs.size());

        std::map<AdamsDeg, int> gen_index;
        for (size_t i = 0; i < gen_names.size(); ++i) {
            ++gen_index[gen_degs[i]];
            if (gen_degs[i].s == 0) {
                gen_names[i] = fmt::format("[{}]", gen_degs[i].stem());
                gen_cells[i] = GenCell{gen_degs[i].stem(), Mon()};
                fmt::print("{}: v_{{{},{}{}}} --> {}\n", name, gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[gen_degs[i]]), gen_names[i]);
            }
        }

        db.begin_transaction();
        AddCellColumn(db, name);
        db.save_gen_cells(fmt::format("{}_AdamsE2", name), gen_names, gen_cells);
        db.end_transaction();
    }

    return 0;
}

bool IsNamed(const Poly& p, const myio::string1d& gen_names)
{
    for (auto& m : p.data)
        for (auto& ge : m)
            if (gen_names[ge.g()].empty())
                return false;
    return true;
}

bool IsNamed(const Mod& p, const GenCell1d& gen_cells, const myio::string1d& gen_names_ring)
{
    for (auto& m : p.data) {
        for (auto& ge : m.m)
            if (gen_names_ring[ge.g()].empty())
                return false;
        if (gen_cells[m.v].cell == -1)
            return false;
    }
    return true;
}

std::string strMon(const Mon& m, const myio::string1d& gen_names)
{
    myio::string1d strs;
    const std::regex no_need_bracket("^(?:(?:\\\\[a-zA-Z]+(?:\\{\\}|)|[a-zA-Z])(?:_\\d|_\\{\\d+\\}|)|(?:\\(|\\[).*(?:\\)|\\]))$"); /* match example: a_1, \alpha_{31} */
    std::smatch match;
    for (auto& ge : m) {
        std::string gen_name = gen_names[ge.g()];
        if (std::regex_search(gen_name, match, no_need_bracket); !match[0].matched)
            gen_name = fmt::format("({})", gen_name);
        int e = (int)ge.e();
        if (e == 1)
            strs.push_back(std::move(gen_name));
        else if (e < 10)
            strs.push_back(fmt::format("{}^{}", gen_name, e));
        else
            strs.push_back(fmt::format("{}^{{{}}}", gen_name, e));
    }
    std::string result = myio::join(" ", strs);
    if (result.empty())
        result = "1";
    return result;
}

std::string StrPoly(const Poly& p, const myio::string1d& gen_names)
{
    std::string result = myio::Tpljoin("+", p.data.begin(), p.data.end(), [&gen_names](const Mon& m) { return strMon(m, gen_names); });
    if (result.empty())
        result = "0";
    return result;
}

std::string StrGenCell(const GenCell& gen_cell, const myio::string1d& gen_names)
{
    std::string strPoly = StrPoly(gen_cell.poly, gen_names);
    const std::regex is_sum("^.*\\+.*$");                                                                     /* match example: a+b */
    const std::regex has_outer_brackets("^\\((?:[^)(]|\\((?:[^)(]|\\((?:[^)(]|\\([^)(]*\\))*\\))*\\))*\\)$"); /* match example: (a(((d))bc)d) */
    std::smatch match;
    if (std::regex_search(strPoly, match, is_sum); match[0].matched) {
        if (std::regex_search(strPoly, match, has_outer_brackets); !match[0].matched)
            return fmt::format("(({})[{}])", strPoly, gen_cell.cell);
        else
            return fmt::format("({}[{}])", strPoly, gen_cell.cell);
    }
    else
        return fmt::format("({}[{}])", strPoly, gen_cell.cell);
}

GenCell TopCell(const Mod& p, const GenCell1d& gen_cells, const RingSp& ring)
{
    GenCell result;
    GenCell1d cells;
    for (auto& m : p.data) {
        if (gen_cells[m.v].poly.data.size() > 10000) {
            fmt::print("{} {}\n", m.v, gen_cells[m.v].cell);
        }
        auto poly = ring.gb.Reduce(m.m * gen_cells[m.v].poly);
        if (poly)
            cells.push_back(GenCell{gen_cells[m.v].cell, std::move(poly)});
    }
    if (!cells.empty()) {
        int top_cell = std::max_element(cells.begin(), cells.end(), [](const GenCell& c1, const GenCell& c2) { return c1.cell < c2.cell; })->cell;
        result.cell = top_cell;
        ut::RemoveIf(cells, [top_cell](const GenCell& c) { return c.cell < top_cell; });
        for (auto& c : cells)
            result.poly += c.poly;
        if (!result.poly)
            result.cell = -1;
    }
    return result;
}

void ResolveNameConflict(myio::string1d& gen_names)
{
    std::set<std::string> set_gen_names;
    for (auto& name : gen_names) {
        if (!name.empty()) {
            if (ut::has(set_gen_names, name)) {
                int index = 2;
                while (ut::has(set_gen_names, fmt::format("{}_{}", name, index)))
                    ++index;
                name += fmt::format("_{}", index);
            }
            set_gen_names.insert(name);
        }
    }
}

int main_rename_gen_pull_back(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;

    Diagram diagram(diagram_name, DeduceFlag::no_op);
    const auto& rings = diagram.GetRings();
    const auto& mods = diagram.GetModules();
    const auto& maps = diagram.GetMaps();
    std::vector<std::vector<std::string>> gen_names_rings(rings.size());
    std::vector<std::vector<std::string>> gen_names_mods(mods.size());
    AdamsDeg2d gen_degs_all(rings.size() + mods.size());
    GenCell2d gen_cells_mods;
    for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
        auto& name = rings[iRing].name;
        auto path = diag_json.at("rings")[iRing].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        gen_names_rings[iRing] = db.load_gen_names(fmt::format("{}_AdamsE2", name));
        gen_degs_all[iRing] = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
    }
    for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
        auto& name = mods[iMod].name;
        auto path = diag_json.at("modules")[iMod].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        std::string table = fmt::format("{}_AdamsE2_generators", name);
        if (!db.has_column(table, "cell")) {
            db.add_column(table, "cell SMALLINT");
            db.add_column(table, "cell_coeff TEXT");
        }
        gen_names_mods[iMod] = db.load_gen_names(fmt::format("{}_AdamsE2", name));
        gen_cells_mods.push_back(db.load_gen_cells(fmt::format("{}_AdamsE2", name)));
        gen_degs_all[rings.size() + iMod] = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
    }

    int count = 0, old_count = -1;
    while (count != old_count) {
        old_count = count;
        for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
            auto& ring = rings[iRing];
            auto& name = rings[iRing].name;
            auto& gen_names = gen_names_rings[iRing];
            auto& gen_degs = gen_degs_all[iRing];

            std::map<AdamsDeg, int> gen_index;
            for (size_t i = 0; i < gen_names.size(); ++i) {
                AdamsDeg deg = gen_degs[i];
                ++gen_index[deg];
                if (!gen_names[i].empty())
                    continue;

                for (size_t iMap : ring.ind_maps) {
                    auto map = (MapRing2Ring*)maps[iMap].get();
                    if (deg.t <= map->t_max) {
                        auto f_gen = rings[map->to].gb.Reduce(map->images[i]);
                        auto& gen_names_tgt = gen_names_rings[map->to];
                        if (f_gen && IsNamed(f_gen, gen_names_tgt)) {
                            gen_names[i] = StrPoly(f_gen, gen_names_tgt);
                            auto old_gen_name = fmt::format("x_{{{},{}{}}}", gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[gen_degs[i]]));
                            fmt::print("({}) {} {} --> {}\n", map->display, name, old_gen_name, gen_names[i]);
                            ++count;
                            break;
                        }
                    }
                }
            }
        }
        for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
            auto& mod = mods[iMod];
            auto& name = mods[iMod].name;
            auto& gen_cells = gen_cells_mods[iMod];
            auto& gen_names = gen_names_mods[iMod];
            auto& gen_degs = gen_degs_all[rings.size() + iMod];

            std::map<AdamsDeg, int> gen_index;
            for (size_t i = 0; i < gen_cells.size(); ++i) {
                AdamsDeg deg = gen_degs[i];
                ++gen_index[deg];
                if (!gen_names[i].empty())
                    continue;

                for (size_t iMap : mod.ind_maps) {
                    auto map1 = maps[iMap].get();
                    size_t iMap_json = 0;
                    while (diag_json.at("maps")[iMap_json].at("name") != map1->name)
                        ++iMap_json;
                    if (!diag_json.at("maps")[iMap_json].contains("type") || diag_json.at("maps")[iMap_json].at("type") != "quo_skeleton")
                        continue;
                    if (deg.t > map1->t_max)
                        continue;
                    if (auto map = dynamic_cast<MapMod2Ring*>(map1); map != nullptr) {
                        if (map->deg.s == 0) {
                            auto& f_gen = map->images[i];
                            auto& gen_names_tgt = gen_names_rings[map->to];
                            if (f_gen && IsNamed(f_gen, gen_names_tgt)) {
                                gen_cells[i] = GenCell{map->deg.s, f_gen};
                                gen_names[i] = StrGenCell(gen_cells[i], gen_names_tgt);
                                auto old_gen_name = fmt::format("v_{{{},{}{}}}", gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[gen_degs[i]]));
                                fmt::print("({}) {} {} --> {}\n", map->display, name, old_gen_name, gen_names[i]);
                                ++count;
                                break;
                            }
                        }
                    }
                    else if (auto map = dynamic_cast<MapMod2Mod*>(map1); map != nullptr) {
                        if (map->deg.s == 0) {
                            auto& f_gen = map->images[i];
                            auto& gen_cells_tgt = gen_cells_mods[map->to];
                            auto& gen_names_ring_tgt = gen_names_rings[mods[map->to].iRing];
                            if (f_gen && IsNamed(f_gen, gen_cells_tgt, gen_names_ring_tgt)) {
                                auto gen_cell = TopCell(f_gen, gen_cells_tgt, rings[mods[map->to].iRing]);
                                if (gen_cell.cell != -1 && IsNamed(gen_cell.poly, gen_names_ring_tgt)) {
                                    gen_cells[i] = std::move(gen_cell);
                                    gen_cells[i].cell += -map->deg.stem();
                                    gen_names[i] = StrGenCell(gen_cells[i], gen_names_ring_tgt);
                                    auto old_gen_name = fmt::format("v_{{{},{}{}}}", gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[gen_degs[i]]));
                                    fmt::print("({}) {} {} --> {}\n", map->display, name, old_gen_name, gen_names[i]);
                                    ++count;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
        auto& name = rings[iRing].name;
        auto path = diag_json.at("rings")[iRing].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        ResolveNameConflict(gen_names_rings[iRing]);
        db.begin_transaction();
        db.execute_cmd(fmt::format("UPDATE {}_AdamsE2_generators SET name=NULL", name));
        db.save_gen_names(fmt::format("{}_AdamsE2", name), gen_names_rings[iRing]);
        db.end_transaction();
    }
    for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
        auto& name = mods[iMod].name;
        auto& gen_cells = gen_cells_mods[iMod];
        auto path = diag_json.at("modules")[iMod].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        ResolveNameConflict(gen_names_mods[iMod]);
        db.begin_transaction();
        db.execute_cmd(fmt::format("UPDATE {}_AdamsE2_generators SET name=NULL, cell=NULL, cell_coeff=NULL", name));
        db.save_gen_cells(fmt::format("{}_AdamsE2", name), gen_names_mods[iMod], gen_cells_mods[iMod]);
        db.end_transaction();
    }
    return 0;
}

int main_rename_gen_push_forward(int argc, char** argv, int& index, const char* desc)
{
    std::string diagram_name = "default";
    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"diagram", &diagram_name}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    using json = nlohmann::json;
    json root_json, diag_json;
    LoadJson(diagram_name, root_json, diag_json);
    std::string db_dir = root_json["diagrams"].contains(diagram_name) ? root_json["diagrams"][diagram_name].get<std::string>() : diagram_name;

    Diagram diagram(diagram_name, DeduceFlag::no_op);
    const auto& rings = diagram.GetRings();
    const auto& mods = diagram.GetModules();
    const auto& maps = diagram.GetMaps();
    std::vector<std::vector<std::string>> gen_names_rings(rings.size());
    std::vector<std::vector<std::string>> gen_names_mods(mods.size());
    AdamsDeg2d gen_degs_all(rings.size() + mods.size());
    GenCell2d gen_cells_mods;
    for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
        auto& name = rings[iRing].name;
        auto path = diag_json.at("rings")[iRing].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        gen_names_rings[iRing] = db.load_gen_names(fmt::format("{}_AdamsE2", name));
        gen_degs_all[iRing] = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
    }
    for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
        auto& name = mods[iMod].name;
        auto path = diag_json.at("modules")[iMod].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        std::string table = fmt::format("{}_AdamsE2_generators", name);
        if (!db.has_column(table, "cell")) {
            db.add_column(table, "cell SMALLINT");
            db.add_column(table, "cell_coeff TEXT");
        }
        gen_names_mods[iMod] = db.load_gen_names(fmt::format("{}_AdamsE2", name));
        gen_cells_mods.push_back(db.load_gen_cells(fmt::format("{}_AdamsE2", name)));
        gen_degs_all[rings.size() + iMod] = db.load_gen_adamsdegs(fmt::format("{}_AdamsE2", name));
    }
    std::vector<std::map<AdamsDeg, int1d>> gen_indices_all(gen_degs_all.size());
    for (size_t iCw = 0; iCw < gen_degs_all.size(); ++iCw)
        for (size_t i = 0; i < gen_degs_all[iCw].size(); ++i)
            gen_indices_all[iCw][gen_degs_all[iCw][i]].push_back((int)i);

    int count = 0, old_count = -1;
    while (count != old_count) {
        old_count = count;
        /*for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
            auto& ring = rings[iRing];
            auto& name = rings[iRing].name;
            auto& gen_names = gen_names_rings[iRing];
            auto& gen_degs = gen_degs_all[iRing];

            std::map<AdamsDeg, int> gen_index;
            for (size_t i = 0; i < gen_names.size(); ++i) {
                AdamsDeg deg = gen_degs[i];
                ++gen_index[deg];
                if (!gen_names[i].empty())
                    continue;

                for (size_t iMap : ring.ind_maps) {
                    auto& map = maps[iMap];
                    if (deg.t <= map->t_max) {
                        auto& f = std::get<MapRing2Ring>(map.map);
                        auto& f_gen = rings[to].gb.Reduce(f.images[i]);
                        auto& gen_names_tgt = gen_names_rings[to];
                        if (f_gen && IsNamed(f_gen, gen_names_tgt)) {
                            gen_names[i] = StrPoly(f_gen, gen_names_tgt);
                            auto& old_gen_name = fmt::format("x_{{{},{}{}}}", gen_degs[i].stem(), gen_degs[i].s, gen_index[gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[gen_degs[i]]));
                            fmt::print("{}: ({}) --> {}\n", name, old_gen_name, gen_names[i]);
                            ++count;
                            break;
                        }
                    }
                }
            }
        }*/
        for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
            auto& mod = mods[iMod];
            auto& name = mods[iMod].name;
            auto& gen_cells = gen_cells_mods[iMod];
            auto& gen_names = gen_names_mods[iMod];
            auto& gen_degs = gen_degs_all[rings.size() + iMod];
            auto& gen_indices = gen_indices_all[rings.size() + iMod];

            for (auto& [deg, indices] : gen_indices) {
                if (std::all_of(indices.begin(), indices.end(), [&gen_names](int i) { return gen_names[i].empty(); }))
                    continue;
                for (size_t iMap : mod.ind_maps) {
                    auto map = dynamic_cast<MapMod2Mod*>(maps[iMap].get());
                    if (map == nullptr)
                        continue;
                    size_t iMap_json = 0;
                    while (diag_json.at("maps")[iMap_json].at("name") != map->name)
                        ++iMap_json;
                    if (!diag_json.at("maps")[iMap_json].contains("type") || diag_json.at("maps")[iMap_json].at("type") != "skeleton")
                        continue;
                    if (deg.t > map->t_max)
                        continue;
                    if (map->deg.s != 0)
                        continue;
                    auto& gen_names_tgt = gen_names_mods[map->to];
                    auto& gen_cells_tgt = gen_cells_mods[map->to];
                    AdamsDeg deg_tgt = deg + map->deg;
                    if (!ut::has(gen_indices_all[rings.size() + map->to], deg_tgt))
                        continue;
                    auto& indices_tgt = gen_indices_all[rings.size() + map->to].at(deg_tgt);
                    if (std::all_of(indices_tgt.begin(), indices_tgt.end(), [&gen_names_tgt](int i) { return !gen_names_tgt[i].empty(); }))
                        continue;
                    int2d fxs, image, kernel, g;
                    for (size_t i = 0; i < mod.basis.at(deg).size(); ++i)
                        fxs.push_back(map->map({int(i)}, deg, diagram));
                    lina::SetLinearMap(fxs, image, kernel, g);

                    auto& mod_tgt = mods[map->to];
                    int gen_index = 0;
                    for (int i : indices_tgt) {
                        ++gen_index;
                        if (!gen_names_tgt[i].empty())
                            continue;
                        MMod mi({}, i);
                        int1d fx = Mod2Indices(mi, mod_tgt.basis.at(deg_tgt));
                        if (lina::Residue(image, fx).empty()) {
                            int1d x = lina::GetImage(image, g, fx);
                            Mod x_alg = Indices2Mod(x, mod.basis.at(deg));
                            if (IsNamed(x_alg, gen_cells, gen_names)) {
                                auto gen_cell = TopCell(x_alg, gen_cells, rings[mod.iRing]);
                                if (gen_cell.cell != -1 && IsNamed(gen_cell.poly, gen_names_rings[mod.iRing])) {
                                    gen_cells_tgt[i] = std::move(gen_cell);
                                    gen_cells_tgt[i].cell += map->deg.stem();
                                    gen_names_tgt[i] = StrGenCell(gen_cells_tgt[i], gen_names_rings[mod.iRing]);
                                    auto old_gen_name = fmt::format("v_{{{},{}{}}}", deg_tgt.stem(), deg_tgt.s, gen_index == 1 ? "" : fmt::format(",{}", gen_index));
                                    fmt::print("({}) {} {} --> {}\n", map->display, name, old_gen_name, gen_names_tgt[i]);
                                    ++count;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*for (size_t iRing = 0; iRing < rings.size(); ++iRing) {
        auto& name = rings[iRing].name;
        auto path = diag_json.at("rings")[iRing].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        ResolveNameConflict(gen_names_rings[iRing]);
        db.begin_transaction();
        db.execute_cmd(fmt::format("UPDATE {}_AdamsE2_generators SET name=NULL", name));
        db.save_gen_names(fmt::format("{}_AdamsE2", name), gen_names_rings[iRing]);
        db.end_transaction();
    }*/
    for (size_t iMod = 0; iMod < mods.size(); ++iMod) {
        auto& name = mods[iMod].name;
        auto& gen_cells = gen_cells_mods[iMod];
        auto path = diag_json.at("modules")[iMod].at("path").get<std::string>();
        MyDB db(db_dir + "/" + path);
        ResolveNameConflict(gen_names_mods[iMod]);
        db.begin_transaction();
        db.execute_cmd(fmt::format("UPDATE {}_AdamsE2_generators SET name=NULL, cell=NULL, cell_coeff=NULL", name));
        db.save_gen_cells(fmt::format("{}_AdamsE2", name), gen_names_mods[iMod], gen_cells_mods[iMod]);
        db.end_transaction();
    }
    return 0;
}

/* Manage generator names */
int main_rename_gen(int argc, char** argv, int& index, const char* desc)
{
    myio::SubCmdArg1d subcmds = {
        {"reset", "reset gen_names of modules", main_rename_gen_reset},
        {"defaultS0", "Set unnamed S0 gen_names by x_{...}", main_rename_gen_defaultS0},
        {"export", "Export current generator name to a json file", main_rename_gen_export},
        {"import", "Import generator name from a json file", main_rename_gen_import},
        {"cell", "Rename generators in zero filtration", main_rename_gen_cell},
        {"pull_back", "Pull-back names along maps", main_rename_gen_pull_back},
        {"push_forward", "push-forward names along maps", main_rename_gen_push_forward},
    };
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, desc, VERSION, subcmds))
        return error;

    return 0;
}
