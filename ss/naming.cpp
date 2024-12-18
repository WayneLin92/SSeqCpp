#include "main.h"

void Category::FillGeneratorNames()
{
    std::map<AdamsDeg, int> gen_index;
    for (auto& ring : rings_) {
        for (size_t i = 0; i < ring.gen_degs.size(); ++i) {
            AdamsDeg deg = ring.gen_degs[i];
            ++gen_index[deg];
            if (ring.gen_names[i].empty())
                ring.gen_names[i] = fmt::format("x_{{{},{}{}}}", ring.gen_degs[i].stem(), ring.gen_degs[i].s, gen_index[ring.gen_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[ring.gen_degs[i]]));
        }
    }
    gen_index.clear();
    for (auto& mod : modules_) {
        for (size_t i = 0; i < mod.v_degs.size(); ++i) {
            AdamsDeg deg = mod.v_degs[i];
            ++gen_index[deg];
            if (mod.v_names[i].empty())
                mod.v_names[i] = fmt::format("v_{{{},{}{}}}", mod.v_degs[i].stem(), mod.v_degs[i].s, gen_index[mod.v_degs[i]] == 1 ? "" : fmt::format(",{}", gen_index[mod.v_degs[i]]));
        }
    }
}

std::string Category::GetName(IndexUniv iCw_x, AdamsDeg deg_x, const int1d& x) const
{
    if (x.empty())
        return "0";
    if (iCw_x.isRing()) {
        if (!rings_[iCw_x.index].basis.has(deg_x) || x.back() >= (int)rings_[iCw_x.index].basis.at(deg_x).size())
            throw RunTimeError(fmt::format("GetName Error. iCw_x={}, deg_x={}, x={}", iCw_x, deg_x, myio::Serialize(x)));
        return Str(Indices2Poly(x, rings_[iCw_x.index].basis.at(deg_x)), rings_[iCw_x.index].gen_names);
    }
    else {
        if (!modules_[iCw_x.index].basis.has(deg_x) || x.back() >= (int)modules_[iCw_x.index].basis.at(deg_x).size())
            throw RunTimeError(fmt::format("GetName Error. iCw_x={}, deg_x={}, x={}", iCw_x, deg_x, myio::Serialize(x)));
        return Str(Indices2Mod(x, modules_[iCw_x.index].basis.at(deg_x)), rings_[modules_[iCw_x.index].iRing].gen_names, modules_[iCw_x.index].v_names);
    }
}

json Category::Sc2Json(IndexUniv iCw, AdamsDeg deg) const
{
    json result = {{"name", GetCwName(iCw)}, {"deg", {deg.stem(), deg.s}}, {"sc", json::array()}};

    auto sc = GetNodesSS(iCw).GetSc4Display(deg);
    for (size_t i = 0; i < sc.levels.size(); ++i)
        result.at("sc").push_back(json{{"level", sc.levels[i]}, {"x", GetName(iCw, deg, sc.basis[i])}, {"known_diff", sc.diffs[i] != NULL_DIFF}});
    return result;
}

std::string Category::ScB(IndexUniv iCw, AdamsDeg deg, int level_max) const
{
    auto sc = GetNodesSS(iCw).GetSc4Display(deg);
    std::string result = "\\vspan\\{";
    for (size_t i = 0; i < sc.levels.size(); ++i) {
        if (sc.levels[i] > level_max)
            break;
        if (i > 0)
            result.push_back(',');
        fmt::format_to(std::back_inserter(result), "{}", GetName(iCw, deg, sc.basis[i]));
    }
    result += "\\}";
    if (result == "\\vspan\\{\\}")
        return "0";
    return result;
}

json Category::Local2Json(IndexUniv iCw_x, AdamsDeg deg_x, int r) const
{
    json result = {{"name", GetCwName(iCw_x)}, {"stem", deg_x.stem()}, {"s", deg_x.s}, {"r", r}};
    result["table"] = {json::array(), json::array()};

    auto& nodes_ss = GetNodesSS(iCw_x);
    for (int k = 0; k < 2; ++k) {
        for (int s = 0; s <= r; ++s) {
            auto deg = deg_x + AdamsDeg(s, s - 1 + k);
            result["table"][k].push_back(json::array());
            if (!nodes_ss.front().has(deg))
                continue;
            auto sc = nodes_ss.GetSc4Display(deg);
            auto& js_sc = result["table"][k].back();
            for (size_t i = 0; i < sc.levels.size(); ++i)
                js_sc.push_back(json{{"level", sc.levels[i]}, {"x", GetName(iCw_x, deg, sc.basis[i])}, {"known_diff", sc.diffs[i] != NULL_DIFF}});
        }
    }
    return result;
}

constexpr auto MyFmtLocalTableLastRow = R"(\multicolumn{4}{|c|}{\text{stem}=<arg>} \\\hline
\multicolumn{4}{|c|}{\TitleCellColor $E_2^{*,*}(<arg>$)}\\\hline
)";
MYFMT(FmtLocalTableLastRow, MyFmtLocalTableLastRow)

std::string Category::Local2Table(IndexUniv iCw, int stem, int s_min, int s_max) const
{
    std::string strTable = R"(\begin{center}\scalebox{1.0}{\begin{tabular}{|c|l|l|l|}\hline
$s$ & Elements & $d_r$ & value\\\hline\hline
)";
    auto& nodes_ss = GetNodesSS(iCw);
    auto GetStrS = [](size_t i, size_t n, int s) {
        if (i == 0)
            return fmt::format("\\multirow{{{}}}{{*}}{{{}}}", n, s);
        return std::string("");
    };
    for (int s = s_max; s >= s_min; --s) {
        auto deg = AdamsDeg(s, stem + s);
        if (nodes_ss.front().has(deg)) {
            auto& sc = nodes_ss.GetRecentSc(deg);
            auto n = sc.levels.size();
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                auto hline = i == n - 1 ? R"(\hline\hline)" : R"(\cline{2-4})";
                if (sc.levels[i] > LEVEL_PERM) {
                    int r = LEVEL_MAX - sc.levels[i];
                    fmt::format_to(std::back_inserter(strTable), R"({} & ${}$ & $d_{{{}}}$ & ${}$ \\{})", GetStrS(i, n, s), GetName(iCw, deg, sc.basis[i]), r, sc.diffs[i] == NULL_DIFF ? "?" : GetName(iCw, deg + AdamsDeg(r, r - 1), sc.diffs[i]), hline);
                }
                else if (sc.levels[i] < LEVEL_MAX / 2) {
                    int r = sc.levels[i];
                    fmt::format_to(std::back_inserter(strTable), R"({} & ${}$ & $d_{{{}}}^{{-1}}$ & ${}$ \\{})", GetStrS(i, n, s), GetName(iCw, deg, sc.basis[i]), r, sc.diffs[i] == NULL_DIFF ? "?" : GetName(iCw, deg - AdamsDeg(r, r - 1), sc.diffs[i]),
                                   hline);
                }
                else
                    fmt::format_to(std::back_inserter(strTable), R"({} & ${}$ & & Permanent \\{})", GetStrS(i, n, s), GetName(iCw, deg, sc.basis[i]), hline);
                strTable.push_back('\n');
            }
        }
        else {
            fmt::format_to(std::back_inserter(strTable), R"({} & \multicolumn{{3}}{{c|}}{{}}\\\hline\hline{})", GetStrS(0, 1, s), "\n");
        }
    }
    fmt::format_to(std::back_inserter(strTable), FmtLocalTableLastRow, stem, GetCwName(iCw));
    fmt::format_to(std::back_inserter(strTable), R"(\end{{tabular}}}}\end{{center}}{})", "\n");
    return strTable;
}

std::tuple<std::string, std::string, std::string, std::string> Category::GetMapName(IndexUniv iCs) const
{
    auto& cofseq = cofseqs_[iCs.index];
    auto map_letter = "fgh";
    auto cofseq_name = fmt::format(R"({} \fto{{f}} {} \fto{{g}} {} \fto{{h}} {})", cofseq.nameCw[0], cofseq.nameCw[1], cofseq.nameCw[2], cofseq.nameCw[0]);
    return std::make_tuple(fmt::format("{}", map_letter[iCs.iTri]), cofseq.nameCw[iCs.iTri], cofseq.nameCw[NextiTri(iCs.iTri)], cofseq_name);
}

json Category::Sc2JsonCofseq(IndexUniv iCs, AdamsDeg deg) const
{
    auto& cofseq = cofseqs_[iCs.index];
    json result = {{"name", cofseq.nameCw[iCs.iTri]}, {"deg", {deg.stem(), deg.s}}, {"sc", json::array()}};
    auto& js_sc = result.at("sc");
    auto sc = cofseq.nodes_ss[iCs.iTri]->GetSc4Display(deg);
    IndexUniv iCw = cofseq.indexCw[iCs.iTri];

    for (size_t i = 0; i < sc.levels.size(); ++i) {
        if (sc.levels[i] >= LEVEL_PERM)
            break;
        js_sc.push_back(json{{"level", sc.levels[i]}, {"x", GetName(iCw, deg, sc.basis[i])}, {"known_diff", false}});
    }

    auto sc_cs = GetSc4Display(cofseq.nodes_cofseq[iCs.iTri], cofseq.nodes_cofseq[PreviTri(iCs.iTri)], deg, cofseq.degMap[PreviTri(iCs.iTri)].stem());
    for (size_t i = 0; i < sc_cs.levels.size(); ++i)
        js_sc.push_back(json{{"level", sc_cs.levels[i] + 2 * LEVEL_MAX}, {"x", GetName(iCw, deg, sc_cs.basis[i])}, {"known_diff", sc_cs.diffs[i] != NULL_DIFF}});

    for (size_t i = 0; i < sc.levels.size(); ++i) {
        if (sc.levels[i] > LEVEL_PERM)
            js_sc.push_back(json{{"level", sc.levels[i]}, {"x", GetName(iCw, deg, sc.basis[i])}, {"known_diff", false}});
    }
    return result;
}

std::string Category::ScBCofseq(IndexUniv iCs, AdamsDeg deg, int level_max) const
{
    auto& cofseq = cofseqs_[iCs.index];
    IndexUniv iCw = cofseq.indexCw[iCs.iTri];
    std::string result = "\\vspan\\{";
    {
        auto sc_ss = cofseq.nodes_ss[iCs.iTri]->GetSc4Display(deg);
        for (size_t i = 0; i < sc_ss.levels.size(); ++i) {
            if (sc_ss.levels[i] > LEVEL_MAX / 2)
                break;
            if (i > 0)
                result.push_back(',');
            fmt::format_to(std::back_inserter(result), "{}", GetName(iCw, deg, sc_ss.basis[i]));
        }
    }
    {
        auto sc_cs = GetSc4Display(cofseq.nodes_cofseq[iCs.iTri], cofseq.nodes_cofseq[PreviTri(iCs.iTri)], deg, cofseq.degMap[PreviTri(iCs.iTri)].stem());
        for (size_t i = 0; i < sc_cs.levels.size(); ++i) {
            if (i > 0)
                result.push_back(',');
            if (sc_cs.levels[i] > level_max)
                break;
            fmt::format_to(std::back_inserter(result), "{}", GetName(iCw, deg, sc_cs.basis[i]));
        }
    }
    result += "\\}";
    if (result == "\\vspan\\{\\}")
        return "0";
    return result;
}

json Category::LocalCofseqToJson(IndexUniv iCs_x, AdamsDeg deg_x, AdamsDeg deg_dx, int r) const
{
    auto [map_letter, name1, name2, _] = GetMapName(iCs_x);
    json result = {{"name", {map_letter, name1, name2}}, {"deg_x", {deg_x.stem(), deg_x.s}}, {"deg_dx", {deg_dx.stem(), deg_dx.s}}};
    result["table"] = {json::array(), json::array()};

    auto& cofseq = cofseqs_[iCs_x.index];
    auto nodes_sses = std::array{cofseq.nodes_ss[iCs_x.iTri], cofseq.nodes_ss[NextiTri(iCs_x.iTri)]};

    auto iCw_x = cofseq.indexCw[iCs_x.iTri];
    auto iCw_dx = cofseq.indexCw[NextiTri(iCs_x.iTri)];
    auto iCws = std::array{iCw_x, iCw_dx};
    auto degs = std::array{deg_x, deg_dx};
    auto iCses = std::array{iCs_x, iCs_x.next()};
    for (int k = 0; k < 2; ++k) {
        for (int s = 0; s <= r; ++s) {
            auto iCw = iCws[k];
            auto deg = AdamsDeg(deg_x.s, degs[k].stem() + deg_x.s) + AdamsDeg(s, s);
            result["table"][k].push_back(json::array());
            if (!nodes_sses[k]->front().has(deg))
                continue;

            auto sc = nodes_sses[k]->GetSc4Display(deg);
            auto& js_sc = result["table"][k].back();
            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] >= LEVEL_PERM)
                    break;
                js_sc.push_back(json{{"level", sc.levels[i]}, {"x", GetName(iCw, deg, sc.basis[i])}, {"known_diff", false}});
            }

            auto sc_cs = GetSc4Display(cofseq.nodes_cofseq[iCses[k].iTri], cofseq.nodes_cofseq[PreviTri(iCses[k].iTri)], deg, cofseq.degMap[PreviTri(iCses[k].iTri)].stem());
            for (size_t i = 0; i < sc_cs.levels.size(); ++i)
                js_sc.push_back(json{{"level", sc_cs.levels[i] + 2 * LEVEL_MAX}, {"x", GetName(iCw, deg, sc_cs.basis[i])}, {"known_diff", sc_cs.diffs[i] != NULL_DIFF}});

            for (size_t i = 0; i < sc.levels.size(); ++i) {
                if (sc.levels[i] > LEVEL_PERM)
                    js_sc.push_back(json{{"level", sc.levels[i]}, {"x", GetName(iCw, deg, sc.basis[i])}, {"known_diff", false}});
            }
        }
    }
    return result;
}

void print_to(std::string& out, const Mon& m, const std::vector<std::string>& gen_names)
{
    if (m.size() == 0) {
        out.push_back('1');
        return;
    }
    for (size_t i = 0; i < m.size(); ++i) {
        fmt::format_to(std::back_inserter(out), "{}", gen_names[m[i].g()]);
        if (m[i].e() >= 10)
            fmt::format_to(std::back_inserter(out), "^{{{}}}", m[i].e());
        else if (m[i].e() >= 2)
            fmt::format_to(std::back_inserter(out), "^{}", m[i].e());
    }
}

void print_to(std::string& out, const MMod& m, const std::vector<std::string>& gen_names, const std::vector<std::string>& v_names)
{
    if (m.m.size())
        print_to(out, m.m, gen_names);
    fmt::format_to(std::back_inserter(out), "{}", v_names[m.v]);
}

void print_to(std::string& out, const Poly& p, const std::vector<std::string>& gen_names)
{
    if (p.data.size() == 0) {
        out.push_back('0');
        return;
    }
    for (size_t i = 0; i < p.data.size(); ++i) {
        print_to(out, p.data[i], gen_names);
        if (i + 1 < p.data.size())
            out.push_back('+');
    }
}

void print_to(std::string& out, const Mod& p, const std::vector<std::string>& gen_names, const std::vector<std::string>& v_names)
{
    if (p.data.size() == 0) {
        out.push_back('0');
        return;
    }
    for (size_t i = 0; i < p.data.size(); ++i) {
        print_to(out, p.data[i], gen_names, v_names);
        if (i + 1 < p.data.size())
            out.push_back('+');
    }
}

std::string Str(const Mon& m, const std::vector<std::string>& gen_names)
{
    std::string result;
    print_to(result, m, gen_names);
    return result;
}

std::string Str(const Poly& p, const std::vector<std::string>& gen_names)
{
    std::string result;
    print_to(result, p, gen_names);
    return result;
}

std::string Str(const MMod& m, const std::vector<std::string>& gen_names, const std::vector<std::string>& v_names)
{
    std::string result;
    print_to(result, m, gen_names, v_names);
    return result;
}

std::string Str(const Mod& p, const std::vector<std::string>& gen_names, const std::vector<std::string>& v_names)
{
    std::string result;
    print_to(result, p, gen_names, v_names);
    return result;
}