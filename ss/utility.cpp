#include "main.h"
#include "mylog.h"
#include <regex>

std::string NullDiffCofseq::Str() const
{
    return fmt::format("x={} r={}\nfirst={} count={}\nfirst_ss={} count_ss={}\n", myio::Serialize(x), r, first, count, first_ss, count_ss);
}

int1d Category::Multiply(IndexUniv iCw_x, AdamsDeg deg_x, const int1d& x, IndexUniv iCw_y, AdamsDeg deg_y, const int1d& y) const
{
    int1d xy;
    if (x.empty() || y.empty())
        return xy;
    AdamsDeg deg_xy = deg_x + deg_y;
    if (iCw_x.isRing() && iCw_y.isRing()) {
        auto& ring = rings_[iCw_x.index];
        Poly alg_x = Indices2Poly(x, ring.basis.at(deg_x));
        Poly alg_y = Indices2Poly(y, ring.basis.at(deg_y));
        Poly alg_xy = ring.gb.Reduce(alg_x * alg_y);
        xy = alg_xy ? Poly2Indices(alg_xy, ring.basis.at(deg_xy)) : int1d{};
    }
    else if (iCw_x.isRing() && !iCw_y.isRing()) {
        auto& ring = rings_[iCw_x.index];
        auto& mod = modules_[iCw_y.index];
        Poly alg_x = Indices2Poly(x, ring.basis.at(deg_x));
        Mod alg_y = Indices2Mod(y, mod.basis.at(deg_y));
        Mod alg_xy = mod.gb.Reduce(alg_x * alg_y);
        xy = alg_xy ? Mod2Indices(alg_xy, mod.basis.at(deg_xy)) : int1d{};
    }
    else if (!iCw_x.isRing() && iCw_y.isRing()) {
        auto& mod = modules_[iCw_x.index];
        auto& ring = rings_[iCw_y.index];
        Mod alg_x = Indices2Mod(x, mod.basis.at(deg_x));
        Poly alg_y = Indices2Poly(y, ring.basis.at(deg_y));
        Mod alg_xy = mod.gb.Reduce(alg_y * alg_x);
        xy = alg_xy ? Mod2Indices(alg_xy, mod.basis.at(deg_xy)) : int1d{};
    }
    else
        throw RunTimeError("Cannot multiply modules with modules");
    return xy;
}

int main_mul(int argc, char** argv, int& index, const char* desc)
{
    std::string cw;
    int stem1, s1, stem2, s2;
    std::string str_x1, str_x2;
    int1d x1, x2;

    std::string cat_name;

    myio::CmdArg1d args = {{"cw", &cw}, {"stem1", &stem1}, {"s1", &s1}, {"x1", &str_x1}, {"stem2", &stem2}, {"s2", &s2}, {"x2", &str_x2}, {"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    AdamsDeg d1(s1, stem1 + s1), d2(s2, stem2 + s2), d3(d1 + d2);
    x1 = myio::Deserialize<int1d>(str_x1);
    x2 = myio::Deserialize<int1d>(str_x2);

    SSFlag flag = SSFlag::no_op;
    Category diagram(cat_name, "", flag);

    if (auto iCw = diagram.GetIndexCwByName(cw); iCw.isRing()) {
        auto& ring = diagram.GetRingByName(cw);
        if (d3.t > ring.t_max) {
            fmt::print("degree out of range");
            return 0;
        }
        Poly poly_x1 = Indices2Poly(x1, ring.basis.at(d1));
        Poly poly_x2 = Indices2Poly(x2, ring.basis.at(d2));
        Poly poly_x3 = ring.gb.Reduce(poly_x1 * poly_x2);
        int1d x3 = Poly2Indices(poly_x3, ring.basis.at(d3));
        fmt::print("{} [{}]*[{}]=[{}]\n", d3, myio::Serialize(x1), myio::Serialize(x2), myio::Serialize(x3));
    }
    else {
        auto& mod = diagram.GetModuleByName(cw);
        auto& ring = diagram.GetRings()[mod.iRing];
        if (d3.t > ring.t_max) {
            fmt::print("degree out of range");
            return 0;
        }
        Poly poly_x1 = Indices2Poly(x1, ring.basis.at(d1));
        Mod mod_x2 = Indices2Mod(x2, mod.basis.at(d2));
        Mod mod_x3 = mod.gb.Reduce(poly_x1 * mod_x2);
        int1d x3 = Mod2Indices(mod_x3, mod.basis.at(d3));
        fmt::print("{} [{}]*[{}]=[{}]\n", d3, myio::Serialize(x1), myio::Serialize(x2), myio::Serialize(x3));
    }

    return 0;
}

int main_copy(int argc, char** argv, int& index, const char* desc)
{
    namespace fs = std::filesystem;
    std::string cat_name1, cat_name2;
    bool with_log = false;

    myio::CmdArg1d args = {{"category1", &cat_name1}, {"category2", &cat_name2}};
    myio::CmdArg1d op_args = {{"with_log", &with_log}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    json js_cat = myio::load_json(fmt::format("{}/ss.json", cat_name1));
    fs::create_directory(cat_name2);

    /* copy rings */
    json& json_rings = js_cat.at("rings");
    for (auto& json_ring : json_rings) {
        std::string abs_path1 = fmt::format("{}/{}", cat_name1, json_ring.at("path").get<std::string>());
        std::string abs_path2 = fmt::format("{}/{}", cat_name2, json_ring.at("path").get<std::string>());
        fs::copy(abs_path1, abs_path2);
    }

    /* copy modules */
    json& json_mods = js_cat.at("modules");
    for (auto& json_mod : json_mods) {
        std::string abs_path1 = fmt::format("{}/{}", cat_name1, json_mod.at("path").get<std::string>());
        std::string abs_path2 = fmt::format("{}/{}", cat_name2, json_mod.at("path").get<std::string>());
        fs::copy(abs_path1, abs_path2);
    }

    /* copy maps */
    json& json_maps = js_cat.at("maps");
    for (auto& json_map : json_maps) {
        std::string abs_path1 = fmt::format("{}/{}", cat_name1, json_map.at("path").get<std::string>());
        std::string abs_path2 = fmt::format("{}/{}", cat_name2, json_map.at("path").get<std::string>());
        fs::copy(abs_path1, abs_path2);
    }

    /* copy cofseq */
    json& json_cofseqs = js_cat.at("cofseqs");
    for (auto& json_cof : json_cofseqs) {
        std::string abs_path1 = fmt::format("{}/{}", cat_name1, json_cof.at("path").get<std::string>());
        std::string abs_path2 = fmt::format("{}/{}", cat_name2, json_cof.at("path").get<std::string>());
        fs::copy(abs_path1, abs_path2);
    }

    /* copy ss.json */
    fs::copy(fmt::format("{}/ss.json", cat_name1), fmt::format("{}/ss.json", cat_name2));

    if (with_log) {
        /* copy dbLog.db */
        fs::copy(fmt::format("{}/log.db", cat_name1), fmt::format("{}/log.db", cat_name2));
    }

    return 0;
}

int main_ut_web(int argc, char** argv, int& index, const char* desc)
{
    namespace fs = std::filesystem;
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    json js_cat = myio::load_json(fmt::format("{}/ss.json", cat_name));

    /* copy rings */
    json& json_rings = js_cat.at("rings");
    for (auto& json_ring : json_rings) {
        std::string name = json_ring.at("name");
        fmt::print("[{}](plot.html?diagram=kervaire&data={}),\n", name, name);
    }
    fmt::print("\n");

    /* copy modules */
    json& json_mods = js_cat.at("modules");
    for (const auto& json_mod : json_mods) {
        std::string name = json_mod.at("name");
        fmt::print("[{}](plot.html?diagram=kervaire&data={}),\n", name, name);
    }
    fmt::print("\n");

    /* copy cofseq */
    json& json_cofseqs = js_cat.at("cofseqs");
    for (auto& json_cof : json_cofseqs) {
        std::string name = json_cof.at("name");
        fmt::print("[{}](plot.html?diagram=kervaire&data={}),\n", name, name);
    }
    fmt::print("\n");

    /* copy maps */
    json& json_maps = js_cat.at("maps");
    for (auto& json_map : json_maps) {
        std::string name = json_map.at("name");
        fmt::print("[{}](plot.html?diagram=kervaire&data={}),\n", name, name);
    }

    return 0;
}

constexpr auto MyFmtE2Gen = R"(
echo ".headers on
.mode csv
select id,name,t-s as stem,s from <arg>_AdamsE2_generators;" | sqlite3 kervaire/<arg> > kervaire_csv/<arg>_AdamsE2_generators.csv
)";
MYFMT(FmtE2Gen, MyFmtE2Gen)

constexpr auto MyFmtE2Rel = R"(
echo ".headers on
.mode csv
select rel,t-s as stem,s from <arg>_AdamsE2_relations;" | sqlite3 kervaire/<arg> > kervaire_csv/<arg>_AdamsE2_relations.csv
)";
MYFMT(FmtE2Rel, MyFmtE2Rel)

constexpr auto MyFmtE2Basis = R"(
echo ".headers on
.mode csv
select id - (select min(id) from <arg>_AdamsE2_basis as t2 where t2.s=t1.s and t1.t=t2.t) as 'index', mon, t-s as stem, s, d2 from <arg>_AdamsE2_basis as t1;" | sqlite3 kervaire/<arg> > kervaire_csv/<arg>_AdamsE2_basis.csv
)";
MYFMT(FmtE2Basis, MyFmtE2Basis)

constexpr auto MyFmtSs = R"(
echo ".headers on
.mode csv
select t-s as stem,s,base,diff,level from <arg>_AdamsE2_ss;" | sqlite3 kervaire/<arg> > kervaire_csv/<arg>_AdamsE2_ss.csv
)";
MYFMT(FmtSs, MyFmtSs)

constexpr auto MyFmtMap = R"(
echo ".headers on
.mode csv
select id,map from map_AdamsE2_<arg>;" | sqlite3 kervaire/<arg> > kervaire_csv/map_<arg>.csv
)";
MYFMT(FmtMap, MyFmtMap)

constexpr auto MyFmtCs = R"(
echo ".headers on
.mode csv
select iC,t-s as stem,s,base,diff,level from cofseq_<arg>;" | sqlite3 kervaire/<arg> > kervaire_csv/cofseq_<arg>.csv
)";
MYFMT(FmtCs, MyFmtCs)

int main_ut_export_csv(int argc, char** argv, int& index, const char* desc)
{
    namespace fs = std::filesystem;
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    json js_cat = myio::load_json(fmt::format("{}/ss.json", cat_name));

    /* copy rings */
    json& json_rings = js_cat.at("rings");
    for (auto& js : json_rings) {
        std::string name = js.at("name");
        std::string path = js.at("path");
        /*fmt::print(FmtE2Gen, name, path, name);
        fmt::print(FmtE2Rel, name, path, name);
        fmt::print(FmtE2Basis, name, name, path, name);*/
        fmt::print(FmtSs, name, path, name);
    }
    fmt::print("\n");

    /* copy modules */
    json& json_mods = js_cat.at("modules");
    for (const auto& js : json_mods) {
        std::string name = js.at("name");
        std::string path = js.at("path");
        /*fmt::print(FmtE2Gen, name, path, name);
        fmt::print(FmtE2Rel, name, path, name);
        fmt::print(FmtE2Basis, name, name, path, name);*/
        fmt::print(FmtSs, name, path, name);
    }
    fmt::print("\n");

    /* copy cofseq */
    json& json_cofseqs = js_cat.at("cofseqs");
    for (auto& js : json_cofseqs) {
        std::string name = js.at("name");
        std::string path = js.at("path");
        //fmt::print(FmtCs, name, path, name);
    }
    fmt::print("\n");

    /* copy maps */
    std::regex is_map_regex("^map_AdamsSS_(\\w+?_(?:to|)_\\w+?)(?:_t\\d+|).db$"); /* match example: map_AdamsSS_RP1_4_to_RP3_4_t169.db */
    std::smatch match;
    json& json_maps = js_cat.at("maps");
    for (auto& js : json_maps) {
        std::string name;
        std::string path = js.at("path");
        if (std::regex_search(path, match, is_map_regex); match[0].matched) {
            name = match[1].str();
        }
        else {
            fmt::print("filename={} not supported.\n", path);
            throw ErrorIdMsg(0x839393b2, "File name is not supported.");
        }
        //fmt::print(FmtMap, name, path, name);
    }

    return 0;
}

int main_ut_E2info(int argc, char** argv, int& index, const char* desc)
{
    namespace fs = std::filesystem;
    std::string cat_name;

    myio::CmdArg1d args = {{"category", &cat_name}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Category category(cat_name, "", SSFlag::naming, false);

    /* print rings */
    fmt::print("{}\n{}\n", R"(\begin{tabular}{|l|c|}\hline)", R"(CW spectra & max(t)\\\hline)");
    int count = 0;
    for (auto& ring : category.GetRings()) {
        fmt::print("\\texttt{{{}}} & {}\\\\\\hline\n", std::regex_replace(ring.name, std::regex("_"), "\\_"), ring.t_max);
        if (++count % 45 == 0) {
            fmt::print("{}\n\n", R"(\end{tabular})");
            fmt::print("{}\n{}\n", R"(\begin{tabular}{|l|c|}\hline)", R"(CW spectra & max(t)\\\hline)");
        }
    }

    /* print modules */
    for (auto& mod : category.GetModules()) {
        fmt::print("\\texttt{{{}}} & {}\\\\\\hline\n", std::regex_replace(mod.name, std::regex("_"), "\\_"), mod.t_max);
        if (++count % 45 == 0) {
            fmt::print("{}\n\n", R"(\end{tabular})");
            fmt::print("{}\n{}\n", R"(\begin{tabular}{|l|c|}\hline)", R"(CW spectra & max(t)\\\hline)");
        }
    }
    fmt::print("{}\n\n", R"(\end{tabular})");
    fmt::print("{}\n{}\n", R"(\begin{tabular}{|l|c|}\hline)", R"(map & max(t)\\\hline)");

    /* print maps */
    count = 0;
    for (auto& map : category.GetMaps()) {
        fmt::print("\\texttt{{{}}} & {}\\\\\\hline\n", std::regex_replace(map->name, std::regex("_"), "\\_"), map->t_max);
        if (++count % 45 == 0) {
            fmt::print("{}\n\n", R"(\end{tabular})");
            fmt::print("{}\n{}\n", R"(\begin{tabular}{|l|c|}\hline)", R"(map & max(t)\\\hline)");
        }
    }
    fmt::print("{}\n\n", R"(\end{tabular})");

    return 0;
}

int main_ut_table(int argc, char** argv, int& index, const char* desc)
{
    namespace fs = std::filesystem;
    std::string cat_name, cw;
    int stem, s_min, s_max;

    myio::CmdArg1d args = {{"category", &cat_name}, {"cw", &cw}, {"stem", &stem}, {"s_min", &s_min}, {"s_max", &s_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    Category category(cat_name, "", SSFlag::naming);
    auto iCw = category.GetIndexCwByName(cw);
    fmt::print("{}\n", category.Local2Table(iCw, stem, s_min, s_max));
    

    return 0;
}

/* Deduce differentials and extensions */
int main_ut(int argc, char** argv, int& index, const char* desc)
{
    myio::SubCmdArg1d subcmds = {
        {"web", "Print info for the webpage", main_ut_web},
        //{"table", "Print local table", main_ut_table},
        {"E2info", "Print E2 info", main_ut_E2info},
        {"export_csv", "Export csv files", main_ut_export_csv},
    };
    if (int error = myio::ParseSubCmd(argc, argv, index, PROGRAM, desc, VERSION, subcmds))
        return error;

    return 0;
}