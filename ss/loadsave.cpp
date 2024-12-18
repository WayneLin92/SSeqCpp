/* @file config.cpp
 * Load settings from json
 */

#include "main.h"
#include "mylog.h"
#include <fstream>
#include <regex>

void LoadJson(const std::string& cat_name, nlohmann::json& root_json, nlohmann::json& diag_json)
{
    using json = nlohmann::json;

    root_json = myio::load_json("ss.json");
    try {
        diag_json = myio::load_json(fmt::format("{}/ss.json", cat_name));

        auto path_log = fmt::format("{}/log.db", cat_name);
        Logger::SetOutDeduce(path_log.c_str());
    }
    catch (nlohmann::detail::exception& e) {
        fmt::print("JsonError({}): {}\n", e.id, e.what());
        throw e;
    }
}

void GetAllDbNames(const std::string& cat_name, std::vector<std::string>& names, std::vector<std::string>& paths, std::vector<int>& isRing, bool log)
{
    using json = nlohmann::json;

    try {
        json js_cat = myio::load_json(fmt::format("{}/ss.json", cat_name));
        if (log) {
            auto path_log = fmt::format("{}/{}", cat_name, js_cat.at("log").get<std::string>());
            Logger::SetOutDeduce(path_log.c_str());
        }

        /* #Load rings */
        json& json_rings = js_cat.at("rings");
        for (auto& json_ring : json_rings) {
            std::string name = json_ring.at("name").get<std::string>(), path = json_ring.at("path").get<std::string>();
            std::string abs_path = fmt::format("{}/{}", cat_name, path);
            names.push_back(name);
            paths.push_back(abs_path);
            isRing.push_back(1);
        }

        /* #Load modules */
        json& json_mods = js_cat.at("modules");
        for (auto& json_mod : json_mods) {
            std::string name = json_mod.at("name").get<std::string>(), path = json_mod.at("path").get<std::string>();
            std::string abs_path = fmt::format("{}/{}", cat_name, path);
            names.push_back(name);
            paths.push_back(abs_path);
            isRing.push_back(0);
        }
    }
    catch (nlohmann::detail::exception& e) {
        fmt::print("JsonError({}): {}\n", e.id, e.what());
        throw e;
    }
}

void get_composition_info(const nlohmann::json& json_maps, const PMap1d& maps, const int1d& indices, const std::string& from, const std::string& to, int t_max_from, int& t_max, AdamsDeg& deg_composition)
{
    t_max = t_max_from;
    int t = t_max;
    for (size_t i = 0; i < indices.size(); ++i) {  // Check if the composition is well-defined
        auto& json_f = json_maps[indices[i]];
        auto& map = maps[indices[i]];
        if (t > map->t_max) {
            t_max += map->t_max - t;
            t = map->t_max;
        }
        t += map->deg.t;
        deg_composition += map->deg;
        if (i == 0)
            ErrorIdMsg::Assert(json_f.at("from").get<std::string>() == from, "Valid composition");
        if (i + 1 == indices.size())
            ErrorIdMsg::Assert(json_f.at("to").get<std::string>() == to, "Valid composition");
        else
            ErrorIdMsg::Assert(json_f.at("to").get<std::string>() == json_maps[indices[i + 1]].at("from").get<std::string>(), "Valid composition");
    }
}

int1d get_compostion(const int1d& x, AdamsDeg deg, const Category& category, const PMap1d& maps, const int1d& indices)
{
    int1d result = x;
    for (size_t i = 0; i < indices.size(); ++i) {
        result = maps[indices[i]]->map(result, deg, category);
        deg += maps[indices[i]]->deg;
    }
    return result;
}

void Category::LoadExclusions(int id_start, SSFlag flag)
{
    std::regex is_cofseq_regex("^(\\w+?__\\w+?__\\w+?):([012])$"); /* match example: X__Y__Z:1 */
    std::smatch match;
    myio::Statement stmt(Logger::db_, fmt::format("SELECT name, s, t, r, x, dx, reason like '%' || 'I' FROM log, exclusions where log.id=exclusions.id and log.id>{}", id_start));
    while (stmt.step() == MYSQLITE_ROW) {
        auto name = stmt.column_str(0);
        auto s = stmt.column_int(1);
        auto t = stmt.column_int(2);
        auto r = stmt.column_int(3);
        int1d x = myio::Deserialize<int1d>(stmt.column_str(4));
        int1d dx = myio::Deserialize<int1d>(stmt.column_str(5));
        auto isDInv = stmt.column_int(6);
        AdamsDeg deg = AdamsDeg(s, t);
        if (isDInv)
            std::swap(x, dx);
        if (std::regex_search(name, match, is_cofseq_regex); match[0].matched) {
            if (flag & SSFlag::cofseq) {
                auto iCof = GetCofSeqIndexByName(match[1].str());
                auto iTri = (size_t)std::stoi(match[2].str());
                auto& nodes_cs = cofseqs_[iCof].nodes_cofseq[iTri];
                auto& e = nodes_cs.exclusions.get(deg, x, !isDInv ? r : -r - 1);
                e.dxs.push_back(std::move(dx));
            }
        }
        else {
            auto iCw = GetIndexCwByName(name);
            auto& nodes_ss = GetNodesSS(iCw);
            auto& e = nodes_ss.exclusions.get(deg, x, !isDInv ? r : -r);
            e.dxs.push_back(std::move(dx));
        }
    }

    for (auto& ring : rings_)
        ring.nodes_ss.exclusions.Sort();
    for (auto& mod : modules_)
        mod.nodes_ss.exclusions.Sort();
    for (auto& cofseq : cofseqs_)
        for (auto& nodes : cofseq.nodes_cofseq)
            nodes.exclusions.Sort();
}

Category::Category(const std::string& cat_root, const std::string& ckpt_name, SSFlag flags, bool log, bool loadD2)
{
    try {  //// TODO: check files in the beginning
        std::string cat_nodes;
        if (ckpt_name.size())
            cat_nodes = fmt::format("{}/checkpoints/{}", cat_root, ckpt_name);

        js_ = myio::load_json(fmt::format("{}/ss.json", cat_root));
        if (log) {
            auto path_log = fmt::format("{}/log.db", cat_root);
            Logger::SetOutDeduce(path_log.c_str());
        }

        /*# Load rings and modules */

        auto& json_rings = js_.at("rings");
        auto& json_mods = js_.at("modules");
        AdamsDeg2d ring_gen_degs;
        AdamsDeg2d module_gen_degs;  //// TODO: Remove module_gen_degs
        rings_.reserve(json_rings.size());
        modules_.reserve(json_mods.size());
        std::vector<std::map<AdamsDeg, int2d>> basis_d2;
        {
            size_t iCw = 0;
            for (auto& json_ring : json_rings) {
                std::string name = json_ring.at("name").get<std::string>(), path = json_ring.at("path").get<std::string>();
                if (myio::get(json_ring, "deduce", "on") == "on")
                    deduce_list_spectra_.push_back(IndexRing(iCw));
                std::string abs_path = fmt::format("{}/{}", cat_root, path);
                std::string table_prefix = fmt::format("{}_AdamsE2", name);

                myio::AssertFileExists(abs_path);
                DBSS db(abs_path);
                DBSS db2;
                if (cat_nodes.size()) {
                    std::string abs_path_nodes = fmt::format("{}/{}", cat_nodes, path);
                    myio::AssertFileExists(abs_path_nodes);
                    db2.open(abs_path_nodes);
                }
                auto& db_nodes = cat_nodes.size() ? db2 : db;
                rings_.push_back({});
                auto& ring = rings_.back();
                ring.name = name;
                if (flags & SSFlag::naming) {
                    ring.gen_degs = db.load_gen_adamsdegs(table_prefix);
                    ring.gen_names = db.load_gen_names(table_prefix);
                }
                ring.basis = db.load_basis(table_prefix);
                if (loadD2)
                    basis_d2.push_back(db.load_basis_d2(table_prefix));
                ring.t_max = db.get_metadata_int("t_max");
                if (ring.t_max <= 0)
                    throw RunTimeError(fmt::format("No t_max. name={}.", ring.name));
                ring.nodes_ss.push_back(db_nodes.load_ss(table_prefix));
                ring.nodes_ss.push_back({});
                ring.degs_ss = ring.nodes_ss.front().arr_degs();
                ring.gb = Groebner(ring.t_max, {}, db.load_gb(table_prefix, DEG_MAX));

                if (flags & SSFlag::pi) {
                    ring.pi_gen_Einf = db.get_column_from_str<Poly>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Poly>);
                    ring.pi_gb = algZ::Groebner(ring.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb(name, DEG_MAX));
                    ring.nodes_pi_basis.reserve(MAX_DEPTH + 1);
                    if (flags & SSFlag::pi_def)
                        db_nodes.load_pi_def(name, ring.pi_gen_defs, ring.pi_gen_def_mons);
                }

                ring_gen_degs.push_back(db.load_gen_adamsdegs(table_prefix));
                ++iCw;
            }

            iCw = 0;
            for (auto& json_mod : json_mods) {
                std::string name = json_mod.at("name").get<std::string>(), path = json_mod.at("path").get<std::string>();
                std::string over = json_mod.at("over").get<std::string>();
                if (myio::get(json_mod, "deduce", "on") == "on")
                    deduce_list_spectra_.push_back(IndexMod(iCw));
                std::string abs_path = fmt::format("{}/{}", cat_root, path);
                std::string table_prefix = fmt::format("{}_AdamsE2", name);

                myio::AssertFileExists(abs_path);
                DBSS db(abs_path);
                DBSS db2;
                if (cat_nodes.size()) {
                    std::string abs_path_nodes = fmt::format("{}/{}", cat_nodes, path);
                    myio::AssertFileExists(abs_path_nodes);
                    db2.open(abs_path_nodes);
                }
                auto& db_nodes = cat_nodes.size() ? db2 : db;
                modules_.push_back({});
                auto& mod = modules_.back();
                mod.name = name;
                if (flags & SSFlag::naming) {
                    mod.v_degs = db.load_gen_adamsdegs(table_prefix);
                    mod.v_names = db.load_gen_names(table_prefix);
                }
                auto indexRing = GetIndexCwByName(over);
                ErrorIdMsg::Assert(indexRing.isRing(), "indexRing.isRing()");
                mod.iRing = indexRing.index;
                auto& ring = rings_[mod.iRing];
                mod.basis = db.load_basis_mod(table_prefix);
                if (loadD2)
                    basis_d2.push_back(db.load_basis_d2(table_prefix));
                mod.t_max = db.get_metadata_int("t_max");
                if (mod.t_max <= 0)
                    throw RunTimeError(fmt::format("No t_max. name={}.", mod.name));
                mod.nodes_ss.push_back(db_nodes.load_ss(table_prefix));
                mod.nodes_ss.push_back({});
                mod.degs_ss = mod.nodes_ss.front().arr_degs();
                Mod1d xs = db.load_gb_mod(table_prefix, DEG_MAX);
                mod.gb = GroebnerMod(&ring.gb, mod.t_max, {}, std::move(xs));

                if (flags & SSFlag::pi) {
                    mod.pi_gen_Einf = db.get_column_from_str<Mod>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Mod>);
                    mod.pi_gb = algZ::GroebnerMod(&ring.pi_gb, mod.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb_mod(name, DEG_MAX));
                    mod.nodes_pi_basis.reserve(MAX_DEPTH + 1);
                    if (flags & SSFlag::pi_def)
                        db_nodes.load_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
                }

                module_gen_degs.push_back(db.load_gen_adamsdegs(table_prefix));
                ring.ind_mods.push_back(modules_.size() - 1);
                ++iCw;
            }
        }

        /*# Load maps */
        std::regex is_map_regex("^map_AdamsSS_(\\w+?_(?:to|)_\\w+?)(?:_t\\d+|).db$"); /* match example: map_AdamsSS_RP1_4_to_RP3_4_t169.db */
        auto& json_maps = js_.at("maps");
        if (flags & SSFlag::cofseq || flags & SSFlag::naming) {
            for (auto& item : js_.at("maps_v2"))
                json_maps.push_back(item);
        }
        maps_.reserve(json_maps.size());
        for (auto& json_map : json_maps) {
            std::string name = json_map.at("name");
            auto display = name;
            auto str_pos = display.find("__", 0, 2);
            display.replace(str_pos, 2, " -> ");
            std::string from = json_map.at("from"), to = json_map.at("to");
            std::unique_ptr<Map> map;

            if (json_map.contains("factor")) {
                auto& strt_factor = json_map.at("factor");
                int stem = strt_factor[0], s = strt_factor[1];
                int1d factor = strt_factor[2].is_number() ? int1d{strt_factor[2].get<int>()} : strt_factor[2].get<int1d>();
                AdamsDeg deg_factor(s, stem + s);
                if (auto ifrom = GetIndexCwByName(from); ifrom.isRing()) {
                    if (auto ito = GetIndexCwByName(to); ito.isRing()) {
                        ErrorIdMsg::Assert(ifrom == ito, "index_from == index_to");
                        Poly poly_factor = factor.size() ? Indices2Poly(factor, rings_[ifrom.index].basis.at(deg_factor)) : Poly();
                        int t_max = rings_[ito.index].t_max - deg_factor.t;
                        map = std::make_unique<MapMulRing2Ring>(name, display, t_max, deg_factor, ifrom.index, std::move(poly_factor));
                    }
                    else {
                        Mod mod_factor = factor.size() ? Indices2Mod(factor, modules_[ito.index].basis.at(deg_factor)) : Mod();
                        int t_max = modules_[ito.index].t_max - deg_factor.t;
                        map = std::make_unique<MapMulRing2Mod>(name, display, t_max, deg_factor, ifrom.index, ito.index, std::move(mod_factor));
                    }
                }
                else {
                    auto ito = GetIndexCwByName(to);
                    ErrorIdMsg::Assert(ifrom == ito, "index_from == index_to");
                    Poly poly_factor = factor.size() ? Indices2Poly(factor, rings_[modules_[ifrom.index].iRing].basis.at(deg_factor)) : Poly();
                    int t_max = modules_[ito.index].t_max - deg_factor.t;
                    map = std::make_unique<MapMulMod2Mod>(name, display, t_max, deg_factor, ifrom.index, std::move(poly_factor));
                }
            }
            else if (json_map.contains("composition")) {
                auto& json_composition = json_map.at("composition");
                int1d indices;
                for (auto& json_f : json_composition) {
                    auto name_map = json_f.get<std::string>();
                    int index_map = GetMapIndexByName(name_map);
                    ErrorIdMsg::Assert(index_map != -1, fmt::format("index_map({}) != -1", name_map));
                    indices.push_back(index_map);
                }
                int t_max_map;
                AdamsDeg deg_map;
                if (auto ifrom = GetIndexCwByName(from); ifrom.isRing()) {
                    auto ito = GetIndexCwByName(to);
                    ErrorIdMsg::Assert(ito.isRing(), fmt::format("{}: ito.isRing()", to));

                    auto& gen_degs = ring_gen_degs[ifrom.index];
                    int2d generators;
                    for (size_t i = 0; i < gen_degs.size(); ++i) {
                        int index = ut::IndexOf(rings_[ifrom.index].basis.at(gen_degs[i]), Mon::Gen((uint32_t)i));
                        generators.push_back({index});
                    }
                    get_composition_info(json_maps, maps_, indices, from, to, rings_[ifrom.index].t_max, t_max_map, deg_map);

                    Poly1d images_map;
                    for (size_t i = 0; i < gen_degs.size(); ++i) {
                        if (gen_degs[i].t > t_max_map)
                            break;
                        int1d fx = get_compostion(generators[i], gen_degs[i], *this, maps_, indices);
                        images_map.push_back(Indices2Poly(fx, rings_[ito.index].basis.at(gen_degs[i] + deg_map)));
                    }
                    map = std::make_unique<MapRing2Ring>(name, display, t_max_map, deg_map, ifrom.index, ito.index, std::move(images_map));
                }
                else {
                    auto& gen_degs = module_gen_degs[ifrom.index];
                    int2d generators;
                    for (size_t i = 0; i < gen_degs.size(); ++i) {
                        int index = ut::IndexOf(modules_[ifrom.index].basis.at(gen_degs[i]), MMod({}, (uint32_t)i));
                        generators.push_back({index});
                    }
                    get_composition_info(json_maps, maps_, indices, from, to, modules_[ifrom.index].t_max, t_max_map, deg_map);

                    if (auto ito = GetIndexCwByName(to); ito.isRing()) {
                        Poly1d images_map;
                        for (size_t i = 0; i < gen_degs.size(); ++i) {
                            if (gen_degs[i].t > t_max_map)
                                images_map.push_back(Mon::Null());
                            else {
                                int1d fx = get_compostion(generators[i], gen_degs[i], *this, maps_, indices);
                                if (fx.empty())
                                    images_map.push_back({});
                                else
                                    images_map.push_back(Indices2Poly(fx, rings_[ito.index].basis.at(gen_degs[i] + deg_map)));
                            }
                        }
                        map = std::make_unique<MapMod2Ring>(name, display, t_max_map, deg_map, ifrom.index, ito.index, std::move(images_map));
                    }
                    else {
                        Mod1d images_map;
                        for (size_t i = 0; i < gen_degs.size(); ++i) {
                            if (gen_degs[i].t > t_max_map)
                                images_map.push_back(MMod::Null());
                            else {
                                int1d fx = get_compostion(generators[i], gen_degs[i], *this, maps_, indices);
                                if (fx.empty())
                                    images_map.push_back({});
                                else
                                    images_map.push_back(Indices2Mod(fx, modules_[ito.index].basis.at(gen_degs[i] + deg_map)));
                            }
                        }
                        map = std::make_unique<MapMod2Mod>(name, display, t_max_map, deg_map, ifrom.index, ito.index, std::move(images_map));
                    }
                }
            }
            else {
                int fil = json_map.contains("fil") ? json_map["fil"].get<int>() : 0;
                int sus = json_map.contains("sus") ? json_map["sus"].get<int>() : 0;
                int t_max = json_map["t_max"].get<int>();
                AdamsDeg deg(fil, fil - sus);

                std::string path = json_map.at("path").get<std::string>();
                std::string abs_path = fmt::format("{}/{}", cat_root, path);
                myio::AssertFileExists(abs_path);
                DBSS db(abs_path);
                std::string table;
                std::smatch match;
                if (std::regex_search(path, match, is_map_regex); match[0].matched) {
                    table = fmt::format("map_AdamsE2_{}", match[1].str());
                }
                else {
                    fmt::print("filename={} not supported.\n", path);
                    throw ErrorIdMsg(0x839393b2, "File name is not supported.");
                }

                if (auto ifrom = GetIndexCwByName(from); ifrom.isRing()) {
                    auto ito = GetIndexCwByName(to);
                    ErrorIdMsg::Assert(ifrom && ifrom.isRing(), fmt::format("Need ifrom && ifrom.isRing(). from={}", from));
                    ErrorIdMsg::Assert(ito && ito.isRing(), fmt::format("Need ito && ito.isRing(). to={}", to));
                    auto images = db.get_column_from_str<Poly>(table, "map", "ORDER BY id", myio::Deserialize<Poly>);
                    map = std::make_unique<MapRing2Ring>(name, display, t_max, deg, ifrom.index, ito.index, std::move(images));
                    rings_[ifrom.index].ind_maps.push_back(maps_.size());
                    rings_[ito.index].ind_maps_prev.push_back(maps_.size());
                }
                else {
                    if (auto ito = GetIndexCwByName(to); ito.isRing()) {
                        ErrorIdMsg::Assert(bool(ifrom), fmt::format("Cannot find from=", from));
                        ErrorIdMsg::Assert(bool(ito), fmt::format("Cannot find to={}", to));
                        auto images = db.get_column_from_str<Poly>(table, "map", "", myio::Deserialize<Poly>);
                        if (!json_map.contains("over")) {
                            map = std::make_unique<MapMod2Ring>(name, display, t_max, deg, ifrom.index, ito.index, std::move(images));
                            //((MapMod2Ring&)(*map)).Verify(*this, ring_gen_degs);  ////
                        }
                        else {
                            std::string over = json_map.at("over").get<std::string>();
                            int index_map = ut::IndexOf(maps_, [&over](const std::unique_ptr<Map>& map) { return map->name == over; });
                            map = std::make_unique<MapMod2RingV2>(name, display, t_max, deg, ifrom.index, ito.index, index_map, std::move(images));
                            //((MapMod2RingV2&)(*map)).Verify(*this, ring_gen_degs);  ////
                        }
                        rings_[ito.index].ind_maps_prev.push_back(maps_.size());
                    }
                    else {
                        auto images = db.get_column_from_str<Mod>(table, "map", "", myio::Deserialize<Mod>);
                        if (!json_map.contains("over")) {
                            map = std::make_unique<MapMod2Mod>(name, display, t_max, deg, ifrom.index, ito.index, std::move(images));
                            //((MapMod2Mod&)(*map)).Verify(*this, ring_gen_degs);  ////
                        }
                        else {
                            std::string over = json_map.at("over").get<std::string>();
                            int index_map = ut::IndexOf(maps_, [&over](const std::unique_ptr<Map>& map) { return map->name == over; });
                            map = std::make_unique<MapMod2ModV2>(name, display, t_max, deg, ifrom.index, ito.index, index_map, std::move(images));
                            //((MapMod2ModV2&)(*map)).Verify(*this, ring_gen_degs);  ////
                        }
                        modules_[ito.index].ind_maps_prev.push_back(maps_.size());
                    }
                    modules_[ifrom.index].ind_maps.push_back(maps_.size());
                }
            }
            maps_.push_back(std::move(map));
        }

        for (auto& json_map : json_maps) {
            std::string from = json_map.at("from");
            if (from == "S0" && json_map.contains("factor")) {
                std::string to = json_map.at("to");
                auto& strt_factor = json_map.at("factor");
                int stem = strt_factor[0], s = strt_factor[1];
                int1d factor = strt_factor[2].is_number() ? int1d{strt_factor[2].get<int>()} : strt_factor[2].get<int1d>();
                AdamsDeg deg_factor(s, stem + s);
                if (auto rt = SetCwDiffGlobal(GetIndexCwByName(to), deg_factor, factor, int1d{}, R_PERM - 1, false, flags))
                    throw RunTimeError(fmt::format("Failed to SetCwDiffGlobal(). map={}.", json_map.at("name").get<std::string>()));
            }
        }

        /*# Load cofseqs */
        if (flags & SSFlag::cofseq) {
            auto& json_cofseqs = js_.at("cofseqs");
            std::vector<std::string> cofseq_maps;
            size_t iCof = 0;
            for (auto& json_cofseq : json_cofseqs) {
                CofSeq cofseq;
                cofseq.name = json_cofseq.at("name").get<std::string>();
                auto path_cofseq = json_cofseq.at("path").get<std::string>();
                auto abs_path_cofseq = fmt::format("{}/{}", cat_nodes.size() ? cat_nodes : cat_root, path_cofseq);
                std::array<std::string, 3> maps_names = {json_cofseq.at("i").get<std::string>(), json_cofseq.at("q").get<std::string>(), json_cofseq.at("d").get<std::string>()};
                if (!json_cofseq.contains("deduce") || json_cofseq.at("deduce").get<std::string>() == "on")
                    deduce_list_cofseq_.push_back(iCof);
                for (size_t iTri = 0; iTri < 3; ++iTri) {
                    size_t index = (size_t)GetMapIndexByName(maps_names[iTri]);
                    ErrorIdMsg::Assert(index != -1, maps_names[iTri] + " not found");
                    const auto& map = maps_[index];
                    cofseq.indexMap[iTri] = index;
                    map->iCof = IndexCof(iCof, iTri);
                    cofseq.degMap[iTri] = map->deg;
                    cofseq.indexCw[iTri] = map->from;
                    GetIndexCofs(cofseq.indexCw[iTri]).push_back(IndexCof(iCof, iTri));
                    cofseq.nameCw[iTri] = GetCwName(cofseq.indexCw[iTri]);
                    cofseq.t_max[iTri] = GetTMax(cofseq.indexCw[iTri]);
                    cofseq.nodes_ss[iTri] = &GetNodesSS(cofseq.indexCw[iTri]);
                    cofseq.nodes_cofseq[iTri].push_back({});
                    cofseq.nodes_cofseq[iTri].push_back({});
                }
                if (myio::FileExists(abs_path_cofseq)) { /* Load cofseq from database */
                    DBSS db(abs_path_cofseq);
                    db.create_cofseq(fmt::format("cofseq_{}", cofseq.name));
                    auto node_cofseq = db.load_cofseq(fmt::format("cofseq_{}", cofseq.name));
                    for (size_t i = 0; i < cofseq.nodes_cofseq.size(); ++i)
                        cofseq.nodes_cofseq[i].front() = std::move(node_cofseq[i]);
                }
                cofseqs_.push_back(cofseq);
                ++iCof;
            }
            if (auto rt = SyncCofseq(flags))
                throw RunTimeError(fmt::format("Failed to SyncCofseq(). Error={}.\n", rt.code));
        }

        /*# Load exclusions */
        if (log && !(flags & SSFlag::no_exclusions))
            LoadExclusions(0, flags);

        /*# Load commutativities */
        if (flags & SSFlag::cofseq && js_.contains("commutativity")) {
            auto& json_comms = js_.at("commutativity");
            for (auto& json_comm : json_comms) {
                auto f0 = json_comm.at("f")[0].get<std::string>();
                auto f1 = json_comm.at("f")[1].get<std::string>();
                int if0 = GetMapIndexByName(f0);
                int if1 = GetMapIndexByName(f1);
                ErrorIdMsg::Assert(if0 >= 0, "if0 >= 0");
                ErrorIdMsg::Assert(if1 >= 0, "if1 >= 0");
                ErrorIdMsg::Assert(maps_[if0]->to == maps_[if1]->from, "f0->to == f1->from");
                ErrorIdMsg::Assert(maps_[if0]->iCof.index != -1, fmt::format("{}->iCof.index != -1", maps_[if0]->name));
                ErrorIdMsg::Assert(maps_[if1]->iCof.index != -1, fmt::format("{}->iCof.index != -1", maps_[if1]->name));
                int ig0 = -1, ig1 = -1;
                if (json_comm.at("g").size() >= 1) {
                    auto g0 = json_comm.at("g")[0].get<std::string>();
                    ig0 = GetMapIndexByName(g0);
                    ErrorIdMsg::Assert(ig0 >= 0, "ig0 >= 0");
                    ErrorIdMsg::Assert(maps_[ig0]->iCof.index != -1, fmt::format("{}->iCof.index != -1", maps_[ig0]->name));
                }
                if (json_comm.at("g").size() >= 2) {
                    auto g1 = json_comm.at("g")[1].get<std::string>();
                    ig1 = GetMapIndexByName(g1);
                    ErrorIdMsg::Assert(ig1 >= 0, "ig1 >= 0");
                    ErrorIdMsg::Assert(maps_[ig1]->iCof.index != -1, fmt::format("{}->iCof.index != -1", maps_[ig1]->name));
                }
                auto name = json_comm.at("name").get<std::string>();
                comms_.push_back(Commutativity{name, if0, if1, ig0, ig1});
            }
        }

        if (loadD2) {
            for (size_t iCw = 0; iCw < rings_.size(); ++iCw) {
                if (basis_d2[iCw].empty())
                    continue;
                auto& name = rings_[iCw].name;
                auto& nodes_ss = rings_[iCw].nodes_ss;
                for (auto& [deg, d2] : basis_d2[iCw]) {
                    for (size_t i = 0; i < d2.size(); ++i) {
                        if (IsNewDiff(nodes_ss, deg, int1d{(int)i}, d2[i], 2)) {
                            Logger::LogDiff(0, EnumReason::d2, name, deg, 2, int1d{(int)i}, d2[i], "", flags);
                            if (auto rt = SetRingDiffGlobal(iCw, deg, int1d{(int)i}, d2[i], 2, true, flags))
                                throw ErrorIdMsg(rt.code, "Failed to SetRingDiffGlobal()");
                        }
                    }
                }
            }
            for (size_t iMod = 0; iMod < modules_.size(); ++iMod) {
                auto iCw = rings_.size() + iMod;
                if (basis_d2[iCw].empty())
                    continue;
                auto& name = modules_[iMod].name;
                auto& nodes_ss = modules_[iMod].nodes_ss;
                for (auto& [deg, d2] : basis_d2[iCw]) {
                    for (size_t i = 0; i < d2.size(); ++i) {
                        if (IsNewDiff(nodes_ss, deg, int1d{(int)i}, d2[i], 2)) {
                            Logger::LogDiff(0, EnumReason::d2, name, deg, 2, int1d{(int)i}, d2[i], "", flags);
                            if (auto rt = SetModuleDiffGlobal(iMod, deg, int1d{(int)i}, d2[i], 2, true, flags))
                                throw ErrorIdMsg(rt.code, "Failed to SetModuleDiffGlobal()");
                        }
                    }
                }
            }
        }
    }
    /*catch (nlohmann::detail::exception& e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }*/
    catch (NoException&) {
    }

    /*if (flags & SSFlag::homotopy)
        UpdateAllPossEinf();*/

    // VersionConvertReorderRels();
}

void Category::SetDeduceList(const std::vector<std::string>& cws)
{
    deduce_list_spectra_.clear();
    for (auto& name : cws) {
        auto iCw = GetIndexCwByName(name);
        deduce_list_spectra_.push_back(iCw);
    }
}

void Category::SetCofDeduceList(const std::vector<std::string>& cofs)
{
    deduce_list_cofseq_.clear();
    for (auto& name : cofs) {
        auto iCw = GetIndexCofByName(name);
        deduce_list_cofseq_.push_back(iCw.index);
    }
}

void Category::LoadNodes(const std::string& cat_name, const std::string& ckpt_name, SSFlag flags)
{
    std::string cat_dir;
    if (ckpt_name.size())
        cat_dir = fmt::format("{}/checkpoints/{}", cat_name, ckpt_name);
    else
        cat_dir = cat_name;
    MyException::Assert(depth_ == 0, "Need depth==0");

    /* Load rings */
    auto& json_rings = js_["rings"];
    for (auto& json_ring : json_rings) {
        std::string name = json_ring["name"].get<std::string>(), path = json_ring["path"].get<std::string>();
        std::string abs_path = fmt::format("{}/{}", cat_dir, path);
        std::string table_prefix = fmt::format("{}_AdamsE2", name);

        DBSS db(abs_path);
        size_t iRing = GetIndexCwByName(name).index;
        rings_[iRing].nodes_ss[0] = db.load_ss(table_prefix);
        rings_[iRing].nodes_ss[1] = {};
    }

    /* Load modules */
    auto& json_mods = js_["modules"];
    for (auto& json_mod : json_mods) {
        std::string name = json_mod["name"].get<std::string>(), path = json_mod["path"].get<std::string>();
        std::string abs_path = fmt::format("{}/{}", cat_dir, path);
        std::string table_prefix = fmt::format("{}_AdamsE2", name);

        DBSS db(abs_path);
        size_t iMod = GetIndexCwByName(name).index;
        modules_[iMod].nodes_ss[0] = db.load_ss(table_prefix);
        modules_[iMod].nodes_ss[1] = {};
    }

    /* Load cofseq */
    if (flags & SSFlag::cofseq) {
        auto& json_cofseqs = js_.at("cofseqs");
        for (size_t i = 0; i < cofseqs_.size(); ++i) {
            auto path_cofseq = json_cofseqs[i].at("path").get<std::string>();
            auto abs_path_cofseq = fmt::format("{}/{}", cat_dir, path_cofseq);
            auto name = json_cofseqs[i].at("name").get<std::string>();

            DBSS db(abs_path_cofseq);
            auto node_cofseq = db.load_cofseq(fmt::format("cofseq_{}", cofseqs_[i].name));
            for (size_t iTri = 0; iTri < 3; ++iTri) {
                cofseqs_[i].nodes_cofseq[iTri].front() = std::move(node_cofseq[iTri]);
            }
        }
    }
}

void Category::SaveNodes(const std::string& cat_root, const std::string& ckpt_name, bool ss_incremental, SSFlag flags)
{
    if (flags & SSFlag::no_save)
        return;
    std::string cat_dir;
    if (ckpt_name.size()) {
        std::filesystem::create_directory(fmt::format("{}/checkpoints", cat_root));
        cat_dir = fmt::format("{}/checkpoints/{}", cat_root, ckpt_name);
        std::filesystem::create_directory(cat_dir);
    }
    else
        cat_dir = cat_root;

    try {
        /* #save rings */
        auto& json_rings = js_["rings"];
        for (auto& json_ring : json_rings) {
            std::string name = json_ring["name"].get<std::string>(), path = json_ring["path"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", cat_dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            db.begin_transaction();
            size_t iRing = GetIndexCwByName(name).index;
            auto& ring = rings_[iRing];
            if (ss_incremental)
                db.update_ss(table_prefix, ring.nodes_ss);
            else
                db.save_ss(table_prefix, ring.nodes_ss);

            if (flags & SSFlag::pi) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);

                db.drop_and_create_pi_generators(name);
                db.save_pi_generators(name, ring.pi_gb.gen_degs(), ring.pi_gen_Einf);
                db.save_pi_gb(name, ring.pi_gb.OutputForDatabase(), GetRingGbEinf(iRing));
                db.save_pi_basis(name, ring.nodes_pi_basis.front());
                if (flags & SSFlag::pi_def) {
                    db.drop_and_create_pi_definitions(name);
                    db.save_pi_def(name, ring.pi_gen_defs, ring.pi_gen_def_mons);
                }
            }
            db.end_transaction();
        }

        /* #save modules */
        auto& json_mods = js_["modules"];
        for (auto& json_mod : json_mods) {
            std::string name = json_mod["name"].get<std::string>(), path = json_mod["path"].get<std::string>();
            std::string over = json_mod["over"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", cat_dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            db.begin_transaction();
            size_t iMod = GetIndexCwByName(name).index;
            auto& mod = modules_[iMod];
            if (ss_incremental)
                db.update_ss(table_prefix, mod.nodes_ss);
            else
                db.save_ss(table_prefix, mod.nodes_ss);

            if (flags & SSFlag::pi) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);
                if (flags & SSFlag::pi_def)
                    db.drop_and_create_pi_generators_mod(name);
                db.save_pi_generators_mod(name, mod.pi_gb.v_degs(), mod.pi_gen_Einf);
                db.save_pi_gb_mod(name, mod.pi_gb.OutputForDatabase(), GetModuleGbEinf(iMod));
                db.save_pi_basis_mod(name, mod.nodes_pi_basis.front());
                if (flags & SSFlag::pi_def) {
                    db.drop_and_create_pi_definitions(name);
                    db.save_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
                }
            }
            db.end_transaction();
        }

        /* #save cofseq */
        if (flags & SSFlag::cofseq) {
            auto& json_cofseqs = js_.at("cofseqs");
            for (size_t i = 0; i < cofseqs_.size(); ++i) {
                auto path_cofseq = json_cofseqs[i].at("path").get<std::string>();
                auto abs_path_cofseq = fmt::format("{}/{}", cat_dir, path_cofseq);
                auto name = json_cofseqs[i].at("name").get<std::string>();
                DBSS db(abs_path_cofseq);
                db.begin_transaction();
                auto table = fmt::format("cofseq_{}", name);
                db.create_cofseq(table);
                db.save_cofseq(table, cofseqs_[i], ss_incremental);
                db.end_transaction();
            }
        }
    }
    catch (nlohmann::detail::exception& e) {
        fmt::print("JsonError({}): {}\n", e.id, e.what());
        throw e;
    }
}

void Category::PrintSummary() const
{
    fmt::print("Summary of changes:     \n");

    /* ss */
    int count = 0;
    for (auto& ring : rings_)
        for (auto _ : ring.nodes_ss.changes().degs())
            ++count;
    for (auto& mod : modules_)
        for (auto _ : mod.nodes_ss.changes().degs())
            ++count;
    fmt::print("  ss: {}\n", count);

    /* cofseq */
    count = 0;
    for (auto& cofseq : cofseqs_)
        for (size_t iTri = 0; iTri < 3; ++iTri)
            for (auto _ : cofseq.nodes_cofseq[iTri].changes().degs())
                ++count;
    fmt::print("  cofseq: {}\n", count);
    fmt::print("\n");
}