/* @file config.cpp
 * Load settings from json
 */

#include "json.h"
#include "main.h"
#include "mylog.h"
#include <filesystem>
#include <fstream>
#include <regex>

void LoadJson(const std::string& diagram_name, nlohmann::json& root_json, nlohmann::json& diag_json)
{
    using json = nlohmann::json;

    root_json = myio::load_json("ss.json");
    try {
        json& diagrams = root_json.at("diagrams");
        std::string dir = diagrams.contains(diagram_name) ? diagrams[diagram_name].get<std::string>() : diagram_name;
        diag_json = myio::load_json(fmt::format("{}/ss.json", dir));
    }
    catch (nlohmann::detail::exception& e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }
}

void GetAllDbNames(const std::string& diagram_name, std::vector<std::string>& names, std::vector<std::string>& paths, std::vector<int>& isRing, bool log)
{
    using json = nlohmann::json;
    json js = myio::load_json("ss.json");

    try {
        json& diagrams = js.at("diagrams");
        std::string dir = diagrams.contains(diagram_name) ? diagrams[diagram_name].get<std::string>() : diagram_name;
        json diagram = myio::load_json(fmt::format("{}/ss.json", dir));
        if (log) {
            auto path_log = fmt::format("{}/{}", dir, diagram.at("log").get<std::string>());
            Logger::SetOutDeduce(path_log.c_str());
        }

        /* #Load rings */
        json& json_rings = diagram.at("rings");
        for (auto& json_ring : json_rings) {
            std::string name = json_ring.at("name").get<std::string>(), path = json_ring.at("path").get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            names.push_back(name);
            paths.push_back(abs_path);
            isRing.push_back(1);
        }

        /* #Load modules */
        json& json_mods = diagram.at("modules");
        for (auto& json_mod : json_mods) {
            std::string name = json_mod.at("name").get<std::string>(), path = json_mod.at("path").get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            names.push_back(name);
            paths.push_back(abs_path);
            isRing.push_back(0);
        }
    }
    catch (nlohmann::detail::exception& e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }
}

Diagram::Diagram(std::string diagram_name, DeduceFlag flag, bool log)
{
    using json = nlohmann::json;
    json js = myio::load_json("ss.json");

    try {
        json& diagrams = js.at("diagrams");
        std::string dir = diagrams.contains(diagram_name) ? diagrams[diagram_name].get<std::string>() : diagram_name;
        json diagram = myio::load_json(fmt::format("{}/ss.json", dir));
        if (log) {
            auto path_log = fmt::format("{}/{}", dir, diagram.at("log").get<std::string>());
            Logger::SetOutDeduce(path_log.c_str());
        }

        /*# Load rings */

        json& json_rings = diagram.at("rings");
        rings_.reserve(json_rings.size());
        size_t iCw = 0;

        for (auto& json_ring : json_rings) {
            std::string name = json_ring.at("name").get<std::string>(), path = json_ring.at("path").get<std::string>();
            if (!json_ring.contains("deduce") || json_ring.at("deduce").get<std::string>() == "on")
                deduce_list_spectra_.push_back(iCw);
            std::string abs_path = fmt::format("{}/{}", dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            myio::AssertFileExists(abs_path);
            DBSS db(abs_path);
            RingSp ring;
            ring.name = name;
            ring.basis = db.load_basis(table_prefix);
            ring.degs_basis_order_by_stem = OrderDegsByStem(ring.basis);
            ring.t_max = ring.basis.rbegin()->first.t;
            ring.nodes_ss = {db.load_ss(table_prefix), {}};
            ring.nodes_ss.reserve(MAX_DEPTH + 3);
            ring.gb = Groebner(ring.t_max, {}, db.load_gb(table_prefix, DEG_MAX));

            if (flag & DeduceFlag::pi) {
                ring.pi_gen_Einf = db.get_column_from_str<Poly>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Poly>);
                ring.pi_gb = algZ::Groebner(ring.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb(name, DEG_MAX), true);
                ring.nodes_pi_basis.reserve(MAX_DEPTH + 2);
                if (flag & DeduceFlag::pi_def)
                    db.load_pi_def(name, ring.pi_gen_defs, ring.pi_gen_def_mons);
            }

            rings_.push_back(std::move(ring));
            ++iCw;
        }

        /*# Load modules */

        json& json_mods = diagram.at("modules");
        modules_.reserve(json_mods.size());
        for (auto& json_mod : json_mods) {
            std::string name = json_mod.at("name").get<std::string>(), path = json_mod.at("path").get<std::string>();
            std::string over = json_mod.at("over").get<std::string>();
            if (!json_mod.contains("deduce") || json_mod.at("deduce").get<std::string>() == "on")
                deduce_list_spectra_.push_back(iCw);
            std::string abs_path = fmt::format("{}/{}", dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            myio::AssertFileExists(abs_path);
            DBSS db(abs_path);
            ModSp mod;
            mod.name = name;
            mod.iRing = (size_t)GetRingIndexByName(over);
            MyException::Assert(mod.iRing != -1, "mod.iRing != -1");
            auto& ring = rings_[mod.iRing];
            mod.basis = db.load_basis_mod(table_prefix);
            mod.degs_basis_order_by_stem = OrderDegsByStem(mod.basis);
            mod.t_max = mod.basis.rbegin()->first.t;
            mod.nodes_ss = {db.load_ss(table_prefix), {}};
            mod.nodes_ss.reserve(MAX_DEPTH + 1);
            Mod1d xs = db.load_gb_mod(table_prefix, DEG_MAX);
            mod.gb = GroebnerMod(&ring.gb, mod.t_max, {}, std::move(xs));

            if (flag & DeduceFlag::pi) {
                mod.pi_gen_Einf = db.get_column_from_str<Mod>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Mod>);
                mod.pi_gb = algZ::GroebnerMod(&ring.pi_gb, mod.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb_mod(name, DEG_MAX), true);
                mod.nodes_pi_basis.reserve(MAX_DEPTH);
                if (flag & DeduceFlag::pi_def)
                    db.load_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
            }

            modules_.push_back(std::move(mod));
            ring.ind_mods.push_back(modules_.size() - 1);
            ++iCw;
        }

        /*# Load maps */
        std::regex is_map_regex("^map_AdamsSS_(\\w+?_to_\\w+?)(?:_t\\d+|).db$"); /* match example: map_AdamsSS_RP1_4_to_RP3_4_t169.db */
        std::smatch match;

        json json_maps = diagram.at("maps");
        if (flag & DeduceFlag::cofseq) {
            for (auto& item : diagram.at("maps_v2"))
                json_maps.push_back(item);
        }
        maps_.reserve(json_maps.size());
        for (auto& json_map : json_maps) {
            auto name = json_map.at("name").get<std::string>();
            auto display = json_map.contains("display") ? json_map["display"].get<std::string>() : name;
            auto t_max = json_map["t_max"].get<int>();
            std::string from = json_map.at("from").get<std::string>(), to = json_map.at("to").get<std::string>();
            std::unique_ptr<Map> map;

            if (json_map.contains("factor")) {
                auto& strt_factor = json_map.at("factor");
                int stem = strt_factor[0].get<int>(), s = strt_factor[1].get<int>(), i_factor = strt_factor[2].get<int>();
                AdamsDeg deg_factor(s, stem + s);
                if (int index_from = GetRingIndexByName(from); index_from != -1) {
                    if (int index_to = GetRingIndexByName(to); index_to != -1) {
                        MyException::Assert(index_from == index_to, "index_from == index_to");
                        Poly factor = rings_[index_from].basis.at(deg_factor)[i_factor];
                        map = std::make_unique<MapMulRing2Ring>(name, display, t_max, deg_factor, (size_t)index_from, std::move(factor));
                    }
                    else {
                        index_to = GetModuleIndexByName(to);
                        Mod factor = modules_[index_to].basis.at(deg_factor)[i_factor];
                        map = std::make_unique<MapMulRing2Mod>(name, display, t_max, deg_factor, (size_t)index_from, (size_t)index_to, std::move(factor));
                    }
                }
                else {
                    index_from = GetModuleIndexByName(from);
                    int index_to = GetModuleIndexByName(to);
                    MyException::Assert(index_from == index_to, "index_from == index_to");
                    Poly factor = rings_[modules_[index_from].iRing].basis.at(deg_factor)[i_factor];
                    map = std::make_unique<MapMulMod2Mod>(name, display, t_max, deg_factor, (size_t)index_from, std::move(factor));
                }
            }
            else {
                int fil = json_map.contains("fil") ? json_map["fil"].get<int>() : 0;
                int sus = json_map.contains("sus") ? json_map["sus"].get<int>() : 0;
                AdamsDeg deg(fil, fil - sus);

                std::string path = json_map.at("path").get<std::string>();
                std::string abs_path = fmt::format("{}/{}", dir, path);
                myio::AssertFileExists(abs_path);
                DBSS db(abs_path);
                std::string table;
                if (std::regex_search(path, match, is_map_regex); match[0].matched) {
                    table = fmt::format("map_AdamsE2_{}", match[1].str());
                }
                else {
                    fmt::print("filename={} not supported.\n", path);
                    throw MyException(0x839393b2, "File name is not supported.");
                }

                if (GetRingIndexByName(from) != -1) {
                    size_t index_from = (size_t)GetRingIndexByName(from);
                    size_t index_to = (size_t)GetRingIndexByName(to);
                    MyException::Assert(index_from != -1 && index_to != -1, "index_from != -1 && index_to != -1");
                    auto images = db.get_column_from_str<Poly>(table, "map", "ORDER BY id", myio::Deserialize<Poly>);
                    map = std::make_unique<MapRing2Ring>(name, display, t_max, deg, index_from, index_to, std::move(images));
                    rings_[index_from].ind_maps.push_back(maps_.size());
                }
                else {
                    size_t index_from = (size_t)GetModuleIndexByName(from);
                    MyException::Assert(index_from != -1, "index_from != -1");
                    if (GetRingIndexByName(to) != -1) {
                        size_t index_to = (size_t)GetRingIndexByName(to);
                        MyException::Assert(index_to != -1, "index_to != -1");
                        auto images = db.get_column_from_str<Poly>(table, "map", "", myio::Deserialize<Poly>);
                        if (!json_map.contains("over")) {
                            map = std::make_unique<MapMod2Ring>(name, display, t_max, deg, index_from, index_to, std::move(images));
                        }
                        else {
                            std::string over = json_map.at("over").get<std::string>();
                            int index_map = ut::IndexOf(maps_, [&over](const std::unique_ptr<Map>& map) { return map->name == over; });
                            map = std::make_unique<MapMod2RingV2>(name, display, t_max, deg, index_from, index_to, index_map, std::move(images));
                        }
                    }
                    else {
                        size_t index_to = (size_t)GetModuleIndexByName(to);
                        MyException::Assert(index_to != -1, "index_to != -1");
                        auto images = db.get_column_from_str<Mod>(table, "map", "", myio::Deserialize<Mod>);
                        if (!json_map.contains("over")) {
                            map = std::make_unique<MapMod2Mod>(name, display, t_max, deg, index_from, index_to, std::move(images));
                        }
                        else {
                            std::string over = json_map.at("over").get<std::string>();
                            int index_map = ut::IndexOf(maps_, [&over](const std::unique_ptr<Map>& map) { return map->name == over; });
                            map = std::make_unique<MapMod2ModV2>(name, display, t_max, deg, index_from, index_to, index_map, std::move(images));
                        }
                    }
                    modules_[index_from].ind_maps.push_back(maps_.size());
                }
            }
            maps_.push_back(std::move(map));
        }

        if (flag & DeduceFlag::cofseq) {
            auto& json_cofseqs = diagram.at("cofseqs");
            std::vector<std::string> cofseq_maps;
            for (auto& json_cofseq : json_cofseqs) {
                CofSeq cofseq;
                cofseq.name = json_cofseq.at("name").get<std::string>();
                auto path_cofseq = json_cofseq.at("path").get<std::string>();
                auto abs_path_cofseq = fmt::format("{}/{}", dir, path_cofseq);
                bool fileExists = myio::FileExists(abs_path_cofseq);
                DBSS db(abs_path_cofseq);
                db.create_cofseq(fmt::format("cofseq_{}", cofseq.name));
                std::array<std::string, 3> maps_names = {json_cofseq.at("i").get<std::string>(), json_cofseq.at("q").get<std::string>(), json_cofseq.at("d").get<std::string>()};
                for (size_t iMap = 0; iMap < maps_names.size(); ++iMap) {
                    size_t index = (size_t)GetMapIndexByName(maps_names[iMap]);
                    const auto& map = maps_[index];
                    cofseq.indexMap[iMap] = index;
                    cofseq.degMap[iMap] = map->deg;
                    cofseq.isRing[iMap] = map->IsFromRing(cofseq.indexCw[iMap]);
                    cofseq.nameCw[iMap] = cofseq.isRing[iMap] ? rings_[cofseq.indexCw[iMap]].name : modules_[cofseq.indexCw[iMap]].name;
                    cofseq.nodes_ss[iMap] = cofseq.isRing[iMap] ? &rings_[cofseq.indexCw[iMap]].nodes_ss : &modules_[cofseq.indexCw[iMap]].nodes_ss;

                    cofseq.nodes_cofseq[iMap].resize(cofseq.nodes_ss[iMap]->size());
                    if (!fileExists) { /* Initialize cofseq */
                        auto& front_cofseq = cofseq.nodes_cofseq[iMap].front();
                        for (auto& [deg, sc] : cofseq.nodes_ss[iMap]->front()) {
                            size_t first_PC = GetFirstIndexOnLevel(sc, LEVEL_PERM);
                            size_t last_PC = GetFirstIndexOnLevel(sc, LEVEL_PERM + 1);
                            for (size_t i = first_PC; i < last_PC; ++i) {
                                front_cofseq[deg].basis.push_back(sc.basis[i]);
                                front_cofseq[deg].diffs.push_back(NULL_DIFF);
                                front_cofseq[deg].levels.push_back(LEVEL_MAX);
                            }
                        }
                    }
                }
                if (fileExists) { /* Load cofseq from database */
                    auto& node_cofseq = db.load_cofseq(fmt::format("cofseq_{}", cofseq.name));
                    for (size_t i = 0; i < cofseq.nodes_cofseq.size(); ++i)
                        cofseq.nodes_cofseq[i].front() = std::move(node_cofseq[i]);
                }
                else { /* Add differentials in cofseq */

                }
                cofseqs_.push_back(std::move(cofseq));
            }
        }
    }
    /*catch (nlohmann::detail::exception& e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }*/
    catch (NoException&) {
    }

    /*if (flag & DeduceFlag::homotopy)
        UpdateAllPossEinf();*/

    // VersionConvertReorderRels();
}

void Diagram::save(std::string diagram_name, DeduceFlag flag)
{
    using json = nlohmann::json;
    json js = myio::load_json("ss.json");

    std::map<std::string, std::string> paths;
    try {
        json& diagrams = js.at("diagrams");
        std::string dir = diagrams.contains(diagram_name) ? diagrams[diagram_name].get<std::string>() : diagram_name;
        json diagram = myio::load_json(fmt::format("{}/ss.json", dir));

        /* #save rings */
        json& json_rings = diagram["rings"];
        for (auto& json_ring : json_rings) {
            std::string name = json_ring["name"].get<std::string>(), path = json_ring["path"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            db.begin_transaction();
            size_t iRing = (size_t)GetRingIndexByName(name);
            auto& ring = rings_[iRing];
            db.update_ss(table_prefix, ring.nodes_ss[1]);

            if (flag & DeduceFlag::pi) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);

                db.drop_and_create_pi_generators(name);
                db.save_pi_generators(name, ring.pi_gb.gen_degs(), ring.pi_gen_Einf);
                db.save_pi_gb(name, ring.pi_gb.OutputForDatabase(), GetRingGbEinf(iRing));
                db.save_pi_basis(name, ring.nodes_pi_basis.front());
                if (flag & DeduceFlag::pi_def) {
                    db.drop_and_create_pi_definitions(name);
                    db.save_pi_def(name, ring.pi_gen_defs, ring.pi_gen_def_mons);
                }
            }
            db.end_transaction();
        }

        /* #save modules */
        json& json_mods = diagram["modules"];
        for (auto& json_mod : json_mods) {
            std::string name = json_mod["name"].get<std::string>(), path = json_mod["path"].get<std::string>();
            std::string over = json_mod["over"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            db.begin_transaction();
            size_t iMod = (size_t)GetModuleIndexByName(name);
            auto& mod = modules_[iMod];
            db.update_ss(table_prefix, mod.nodes_ss[1]);

            if (flag & DeduceFlag::pi) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);
                if (flag & DeduceFlag::pi_def)
                    db.drop_and_create_pi_generators_mod(name);
                db.save_pi_generators_mod(name, mod.pi_gb.v_degs(), mod.pi_gen_Einf);
                db.save_pi_gb_mod(name, mod.pi_gb.OutputForDatabase(), GetModuleGbEinf(iMod));
                db.save_pi_basis_mod(name, mod.nodes_pi_basis.front());
                if (flag & DeduceFlag::pi_def) {
                    db.drop_and_create_pi_definitions(name);
                    db.save_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
                }
            }
            db.end_transaction();
        }

        /* #save cofseq */
        if (flag & DeduceFlag::cofseq) {
            auto& json_cofseqs = diagram.at("cofseqs");
            std::vector<std::string> cofseq_maps;
            for (size_t i = 0; i < cofseqs_.size(); ++i) {
                auto path_cofseq = json_cofseqs[i].at("path").get<std::string>();
                auto abs_path_cofseq = fmt::format("{}/{}", dir, path_cofseq);
                auto name = json_cofseqs[i].at("name").get<std::string>();
                DBSS db(abs_path_cofseq);
                auto table = fmt::format("cofseq_{}", name);
                db.create_cofseq(table);
                db.save_cofseq(table, cofseqs_[i]);
            }
        }
    }
    catch (nlohmann::detail::exception& e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }
}
