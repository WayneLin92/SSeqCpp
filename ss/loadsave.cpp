/* @file config.cpp
 * Load settings from json
 */

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

        auto path_log = fmt::format("{}/{}", dir, diag_json.at("log").get<std::string>());
        Logger::SetOutDeduce(path_log.c_str());
    }
    catch (nlohmann::detail::exception& e) {
        fmt::print("JsonError({}): {}\n", e.id, e.what());
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
            MyException::Assert(json_f.at("from").get<std::string>() == from, "Valid composition");
        if (i + 1 == indices.size())
            MyException::Assert(json_f.at("to").get<std::string>() == to, "Valid composition");
        else
            MyException::Assert(json_f.at("to").get<std::string>() == json_maps[indices[i + 1]].at("from").get<std::string>(), "Valid composition");
    }
}

int1d get_compostion(const int1d& x, AdamsDeg deg, const Diagram& diagram, const PMap1d& maps, const int1d& indices)
{
    int1d result = x;
    for (size_t i = 0; i < indices.size(); ++i) {
        result = maps[indices[i]]->map(result, deg, diagram);
        deg += maps[indices[i]]->deg;
    }
    return result;
}

Diagram::Diagram(std::string diagram_name, SSFlag flag, bool log, bool loadD2)
{
    try {
        js_ = myio::load_json(fmt::format("{}/ss.json", diagram_name));
        if (log) {
            auto path_log = fmt::format("{}/{}", diagram_name, js_.at("log").get<std::string>());
            Logger::SetOutDeduce(path_log.c_str());
        }

        /*# Load rings and modules */

        auto& json_rings = js_.at("rings");
        auto& json_mods = js_.at("modules");
        AdamsDeg2d ring_gen_degs;
        AdamsDeg2d module_gen_degs;
        rings_.reserve(json_rings.size());
        modules_.reserve(json_mods.size());
        std::vector<std::map<AdamsDeg, int2d>> basis_d2;
        {
            size_t iCw = 0;
            for (auto& json_ring : json_rings) {
                std::string name = json_ring.at("name").get<std::string>(), path = json_ring.at("path").get<std::string>();
                if (myio::get(json_ring, "deduce", "on") == "on")
                    deduce_list_spectra_.push_back(IndexRing(iCw));
                std::string abs_path = fmt::format("{}/{}", diagram_name, path);
                std::string table_prefix = fmt::format("{}_AdamsE2", name);

                myio::AssertFileExists(abs_path);
                DBSS db(abs_path);
                RingSp ring;
                ring.name = name;
                ring.basis = db.load_basis(table_prefix);
                if (loadD2)
                    basis_d2.push_back(db.load_basis_d2(table_prefix));
                ring.degs_basis_order_by_stem = OrderDegsByStem(ring.basis);
                ring.t_max = ring.basis.rbegin()->first.t;
                ring.nodes_ss = {db.load_ss(table_prefix), {}};
                ring.nodes_ss.reserve(MAX_DEPTH + 2);
                ring.gb = Groebner(ring.t_max, {}, db.load_gb(table_prefix, DEG_MAX));

                if (flag & SSFlag::pi) {
                    ring.pi_gen_Einf = db.get_column_from_str<Poly>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Poly>);
                    ring.pi_gb = algZ::Groebner(ring.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb(name, DEG_MAX));
                    ring.nodes_pi_basis.reserve(MAX_DEPTH + 1);
                    if (flag & SSFlag::pi_def)
                        db.load_pi_def(name, ring.pi_gen_defs, ring.pi_gen_def_mons);
                }

                ring_gen_degs.push_back(db.load_gen_adamsdegs(table_prefix));
                rings_.push_back(std::move(ring));
                ++iCw;
            }

            iCw = 0;
            for (auto& json_mod : json_mods) {
                std::string name = json_mod.at("name").get<std::string>(), path = json_mod.at("path").get<std::string>();
                std::string over = json_mod.at("over").get<std::string>();
                if (myio::get(json_mod, "deduce", "on") == "on")
                    deduce_list_spectra_.push_back(IndexMod(iCw));
                std::string abs_path = fmt::format("{}/{}", diagram_name, path);
                std::string table_prefix = fmt::format("{}_AdamsE2", name);

                myio::AssertFileExists(abs_path);
                DBSS db(abs_path);
                ModSp mod;
                mod.name = name;
                auto indexRing = GetIndexCwByName(over);
                MyException::Assert(indexRing.isRing, "indexRing.isRing");
                mod.iRing = indexRing.index;
                auto& ring = rings_[mod.iRing];
                mod.basis = db.load_basis_mod(table_prefix);
                if (loadD2)
                    basis_d2.push_back(db.load_basis_d2(table_prefix));
                mod.degs_basis_order_by_stem = OrderDegsByStem(mod.basis);
                mod.t_max = mod.basis.rbegin()->first.t;
                mod.nodes_ss = {db.load_ss(table_prefix), {}};
                mod.nodes_ss.reserve(MAX_DEPTH + 2);
                Mod1d xs = db.load_gb_mod(table_prefix, DEG_MAX);
                mod.gb = GroebnerMod(&ring.gb, mod.t_max, {}, std::move(xs));

                if (flag & SSFlag::pi) {
                    mod.pi_gen_Einf = db.get_column_from_str<Mod>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Mod>);
                    mod.pi_gb = algZ::GroebnerMod(&ring.pi_gb, mod.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb_mod(name, DEG_MAX));
                    mod.nodes_pi_basis.reserve(MAX_DEPTH + 1);
                    if (flag & SSFlag::pi_def)
                        db.load_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
                }

                modules_.push_back(std::move(mod));
                module_gen_degs.push_back(db.load_gen_adamsdegs(table_prefix));
                ring.ind_mods.push_back(modules_.size() - 1);
                ++iCw;
            }
        }

        /*# Load maps */
        std::regex is_map_regex("^map_AdamsSS_(\\w+?_(?:to|)_\\w+?)(?:_t\\d+|).db$"); /* match example: map_AdamsSS_RP1_4_to_RP3_4_t169.db */
        std::smatch match;

        auto& json_maps = js_.at("maps");
        if (flag & SSFlag::cofseq || flag & SSFlag::naming) {
            for (auto& item : js_.at("maps_v2"))
                json_maps.push_back(item);
        }
        maps_.reserve(json_maps.size());
        for (auto& json_map : json_maps) {
            auto name = json_map.at("name").get<std::string>();
            auto display = json_map.contains("display") ? json_map["display"].get<std::string>() : name;
            std::string from = json_map.at("from").get<std::string>(), to = json_map.at("to").get<std::string>();
            std::unique_ptr<Map> map;

            if (json_map.contains("factor")) {
                auto& strt_factor = json_map.at("factor");
                int stem = strt_factor[0].get<int>(), s = strt_factor[1].get<int>();
                int i_factor = strt_factor.size() > 2 ? strt_factor[2].get<int>() : -1;
                AdamsDeg deg_factor(s, stem + s);
                if (auto ifrom = GetIndexCwByName(from); ifrom.isRing) {
                    if (auto ito = GetIndexCwByName(to); ito.isRing) {
                        MyException::Assert(ifrom == ito, "index_from == index_to");
                        Poly factor = i_factor != -1 ? rings_[ifrom.index].basis.at(deg_factor)[i_factor] : Poly{};
                        int t_max = rings_[ito.index].t_max - deg_factor.t;
                        map = std::make_unique<MapMulRing2Ring>(name, display, t_max, deg_factor, ifrom.index, std::move(factor));
                    }
                    else {
                        Mod factor = i_factor != -1 ? modules_[ito.index].basis.at(deg_factor)[i_factor] : Mod{};
                        int t_max = modules_[ito.index].t_max - deg_factor.t;
                        map = std::make_unique<MapMulRing2Mod>(name, display, t_max, deg_factor, ifrom.index, ito.index, std::move(factor));
                    }
                }
                else {
                    auto ito = GetIndexCwByName(to);
                    MyException::Assert(ifrom == ito, "index_from == index_to");
                    Poly factor = rings_[modules_[ifrom.index].iRing].basis.at(deg_factor)[i_factor];
                    int t_max = modules_[ito.index].t_max - deg_factor.t;
                    map = std::make_unique<MapMulMod2Mod>(name, display, t_max, deg_factor, ifrom.index, std::move(factor));
                }
            }
            else if (json_map.contains("composition")) {
                auto& json_composition = json_map.at("composition");
                int1d indices;
                for (auto& json_f : json_composition) {
                    int index_map = GetMapIndexByName(json_f.get<std::string>());
                    MyException::Assert(index_map != -1, "index_map != -1");
                    indices.push_back(index_map);
                }
                int t_max_map;
                AdamsDeg deg_map;
                if (auto ifrom = GetIndexCwByName(from); ifrom.isRing) {
                    auto ito = GetIndexCwByName(to);
                    MyException::Assert(ito.isRing, fmt::format("{}: ito.isRing", to));

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

                    if (auto ito = GetIndexCwByName(to); ito.isRing) {
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
                std::string abs_path = fmt::format("{}/{}", diagram_name, path);
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

                if (auto ifrom = GetIndexCwByName(from); ifrom.isRing) {
                    auto ito = GetIndexCwByName(to);
                    MyException::Assert(ito.isRing, fmt::format("{}: ito.isRing", to));
                    auto images = db.get_column_from_str<Poly>(table, "map", "ORDER BY id", myio::Deserialize<Poly>);
                    map = std::make_unique<MapRing2Ring>(name, display, t_max, deg, ifrom.index, ito.index, std::move(images));
                    rings_[ifrom.index].ind_maps.push_back(maps_.size());
                    rings_[ito.index].ind_maps_prev.push_back(maps_.size());
                }
                else {
                    if (auto ito = GetIndexCwByName(to); ito.isRing) {
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

        /*# Load cofseqs */
        if (flag & SSFlag::cofseq) {
            auto& json_cofseqs = js_.at("cofseqs");
            std::vector<std::string> cofseq_maps;
            size_t iCof = 0;
            for (auto& json_cofseq : json_cofseqs) {
                CofSeq cofseq;
                cofseq.name = json_cofseq.at("name").get<std::string>();
                auto path_cofseq = json_cofseq.at("path").get<std::string>();
                auto abs_path_cofseq = fmt::format("{}/{}", diagram_name, path_cofseq);
                std::array<std::string, 3> maps_names = {json_cofseq.at("i").get<std::string>(), json_cofseq.at("q").get<std::string>(), json_cofseq.at("d").get<std::string>()};
                if (!json_cofseq.contains("deduce") || json_cofseq.at("deduce").get<std::string>() == "on")
                    deduce_list_cofseq_.push_back(iCof);
                for (size_t iCs = 0; iCs < 3; ++iCs) {
                    size_t index = (size_t)GetMapIndexByName(maps_names[iCs]);
                    MyException::Assert(index != -1, maps_names[iCs] + " not found");
                    const auto& map = maps_[index];
                    cofseq.indexMap[iCs] = index;
                    map->ind_cof = IndexCof{iCof, iCs};
                    cofseq.degMap[iCs] = map->deg;
                    cofseq.indexCw[iCs] = map->from;
                    GetIndexCof(cofseq.indexCw[iCs]).push_back(IndexCof{iCof, iCs});
                    cofseq.nameCw[iCs] = GetCwName(cofseq.indexCw[iCs]);
                    cofseq.t_max[iCs] = GetTMax(cofseq.indexCw[iCs]);
                    cofseq.nodes_ss[iCs] = &GetSS(cofseq.indexCw[iCs]);
                    cofseq.nodes_cofseq[iCs].resize(cofseq.nodes_ss[iCs]->size());
                    cofseq.nodes_cofseq[iCs].reserve(cofseq.nodes_ss[iCs]->capacity());
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
            SyncCofseq(flag);

            if (flag & SSFlag::depth_ss_cofseq)
                ut::RemoveIf(deduce_list_spectra_, [&](IndexCw iCw) {
                    const auto& ind_cofs = GetIndexCof(iCw);
                    return ind_cofs.empty();
                });
        }

        /*# Load commutativities */
        if (flag & SSFlag::cofseq && js_.contains("commutativity")) {
            auto& json_comms = js_.at("commutativity");
            for (auto& json_comm : json_comms) {
                auto f0 = json_comm.at("f")[0].get<std::string>();
                auto f1 = json_comm.at("f")[1].get<std::string>();
                int if0 = GetMapIndexByName(f0);
                int if1 = GetMapIndexByName(f1);
                MyException::Assert(if0 >= 0, "if0 >= 0");
                MyException::Assert(if1 >= 0, "if1 >= 0");
                MyException::Assert(maps_[if0]->to == maps_[if1]->from, "f0->to == f1->from");
                int ig0 = -1, ig1 = -1;
                if (json_comm.at("g").size() >= 1) {
                    auto g0 = json_comm.at("g")[0].get<std::string>();
                    ig0 = GetMapIndexByName(g0);
                    MyException::Assert(ig0 >= 0, "ig0 >= 0");
                }
                if (json_comm.at("g").size() >= 2) {
                    auto g1 = json_comm.at("g")[1].get<std::string>();
                    ig1 = GetMapIndexByName(g1);
                    MyException::Assert(ig1 >= 0, "ig1 >= 0");
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
                            Logger::LogDiff(0, EnumReason::d2, name, deg, int1d{(int)i}, d2[i], 2);
                            SetRingDiffGlobal(iCw, deg, int1d{(int)i}, d2[i], 2, true, flag);
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
                            Logger::LogDiff(0, EnumReason::d2, name, deg, int1d{(int)i}, d2[i], 2);
                            SetModuleDiffGlobal(iMod, deg, int1d{(int)i}, d2[i], 2, true, flag);
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

    /*if (flag & SSFlag::homotopy)
        UpdateAllPossEinf();*/

    // VersionConvertReorderRels();
}

void Diagram::SetDeduceList(const std::vector<std::string>& cws)
{
    deduce_list_spectra_.clear();
    for (auto& name : cws) {
        auto iCw = GetIndexCwByName(name);
        deduce_list_spectra_.push_back(iCw);
    }
}

void Diagram::save(std::string diagram_name, SSFlag flag)
{
    std::map<std::string, std::string> paths;
    try {
        /* #save rings */
        auto& json_rings = js_["rings"];
        for (auto& json_ring : json_rings) {
            std::string name = json_ring["name"].get<std::string>(), path = json_ring["path"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", diagram_name, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            db.begin_transaction();
            size_t iRing = GetIndexCwByName(name).index;
            auto& ring = rings_[iRing];
            db.update_ss(table_prefix, ring.nodes_ss[1]);

            if (flag & SSFlag::pi) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);

                db.drop_and_create_pi_generators(name);
                db.save_pi_generators(name, ring.pi_gb.gen_degs(), ring.pi_gen_Einf);
                db.save_pi_gb(name, ring.pi_gb.OutputForDatabase(), GetRingGbEinf(iRing));
                db.save_pi_basis(name, ring.nodes_pi_basis.front());
                if (flag & SSFlag::pi_def) {
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
            std::string abs_path = fmt::format("{}/{}", diagram_name, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            db.begin_transaction();
            size_t iMod = GetIndexCwByName(name).index;
            auto& mod = modules_[iMod];
            db.update_ss(table_prefix, mod.nodes_ss[1]);

            if (flag & SSFlag::pi) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);
                if (flag & SSFlag::pi_def)
                    db.drop_and_create_pi_generators_mod(name);
                db.save_pi_generators_mod(name, mod.pi_gb.v_degs(), mod.pi_gen_Einf);
                db.save_pi_gb_mod(name, mod.pi_gb.OutputForDatabase(), GetModuleGbEinf(iMod));
                db.save_pi_basis_mod(name, mod.nodes_pi_basis.front());
                if (flag & SSFlag::pi_def) {
                    db.drop_and_create_pi_definitions(name);
                    db.save_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
                }
            }
            db.end_transaction();
        }

        /* #save cofseq */
        if (flag & SSFlag::cofseq) {
            auto& json_cofseqs = js_.at("cofseqs");
            std::vector<std::string> cofseq_maps;
            for (size_t i = 0; i < cofseqs_.size(); ++i) {
                auto path_cofseq = json_cofseqs[i].at("path").get<std::string>();
                auto abs_path_cofseq = fmt::format("{}/{}", diagram_name, path_cofseq);
                auto name = json_cofseqs[i].at("name").get<std::string>();
                DBSS db(abs_path_cofseq);
                db.begin_transaction();
                auto table = fmt::format("cofseq_{}", name);
                db.create_cofseq(table);
                db.save_cofseq(table, cofseqs_[i]);
                db.end_transaction();
            }
        }
    }
    catch (nlohmann::detail::exception& e) {
        fmt::print("JsonError({}): {}\n", e.id, e.what());
        throw e;
    }
}
