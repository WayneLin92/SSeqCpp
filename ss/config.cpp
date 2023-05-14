/* @file config.cpp
 * Load settings from json
 */

#include "json.h"
#include "main.h"
#include "mylog.h"
#include <fstream>

void GetAllDbNames(const std::string& diagram_name, std::vector<std::string>& names, std::vector<std::string>& paths, std::vector<int>& isRing, bool log)
{
    using json = nlohmann::json;
    json js;
    {
        std::ifstream ifs("ss.json");
        if (ifs.is_open())
            ifs >> js;
        else {
            Logger::LogException(0, 0xb8525e9bU, "File ss.json not found\n");
            throw MyException(0xb8525e9bU, "File ss.json not found");
        }
    }

    try {
        json& diagram = js["diagrams"][diagram_name];
        std::string dir = diagram["dir"].get<std::string>();
        if (log) {
            auto path_log = fmt::format("{}/{}", dir, diagram["log"].get<std::string>());
            Logger::SetOutDeduce(path_log.c_str());
        }

        /* #Load rings */
        json& json_rings = diagram["rings"];
        for (auto& json_ring : json_rings) {
            std::string name = json_ring["name"].get<std::string>(), path = json_ring["path"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            names.push_back(name);
            paths.push_back(abs_path);
            isRing.push_back(1);
        }

        /* #Load modules */
        json& json_mods = diagram["modules"];
        for (auto& json_mod : json_mods) {
            std::string name = json_mod["name"].get<std::string>(), path = json_mod["path"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            names.push_back(name);
            paths.push_back(abs_path);
            isRing.push_back(0);
        }
    }
    catch (nlohmann::detail::exception e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }
}

Diagram::Diagram(std::string diagram_name, DeduceFlag flag, bool log)
{
    using json = nlohmann::json;
    json js;
    {
        std::ifstream ifs("ss.json");
        if (ifs.is_open())
            ifs >> js;
        else {
            Logger::LogException(0, 0xb8525e9bU, "File ss.json not found\n");
            throw MyException(0xb8525e9bU, "File ss.json not found");
        }
    }

    try {
        if (diagram_name == "default")
            diagram_name = js["diagrams"]["default"].get<std::string>();
        json& diagram = js["diagrams"][diagram_name];
        std::string dir = diagram["dir"].get<std::string>();
        if (log) {
            auto path_log = fmt::format("{}/{}", dir, diagram["log"].get<std::string>());
            Logger::SetOutDeduce(path_log.c_str());
        }

        /* #Load rings */
        json& json_rings = diagram["rings"];
        rings_.reserve(json_rings.size());
        for (auto& json_ring : json_rings) {
            std::string name = json_ring["name"].get<std::string>(), path = json_ring["path"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            RingSp ring;
            ring.name = name;
            ring.basis = db.load_basis(table_prefix);
            ring.degs_basis_order_by_stem = OrderDegsByStem(ring.basis);
            ring.t_max = ring.basis.rbegin()->first.t;
            ring.nodes_ss = {db.load_basis_ss(table_prefix), {}};
            ring.nodes_ss.reserve(MAX_DEPTH + 3);
            ring.gb = Groebner(ring.t_max, {}, db.load_gb(table_prefix, DEG_MAX));

            if (flag & DeduceFlag::homotopy) {
                ring.pi_gen_Einf = db.get_column_from_str<Poly>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Poly>);
                ring.pi_gb = algZ::Groebner(ring.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb(name, DEG_MAX), true);
                ring.nodes_pi_basis.reserve(MAX_DEPTH + 2);
                if (flag & DeduceFlag::homotopy_def)
                    db.load_pi_def(name, ring.pi_gen_defs, ring.pi_gen_def_mons);
            }

            rings_.push_back(std::move(ring));
        }

        /* #Load modules */
        json& json_mods = diagram["modules"];
        modules_.reserve(json_mods.size());
        for (auto& json_mod : json_mods) {
            std::string name = json_mod["name"].get<std::string>(), path = json_mod["path"].get<std::string>();
            std::string over = json_mod["over"].get<std::string>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            std::string table_prefix = fmt::format("{}_AdamsE2", name);

            DBSS db(abs_path);
            ModSp mod;
            mod.name = name;
            mod.iRing = (size_t)GetRingIndexByName(over);
            MyException::Assert(mod.iRing != -1, "mod.iRing != -1");
            auto& ring = rings_[mod.iRing];
            mod.basis = db.load_basis_mod(table_prefix);
            mod.degs_basis_order_by_stem = OrderDegsByStem(mod.basis);
            mod.t_max = mod.basis.rbegin()->first.t;
            mod.nodes_ss = {db.load_basis_ss(table_prefix), {}};
            mod.nodes_ss.reserve(MAX_DEPTH + 1);
            Mod1d xs = db.load_gb_mod(table_prefix, DEG_MAX);
            mod.gb = GroebnerMod(&ring.gb, mod.t_max, {}, std::move(xs));

            if (flag & DeduceFlag::homotopy) {
                mod.pi_gen_Einf = db.get_column_from_str<Mod>(name + "_pi_generators", "Einf", "ORDER BY id", myio::Deserialize<Mod>);
                mod.pi_gb = algZ::GroebnerMod(&ring.pi_gb, mod.t_max, db.load_pi_gen_adamsdegs(name), db.load_pi_gb_mod(name, DEG_MAX), true);
                mod.nodes_pi_basis.reserve(MAX_DEPTH);
                if (flag & DeduceFlag::homotopy_def)
                    db.load_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
            }

            modules_.push_back(std::move(mod));
            ring.ind_mods.push_back(modules_.size() - 1);
        }

        /* #Load maps */
        json& json_maps = diagram["maps"];
        maps_.reserve(json_maps.size());
        for (auto& json_map : json_maps) {
            std::string name = json_map["name"].get<std::string>(), path = json_map["path"].get<std::string>();
            int sus = json_map["sus"].get<int>();
            std::string abs_path = fmt::format("{}/{}", dir, path);
            DBSS db(abs_path);
            std::string from = json_map["from"].get<std::string>(), to = json_map["to"].get<std::string>();
            std::string table = fmt::format("map_AdamsE2_{}_to_{}", from, to);
            Map map;
            map.name = name;
            if (GetRingIndexByName(from) != -1) {
                size_t index_from = (size_t)GetRingIndexByName(from);
                size_t index_to = (size_t)GetRingIndexByName(to);
                MyException::Assert(index_from != -1 && index_to != -1, "index_from != -1 && index_to != -1");
                auto images = db.get_column_from_str<Poly>(table, "map", "ORDER BY id", myio::Deserialize<Poly>);
                map.map = MapRing2Ring{index_from, index_to, std::move(images)};
                rings_[index_from].ind_maps.push_back(maps_.size());
            }
            else {
                size_t index_from = (size_t)GetModuleIndexByName(from);
                MyException::Assert(index_from != -1, "index_from != -1");
                if (GetRingIndexByName(to) != -1) {
                    size_t index_to = (size_t)GetRingIndexByName(to);
                    MyException::Assert(index_to != -1, "index_to != -1");
                    auto images = db.get_column_from_str<Poly>(table, "map", "", myio::Deserialize<Poly>);
                    map.map = MapMod2Ring{index_from, index_to, sus, std::move(images)};
                }
                else {
                    size_t index_to = (size_t)GetModuleIndexByName(to);
                    MyException::Assert(index_to != -1, "index_to != -1");
                    auto images = db.get_column_from_str<Mod>(table, "map", "", myio::Deserialize<Mod>);
                    map.map = MapMod2Mod{index_from, index_to, sus, std::move(images)};
                }
                modules_[index_from].ind_maps.push_back(maps_.size());
            }
            maps_.push_back(std::move(map));
        }
    }
    catch (nlohmann::detail::exception e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }

    /*if (flag & DeduceFlag::homotopy)
        UpdateAllPossEinf();*/

    // VersionConvertReorderRels();
}

void Diagram::save(std::string diagram_name, DeduceFlag flag)
{
    using json = nlohmann::json;
    json js;
    {
        std::ifstream ifs("ss.json");
        if (ifs.is_open())
            ifs >> js;
        else {
            Logger::LogException(0, 0x10ebb3c3, "File ss.json not found\n");
            throw MyException(0x10ebb3c3, "File ss.json not found");
        }
    }
    std::map<std::string, std::string> paths;
    try {
        if (diagram_name == "default")
            diagram_name = js["diagrams"]["default"].get<std::string>();
        json& diagram = js["diagrams"][diagram_name];
        std::string dir = diagram["dir"].get<std::string>();

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
            db.update_basis_ss(table_prefix, ring.nodes_ss[1]);

            if (flag & DeduceFlag::homotopy) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);

                db.drop_and_create_pi_generators(name);
                db.save_pi_generators(name, ring.pi_gb.gen_degs(), ring.pi_gen_Einf);
                db.save_pi_gb(name, ring.pi_gb.OutputForDatabase(), GetRingGbEinf(iRing));
                db.save_pi_basis(name, ring.nodes_pi_basis.front());
                if (flag & DeduceFlag::homotopy_def) {
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
            db.update_basis_ss(table_prefix, mod.nodes_ss[1]);

            if (flag & DeduceFlag::homotopy) {
                db.drop_and_create_pi_relations(name);
                db.drop_and_create_pi_basis(name);
                if (flag & DeduceFlag::homotopy_def)
                    db.drop_and_create_pi_generators_mod(name);
                db.save_pi_generators_mod(name, mod.pi_gb.v_degs(), mod.pi_gen_Einf);
                db.save_pi_gb_mod(name, mod.pi_gb.OutputForDatabase(), GetModuleGbEinf(iMod));
                db.save_pi_basis_mod(name, mod.nodes_pi_basis.front());
                if (flag & DeduceFlag::homotopy_def) {
                    db.drop_and_create_pi_definitions(name);
                    db.save_pi_def(name, mod.pi_gen_defs, mod.pi_gen_def_mons);
                }
            }
            db.end_transaction();
        }

        ////TODO: save pi_map
    }
    catch (nlohmann::detail::exception e) {
        Logger::LogException(0, e.id, "{}\n", e.what());
        throw e;
    }
}
