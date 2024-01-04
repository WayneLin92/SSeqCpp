#include "algebras/database.h"
#include "algebras/myio.h"
#include "algebras/utility.h"
#include "main.h"
#include <array>
#include <filesystem>
#include <fmt/core.h>
#include <map>
#include <regex>

int get_db_t_verified(const myio::Database& db);

void create_db_version(const myio::Database& db)
{
    db.execute_cmd("CREATE TABLE IF NOT EXISTS version (id INTEGER PRIMARY KEY, name TEXT, value);");
}

int get_db_t_max(const myio::Database& db)
{
    try {
        if (db.has_table("version")) {
            try {
                return db.get_int("select value from version where id=817812698");
            }
            catch (MyException&) {
                return -1;
            }
        }
    }
    catch (MyException&) {
        return -3;
    }
    return -2;
}

void set_db_t_max(const myio::Database& db, int t_max)
{
    myio::Statement stmt(db, "INSERT INTO version (id, name, value) VALUES (817812698, \"t_max\", ?1) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
    stmt.bind_and_step(t_max);
}

void set_db_time(const myio::Database& db)
{
    myio::Statement stmt(db, "INSERT INTO version (id, name, value) VALUES (1954841564, \"timestamp\", unixepoch()) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
    stmt.step_and_reset();
}

class DbAdamsUt : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsUt(const std::string& filename) : Database(filename) {}

    int get_timestamp()
    {
        try {
            if (has_table("version")) {
                try {
                    return get_int("select value from version where id=1954841564");
                }
                catch (MyException&) {
                    return -1;
                }
            }
        }
        catch (MyException&) {
            return -3;
        }
        return -2;
    }
};

void UtStatus(const std::string& dir, bool sorted)
{
    std::regex is_Adams_res_regex("^(\\w+)_Adams_res.db$");                     /* match example: C2h4_Adams_res.db */
    std::regex is_Adams_res_prod_regex("^(\\w+)_Adams_res_prod.db$");           /* match example: C2h4_Adams_res_prod.db */
    std::regex is_AdamsSS_regex("^(\\w+)_AdamsSS.db$");                         /* match example: C2h4_AdamsSS.db */
    std::regex is_map_res_regex("^map(?:_\\w+|)_Adams_res_(\\w+__\\w+).db$"); /* match example: map_Adams_res_C2__S0.db */
    std::regex is_map_SS_regex("^map(?:_\\w+|)_AdamsSS_(\\w+__\\w+).db$");    /* match example: map_AdamsSS_C2__S0.db */
    std::smatch match;

    std::map<std::string, int> timestamps_spectra;
    std::map<std::string, std::array<std::string, 3>> table_spectra;
    std::map<std::string, std::array<std::string, 3>> table_color_spectra;
    std::map<std::string, int> timestamps_maps;
    std::map<std::string, std::array<std::string, 2>> table_maps;
    std::map<std::string, std::array<std::string, 2>> table_color_maps;

    int current_timestamp;
    {
        DbAdamsUt db("");
        current_timestamp = db.get_int("SELECT unixepoch();");
    }

    constexpr std::string_view green = "\033[38;2;0;255;0m";
    constexpr std::string_view blue = "\033[38;2;50;128;255m";
    constexpr std::string_view light_red = "\033[38;2;255;200;200m";
    constexpr std::string_view white = "";
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::string filename = entry.path().filename().string();
        std::string filepath = entry.path().string();

        {
            std::string name;
            int index;
            if (std::regex_search(filename, match, is_Adams_res_regex); match[0].matched) {
                name = match[1].str();
                index = 0;
            }
            else if (std::regex_search(filename, match, is_Adams_res_prod_regex); match[0].matched) {
                name = match[1].str();
                index = 1;
            }
            else if (std::regex_search(filename, match, is_AdamsSS_regex); match[0].matched) {
                name = match[1].str();
                index = 2;
            }
            if (!name.empty()) {
                DbAdamsUt db(filepath);
                table_spectra[name][index] = fmt::format("{}", get_db_t_max(db));
                int timestamp = db.get_timestamp();
                timestamps_spectra[name] = std::max(timestamps_spectra[name], timestamp);
                if (current_timestamp - timestamp < (3600 * 6))
                    table_color_spectra[name][index] = green;
                else if (current_timestamp - timestamp < (3600 * 24))
                    table_color_spectra[name][index] = blue;
            }
        }

        {
            std::string name;
            int index;
            if (std::regex_search(filename, match, is_map_res_regex); match[0].matched) {
                name = match[1].str();
                index = 0;
            }
            else if (std::regex_search(filename, match, is_map_SS_regex); match[0].matched) {
                name = match[1].str();
                index = 1;
            }
            if (!name.empty()) {
                DbAdamsUt db(filepath);
                table_maps[name][index] = fmt::format("{}", get_db_t_max(db));
                int timestamp = db.get_timestamp();
                timestamps_maps[name] = std::max(timestamps_maps[name], timestamp);
                if (current_timestamp - timestamp < (3600 * 6))
                    table_color_maps[name][index] = green;
                else if (current_timestamp - timestamp < (3600 * 24))
                    table_color_maps[name][index] = blue;
            }
        }
    }
    auto names_spectra = ut::get_keys(table_spectra);
    auto names_maps = ut::get_keys(table_maps);
    if (sorted) {
        std::sort(names_spectra.begin(), names_spectra.end(), [&timestamps_spectra](const std::string& name1, const std::string& name2) { return timestamps_spectra.at(name1) > timestamps_spectra.at(name2); });
        std::sort(names_maps.begin(), names_maps.end(), [&timestamps_maps](const std::string& name1, const std::string& name2) { return timestamps_maps.at(name1) > timestamps_maps.at(name2); });
    }

    std::array<size_t, 4> spectra_widths = {7, 3, 4, 6};
    for (auto& cw : names_spectra) {
        auto& t_maxes = table_spectra.at(cw);
        spectra_widths[0] = std::max(spectra_widths[0], cw.size());
        spectra_widths[1] = std::max(spectra_widths[1], t_maxes[0].size());
        spectra_widths[2] = std::max(spectra_widths[2], t_maxes[1].size());
        spectra_widths[3] = std::max(spectra_widths[3], t_maxes[2].size());
    }
    auto fs = fmt::format("| {{}}{{:{}}}\033[0m | {{}}{{:>{}}}\033[0m | {{}}{{:>{}}}\033[0m | {{}}{{:>{}}}\033[0m |\n", spectra_widths[0], spectra_widths[1], spectra_widths[2], spectra_widths[3]);
    fmt::print(fs, light_red, "spectra", light_red, "res", light_red, "prod", light_red, "export");
    fmt::print(fs, white, "-------", white, "---", white, "----", white, "------");
    for (auto& cw : names_spectra) {
        auto& t_maxes = table_spectra.at(cw);
        const auto& colors = table_color_spectra[cw];
        fmt::print(fs, white, cw, colors[0], t_maxes[0], colors[1], t_maxes[1], colors[2], t_maxes[2]);
    }
    fmt::print("----------------------------------------------\n");

    std::array<size_t, 3> maps_widths = {3, 3, 6};
    for (auto& name : names_maps) {
        auto& t_maxes = table_maps.at(name);
        maps_widths[0] = std::max(maps_widths[0], name.size());
        maps_widths[1] = std::max(maps_widths[1], t_maxes[0].size());
        maps_widths[2] = std::max(maps_widths[2], t_maxes[1].size());
    }
    auto fs_map = fmt::format("| {{}}{{:{}}}\033[0m | {{}}{{:>{}}}\033[0m | {{}}{{:>{}}}\033[0m |\n", maps_widths[0], maps_widths[1], maps_widths[2]);
    fmt::print(fs_map, light_red, "map", light_red, "res", light_red, "export");
    fmt::print(fs_map, white, "---", white, "---", white, "------");
    for (auto& name : names_maps) {
        auto& t_maxes = table_maps.at(name);
        const auto& colors = table_color_maps[name];
        fmt::print(fs_map, white, name, colors[0], t_maxes[0], colors[1], t_maxes[1]);
    }
    fmt::print("----------------------------------------------\n");
}

void UtVerifyStatus(const std::string& dir, bool sorted)
{
    std::regex is_map_res_regex("^map(?:_\\w+|)_Adams_res_(\\w+__\\w+).db$"); /* match example: map_Adams_res_C2__S0.db */
    std::smatch match;

    std::map<std::string, int> timestamps_maps;
    std::map<std::string, std::array<std::string, 2>> table_maps;
    std::map<std::string, std::array<std::string, 2>> table_color_maps;

    int current_timestamp;
    {
        DbAdamsUt db("");
        current_timestamp = db.get_int("SELECT unixepoch();");
    }

    constexpr std::string_view green = "\033[38;2;0;255;0m";
    constexpr std::string_view blue = "\033[38;2;50;128;255m";
    constexpr std::string_view light_red = "\033[38;2;255;200;200m";
    constexpr std::string_view white = "";
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::string filename = entry.path().filename().string();
        std::string filepath = entry.path().string();

        {
            std::string name;
            if (std::regex_search(filename, match, is_map_res_regex); match[0].matched) {
                name = match[1].str();
            }
            if (!name.empty()) {
                DbAdamsUt db(filepath);
                table_maps[name][0] = fmt::format("{}", get_db_t_max(db));
                table_maps[name][1] = fmt::format("{}", get_db_t_verified(db));
                int timestamp = db.get_timestamp();
                timestamps_maps[name] = timestamp;
                if (current_timestamp - timestamp < (3600 * 6))
                    table_color_maps[name][1] = green;
                else if (current_timestamp - timestamp < (3600 * 24))
                    table_color_maps[name][1] = blue;
            }
        }
    }
    auto names_maps = ut::get_keys(table_maps);
    if (sorted) {
        std::sort(names_maps.begin(), names_maps.end(), [&timestamps_maps](const std::string& name1, const std::string& name2) { return timestamps_maps.at(name1) > timestamps_maps.at(name2); });
    }

    std::array<size_t, 3> maps_widths = {3, 3, 6};
    for (auto& name : names_maps) {
        auto& t_maxes = table_maps.at(name);
        maps_widths[0] = std::max(maps_widths[0], name.size());
        maps_widths[1] = std::max(maps_widths[1], t_maxes[0].size());
        maps_widths[2] = std::max(maps_widths[2], t_maxes[1].size());
    }
    auto fs_map = fmt::format("| {{}}{{:{}}}\033[0m | {{}}{{:>{}}}\033[0m | {{}}{{:>{}}}\033[0m |\n", maps_widths[0], maps_widths[1], maps_widths[2]);
    fmt::print(fs_map, light_red, "map", light_red, "res", light_red, "verify");
    fmt::print(fs_map, white, "---", white, "---", white, "------");
    for (auto& name : names_maps) {
        auto& t_maxes = table_maps.at(name);
        const auto& colors = table_color_maps[name];
        fmt::print(fs_map, white, name, colors[0], t_maxes[0], colors[1], t_maxes[1]);
    }
    fmt::print("----------------------------------------------\n");
}

void UtExport(const std::string& dir)
{
    std::regex is_Adams_res_prod_regex("^(\\w+)_Adams_res_prod.db$");           /* match example: C2h4_Adams_res_prod.db */
    std::regex is_AdamsSS_regex("^(\\w+)_AdamsSS.db$");                         /* match example: C2h4_AdamsSS.db */
    std::regex is_map_res_regex("^map(?:_\\w+|)_Adams_res_(\\w+__\\w+).db$"); /* match example: map_Adams_res_C2__S0.db */
    std::regex is_map_SS_regex("^map(?:_\\w+|)_AdamsSS_(\\w+__\\w+).db$");    /* match example: map_AdamsSS_C2__S0.db */
    std::smatch match;

    std::map<std::string, std::array<int, 2>> table_spectra, table_maps;

    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::string filename = entry.path().filename().string();
        std::string filepath = entry.path().string();

        {
            std::string name;
            int index;
            if (std::regex_search(filename, match, is_Adams_res_prod_regex); match[0].matched) {
                name = match[1].str();
                index = 0;
            }
            else if (std::regex_search(filename, match, is_AdamsSS_regex); match[0].matched) {
                name = match[1].str();
                index = 1;
            }
            if (!name.empty()) {
                DbAdamsUt db(filepath);
                table_spectra[name][index] = get_db_t_max(db);
            }
        }

        {
            std::string name;
            int index;
            if (std::regex_search(filename, match, is_map_res_regex); match[0].matched) {
                name = match[1].str();
                index = 0;
            }
            else if (std::regex_search(filename, match, is_map_SS_regex); match[0].matched) {
                name = match[1].str();
                index = 1;
            }
            if (!name.empty()) {
                DbAdamsUt db(filepath);
                table_maps[name][index] = get_db_t_max(db);
            }
        }
    }
    auto names_spectra = ut::get_keys(table_spectra);
    auto names_maps = ut::get_keys(table_maps);
    ut::RemoveIf(names_spectra, [&table_spectra](std::string& cw) { return table_spectra.at(cw)[0] <= table_spectra.at(cw)[1]; });
    ut::RemoveIf(names_maps, [&table_maps](std::string& map) { return table_maps.at(map)[0] <= table_maps.at(map)[1]; });

    for (auto& cw : names_spectra)
        fmt::print("./Adams export_mod {} S0 200\n", cw);
    for (auto& map : names_maps)
        fmt::print("./Adams export_map {} 200\n", std::regex_replace(map, std::regex("__"), " "));
    fmt::print("\nrm -r download\nmkdir download\n");
    fmt::print("cp ");
    for (auto& cw : names_spectra)
        fmt::print("{}_AdamsSS.db ", cw);
    fmt::print("download\n");
    fmt::print("cp ");
    for (auto& map : names_maps)
        fmt::print("map_AdamsSS_{}.db ", map);
    fmt::print("download\n");
}

void UtRename(const std::string& old, const std::string& new_)
{
    std::regex is_db_regex(".db$"); /* match example: *.db */
    std::smatch match;

    std::vector<std::string> matched_filenames;
    std::string path = ".";
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        std::string filename = entry.path().filename().string();
        if (filename.find(old) != std::string::npos)
            matched_filenames.push_back(filename);
    }
    fmt::print("Affected items:\n");
    for (const auto& fn : matched_filenames) {
        fmt::print("{}\n", fn);
        if (std::regex_search(fn, match, is_db_regex); match[0].matched) {
            DbAdamsUt db(fn);
            auto tables = db.get_column_str("sqlite_master", "name", fmt::format("WHERE INSTR(name,\"{}\") AND NOT INSTR(name,\"sqlite\")", old));
            for (const auto& table : tables)
                fmt::print("  {}\n", table);
            if (db.has_table("version")) {
                auto names = db.get_column_str("version", "name", fmt::format("WHERE INSTR(value,\"{}\")", old));
                auto values = db.get_column_str("version", "value", fmt::format("WHERE INSTR(value,\"{}\")", old));
                for (size_t i = 0; i < names.size(); ++i)
                    fmt::print("  version: {}={}\n", names[i], values[i]);
            }
        }
    }
    if (myio::UserConfirm()) {
        for (const auto& fn : matched_filenames) {
            std::string fn_new = std::regex_replace(fn, std::regex(old), new_);
            if (std::rename(fn.c_str(), fn_new.c_str()) != 0) {
                fmt::print("Failed: {} --> {}\n", fn, fn_new);
                continue;
            }
            fmt::print("{} --> {}\n", fn, fn_new);

            if (std::regex_search(fn, match, is_db_regex); match[0].matched) {
                DbAdamsUt db(fn_new);
                auto tables = db.get_column_str("sqlite_master", "name", fmt::format("WHERE INSTR(name,\"{}\") AND NOT INSTR(name,\"sqlite\")", old));
                for (const auto& table : tables) {
                    std::string table_new = std::regex_replace(table, std::regex(old), new_);
                    db.rename_table(table, table_new);
                    fmt::print("  {} --> {}\n", table, table_new);
                }
            }
        }
    }
}

void UtAppendTmaxToFilename(const std::string& dir)
{
    std::regex is_AdamsSS_db_regex("(?!\\w+_t\\d+.db)(\\w*AdamsSS\\w*).db$"); /* match example: *.db */
    std::smatch match;

    std::vector<std::string> old, new_;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::string filename = entry.path().filename().string();
        if (std::regex_search(filename, match, is_AdamsSS_db_regex); match[0].matched) {
            old.push_back(dir + "/" + filename);
            DbAdamsUt db(dir + "/" + filename);
            int t_max = get_db_t_max(db);
            new_.push_back(dir + "/" + fmt::format("{}_t{}.db", match[1].str(), t_max));
        }
    }
    for (size_t i = 0; i < old.size(); ++i) {
        if (std::rename(old[i].c_str(), new_[i].c_str()) != 0) {
            fmt::print("Failed: {} --> {}\n", old[i], new_[i]);
        }
    }
}

void UtPrintSSJson(const std::string& dir)
{
    std::regex is_AdamsSS_regex("^(\\w+)_AdamsSS_t\\d+.db$");              /* match example: C2h4_AdamsSS_t100.db */
    std::regex is_map_SS_regex("^map_AdamsSS_(\\w+)__(\\w+)_t\\d+.db$"); /* match example: map_AdamsSS_C2__S0.db */
    std::smatch match, match1;

    std::map<std::string, std::array<int, 2>> table_maps; /* name -> [t_max, index] */
    std::vector<std::string> outputs;
    std::vector<std::string> paths;
    std::vector<int> toBeRemoved;

    int index = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::string filename = entry.path().filename().string();
        if (std::regex_search(filename, match, is_AdamsSS_regex); match[0].matched) {
            // DbAdamsUt db(dir + "/" + filename);
            // fmt::print("{{ \"name\": \"{}\", \"path\": \"{}\", \"over\": \"S0\", \"deduce\": \"on\" }},\n", match[1].str(), filename);  ////
        }
        else if (std::regex_search(filename, match, is_map_SS_regex); match[0].matched) {
            try {
                DbAdamsUt db(dir + "/" + filename);
                int t_max = db.get_int("select value from version where id=817812698");
                auto from = db.get_str("select value from version where id=446174262");
                auto to = db.get_str("select value from version where id=1713085477");
                int sus = db.get_int("select value from version where id=1585932889");
                int fil = db.get_int("select value from version where id=651971502");
                std::string str_sus = sus ? fmt::format(", \"sus\": {}", sus) : "";
                std::string str_fil = fil ? fmt::format(", \"fil\": {}", fil) : "";
                std::string name = fmt::format("{}__{}", match[1].str(), match[2].str());
                if (ut::has(table_maps, name)) {
                    if (table_maps.at(name)[0] < t_max) {
                        toBeRemoved.push_back(table_maps.at(name)[1]);
                        table_maps.at(name) = {t_max, index};
                    }
                }
                else
                    table_maps[name] = {t_max, index};
                std::string tpl = "{{ \"name\": \"{}\", \"display\": \"{} -> {}\", \"path\": \"{}\", \"from\": \"{}\", \"to\": \"{}\"{}{}, \"t_max\": {} }},\n";
                outputs.push_back(fmt::format(tpl, name, match[1].str(), match[2].str(), filename, from, to, str_sus, str_fil, t_max));  ////
                paths.push_back(filename);
                ++index;
            }
            catch (MyException&) {
                fmt::print("Keep {}\n", filename);
            }
        }
    }

    for (auto& [name, t_max_index] : table_maps)
        fmt::print("{}", outputs[t_max_index[1]]);
    fmt::print("\n");
    for (int i : toBeRemoved)
        fmt::print("Remove {}\n", paths[i]);
}

void UtAddFromTo(const std::string& dir)
{
    std::regex is_map_regex("^map_Adams(?:_res|SS)_(\\w+)__(\\w+?)(?:_t\\d+|).db$"); /* match example: map_Adams_res_C2__S0.db */

    std::vector<std::string> matched_filenames;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::string filename = entry.path().filename().string();
        std::smatch match;

        if (std::regex_search(filename, match, is_map_regex); match[0].matched) {
            DbAdamsUt db(dir + "/" + filename);
            if (db.has_table("version")) {
                try {
                    db.get_str("select value from version where id=446174262"); /* from */
                }
                catch (MyException&) {
                    myio::Statement stmt(db, "INSERT INTO version (id, name, value) VALUES (?1, ?2, ?3) ON CONFLICT(id) DO UPDATE SET value=excluded.value;");
                    stmt.bind_and_step(446174262, std::string("from"), match[1].str());
                    stmt.bind_and_step(1713085477, std::string("to"), match[2].str());
                    fmt::print("{}: {} -> {}\n", filename, match[1].str(), match[2].str());
                }
            }
        }
    }
}

void UtAddTMax(const std::string& db_filename, int t_max)
{
    myio::AssertFileExists(db_filename);
    DbAdamsUt db(db_filename);
    if (!db.has_table("version"))
        create_db_version(db);
    set_db_t_max(db, t_max);
}

int main_status(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";
    myio::string1d options;

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}, {"options", &options}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    bool sorted = false;
    for (auto& op : options) {
        if (op == "sorted")
            sorted = true;
    }

    UtStatus(dir, sorted);
    return 0;
}

int main_verify_status(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";
    myio::string1d options;

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}, {"options", &options}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;
    bool sorted = false;
    for (auto& op : options) {
        if (op == "sorted")
            sorted = true;
    }

    UtVerifyStatus(dir, sorted);
    return 0;
}

int main_ut_export(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";
    myio::string1d options;

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtExport(dir);
    return 0;
}

int main_ut_rename(int argc, char** argv, int& index, const char* desc)
{
    std::string old, new_;

    myio::CmdArg1d args = {{"old", &old}, {"new", &new_}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtRename(old, new_);
    return 0;
}

int main_ut_app_t_max(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtAppendTmaxToFilename(dir);
    return 0;
}

int main_ut_ss_json(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtPrintSSJson(dir);
    return 0;
}

int main_ut_add_from_to(int argc, char** argv, int& index, const char* desc)
{
    std::string dir = ".";

    myio::CmdArg1d args = {};
    myio::CmdArg1d op_args = {{"dir", &dir}};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtAddFromTo(dir);
    return 0;
}

int main_ut_add_t_max(int argc, char** argv, int& index, const char* desc)
{
    std::string db_filename;
    int t_max = 0;

    myio::CmdArg1d args = {{"db", &db_filename}, {"t_max", &t_max}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtAddTMax(db_filename, t_max);
    return 0;
}

int main_ut(int argc, char** argv, int& index, const char* desc)
{
    myio::SubCmdArg1d subcmds = {
        {"export", "Print commands for exporting new results", main_ut_export},
        {"rename", "For compatibility: rename files", main_ut_rename},
        {"app_t_max", "append t_max to filenames", main_ut_app_t_max},
        {"ss_json", "print info for ss.json", main_ut_ss_json},
        {"add_from_to", "For compatibility: add from to info to databases of maps", main_ut_add_from_to},
        {"add_t_max", "For compatibility: add t_max info to the database", main_ut_add_t_max},
    };
    if (int error = myio::LoadSubCmd(argc, argv, index, PROGRAM, desc, VERSION, subcmds))
        return error;
    return 0;
}
