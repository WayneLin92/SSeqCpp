#include "algebras/database.h"
#include "main.h"
#include <array>
#include <filesystem>
#include <fmt/core.h>
#include <map>
#include <regex>

class DbAdamsUt : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsUt(const std::string& filename) : Database(filename) {}

    int get_t_max()
    {
        if (has_table("version")) {
            try {
                return get_int("select value from version where id=817812698");
            }
            catch (MyException&) {
                return -1;
            }
        }
        return -2;
    }

    int get_timestamp()
    {
        if (has_table("version")) {
            try {
                return get_int("select value from version where id=1954841564");
            }
            catch (MyException&) {
                return -1;
            }
        }
        return -2;
    }
};

void UtStatus(const std::string& dir)
{
    std::regex is_Adams_regex("^(?:\\w+_Adams(?:_res(?:_prod|)|SS).db|)$");   /* match example: C2h4_Adams_res.db, C2h4_Adams_res_prod.db */
    std::regex is_Adams_res_regex("^(\\w+)_Adams_res.db$");                   /* match example: C2h4_Adams_res.db */
    std::regex is_Adams_res_prod_regex("^(\\w+)_Adams_res_prod.db$");         /* match example: C2h4_Adams_res_prod.db */
    std::regex is_AdamsSS_regex("^(\\w+)_AdamsSS.db$");                       /* match example: C2h4_AdamsSS.db */
    std::regex is_map_res_regex("^map(?:_\\w+|)_Adams_res_(\\w+_to_\\w+).db$"); /* match example: map_Adams_res_C2_to_S0.db */
    std::regex is_map_SS_regex("^map(?:_\\w+|)_AdamsSS_(\\w+_to_\\w+).db$");    /* match example: map_AdamsSS_C2_to_S0.db */

    std::map<std::string, std::array<std::string, 3>> table_spectra;
    std::map<std::string, std::array<std::string, 3>> table_color_spectra;
    std::map<std::string, std::array<std::string, 2>> table_maps;
    std::map<std::string, std::array<std::string, 2>> table_color_maps;

    int current_timestamp;
    {
        DbAdamsUt db("");
        current_timestamp = db.get_int("SELECT unixepoch();");
    }

    constexpr std::string_view green = "\033[38;2;0;255;0m";
    constexpr std::string_view light_green = "\033[38;2;50;255;200m";
    constexpr std::string_view light_red = "\033[38;2;255;200;200m";
    constexpr std::string_view white = "";
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        std::string filename = entry.path().filename().string();
        std::smatch match;

        if (std::regex_search(filename, match, is_Adams_res_regex); match[0].matched) {
            DbAdamsUt db(filename);
            table_spectra[match[1].str()][0] = fmt::format("{}", db.get_t_max());
            if (current_timestamp - db.get_timestamp() < (3600 * 6))
                table_color_spectra[match[1].str()][0] = green;
            else if (current_timestamp - db.get_timestamp() < (3600 * 24))
                table_color_spectra[match[1].str()][0] = light_green;
        }
        else if (std::regex_search(filename, match, is_Adams_res_prod_regex); match[0].matched) {
            DbAdamsUt db(filename);
            table_spectra[match[1].str()][1] = fmt::format("{}", db.get_t_max());
            if (current_timestamp - db.get_timestamp() < (3600 * 6))
                table_color_spectra[match[1].str()][1] = green;
            else if (current_timestamp - db.get_timestamp() < (3600 * 24))
                table_color_spectra[match[1].str()][1] = light_green;
        }
        else if (std::regex_search(filename, match, is_AdamsSS_regex); match[0].matched) {
            DbAdamsUt db(filename);
            table_spectra[match[1].str()][2] = fmt::format("{}", db.get_t_max());
            if (current_timestamp - db.get_timestamp() < (3600 * 6))
                table_color_spectra[match[1].str()][2] = green;
            else if (current_timestamp - db.get_timestamp() < (3600 * 24))
                table_color_spectra[match[1].str()][2] = light_green;
        }
        else if (std::regex_search(filename, match, is_map_res_regex); match[0].matched) {
            DbAdamsUt db(filename);
            table_maps[match[1].str()][0] = fmt::format("{}", db.get_t_max());
            if (current_timestamp - db.get_timestamp() < (3600 * 6))
                table_color_maps[match[1].str()][0] = green;
            else if (current_timestamp - db.get_timestamp() < (3600 * 24))
                table_color_maps[match[1].str()][0] = light_green;
        }
        else if (std::regex_search(filename, match, is_map_SS_regex); match[0].matched) {
            DbAdamsUt db(filename);
            table_maps[match[1].str()][1] = fmt::format("{}", db.get_t_max());
            if (current_timestamp - db.get_timestamp() < (3600 * 6))
                table_color_maps[match[1].str()][1] = green;
            else if (current_timestamp - db.get_timestamp() < (3600 * 24))
                table_color_maps[match[1].str()][1] = light_green;
        }
    }

    std::array<size_t, 4> spectra_widths = {7, 3, 4, 6};
    for (auto& [cw, t_maxes] : table_spectra) {
        spectra_widths[0] = std::max(spectra_widths[0], cw.size());
        spectra_widths[1] = std::max(spectra_widths[1], t_maxes[0].size());
        spectra_widths[2] = std::max(spectra_widths[2], t_maxes[1].size());
        spectra_widths[3] = std::max(spectra_widths[3], t_maxes[2].size());
    }
    auto fs = fmt::format("| {{}}{{:{}}}\033[0m | {{}}{{:>{}}}\033[0m | {{}}{{:>{}}}\033[0m | {{}}{{:>{}}}\033[0m |\n", spectra_widths[0], spectra_widths[1], spectra_widths[2], spectra_widths[3]);
    fmt::print(fs, light_red, "spectra", light_red, "res", light_red, "prod", light_red, "export");
    fmt::print(fs, white, "-------", white, "---", white, "----", white, "------");
    for (auto& [cw, t_maxes] : table_spectra) {
        const auto& colors = table_color_spectra[cw];
        fmt::print(fs, white, cw, colors[0], t_maxes[0], colors[1], t_maxes[1], colors[2], t_maxes[2]);
    }

    fmt::print("------------------------------------------------\n");

    std::array<size_t, 3> maps_widths = {3, 3, 6};
    for (auto& [name, t_maxes] : table_maps) {
        maps_widths[0] = std::max(maps_widths[0], name.size());
        maps_widths[1] = std::max(maps_widths[1], t_maxes[0].size());
        maps_widths[2] = std::max(maps_widths[2], t_maxes[1].size());
    }
    auto fs_map = fmt::format("| {{}}{{:{}}}\033[0m | {{}}{{:>{}}}\033[0m | {{}}{{:>{}}}\033[0m |\n", maps_widths[0], maps_widths[1], maps_widths[2]);
    fmt::print(fs_map, light_red, "map", light_red, "res", light_red, "export");
    fmt::print(fs_map, white, "---", white, "---", white, "------");
    for (auto& [name, t_maxes] : table_maps) {
        const auto& colors = table_color_maps[name];
        fmt::print(fs_map, white, name, colors[0], t_maxes[0], colors[1], t_maxes[1]);
    }
}

void UtRename(const std::string& old, const std::string& new_)
{
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
        DbAdamsUt db(fn);
        auto tables = db.get_column_str("sqlite_master", "name", fmt::format("WHERE INSTR(name,\"{}\") AND NOT INSTR(name,\"sqlite\")", old));
        for (const auto& table : tables)
            fmt::print("  {}\n", table);
    }
    if (myio::UserConfirm()) {
        for (const auto& fn : matched_filenames) {
            std::string fn_new = std::regex_replace(fn, std::regex(old), new_);
            if (std::rename(fn.c_str(), fn_new.c_str()) != 0) {
                fmt::print("Failed: {} --> {}\n", fn, fn_new);
                continue;
            }
            fmt::print("{} --> {}\n", fn, fn_new);

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

void UtAddFromTo()
{
    std::regex is_map_regex("^map(?:_\\w+|)_Adams(?:_res|SS)_(\\w+)_to_(\\w+).db$"); /* match example: map_Adams_res_C2_to_S0.db */

    std::vector<std::string> matched_filenames;
    std::string path = ".";
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        std::string filename = entry.path().filename().string();
        std::smatch match;

        if (std::regex_search(filename, match, is_map_regex); match[0].matched) {
            DbAdamsUt db(filename);
            if (db.has_table("version")) {
                try {
                    db.get_str("select value from version where id=446174262");
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

int main_status(int argc, char** argv, int& index, const char* desc)
{
    std::string dir;

    myio::CmdArg1d args = {{"dir", &dir}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtStatus(dir);
    return 0;
}

int main_rename(int argc, char** argv, int& index, const char* desc)
{
    std::string old, new_;

    myio::CmdArg1d args = {{"old", &old}, {"new", &new_}};
    myio::CmdArg1d op_args = {};
    if (int error = myio::LoadCmdArgs(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    UtRename(old, new_);
    return 0;
}

int main_add_from_to(int, char**, int&, const char*)
{
    UtAddFromTo();
    return 0;
}
