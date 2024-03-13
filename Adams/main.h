#ifndef MAIN_H
#define MAIN_H
#include <string>
#include <fmt/format.h>

inline const char* PROGRAM = "Adams";
inline const char* VERSION = "Version:\n  3.4.0 (2024-01-09)";
inline constexpr int DB_ADAMS_VERSION = 3;
inline constexpr std::string_view DB_VERSION_NOTES_2 = "Add t_max in version table. Change products table.";
inline constexpr std::string_view DB_VERSION_NOTES = "Add fil,from,to in version table of maps.";

void DbResVersionConvert(const char* db_filename);
namespace myio {
class Database;
}
void create_db_version(const myio::Database& db);
int get_db_t_max(const myio::Database& db);
void set_db_t_max(const myio::Database& db, int t_max);
void set_db_over(const myio::Database& db, const std::string& over);
void set_db_d2_t_max(const myio::Database& db, int t_max);
void set_db_time(const myio::Database& db);
bool IsAdamsRunning(const std::string& cmd_prefix);

/* local id for a resolution row */
inline constexpr int LOC_V_BITS = 19;
struct LocId
{
    int s, v;

    explicit LocId(int id) : s(id >> LOC_V_BITS), v(id % (1 << LOC_V_BITS)) {}
    LocId(int s, int v) : s(s), v(v) {}
    int id() const
    {
        return (s << LOC_V_BITS) | v;
    }
};

namespace ut {
inline int get(const nlohmann::json& js, std::string key, int default_)
{
    return js.contains(key) ? js.at(key).get<int>() : default_;
}
}  // namespace ut

#endif