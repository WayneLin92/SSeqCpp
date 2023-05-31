#ifndef MAIN_H
#define MAIN_H
#include <string>

inline const char* PROGRAM = "Adams";
inline const char* VERSION = "Version:\n  3.1.1 (2023-05-29)";
inline constexpr int DB_ADAMS_VERSION = 3;
inline constexpr std::string_view DB_VERSION_NOTES_2 = "Add t_max in version table. Change products table.";
inline constexpr std::string_view DB_VERSION_NOTES = "Add fil,from,to in version table of maps.";

void DbResVersionConvert(const char* db_filename);

//#define MYDEPLOY

/* local id for a resolution row */
inline constexpr int LOC_V_BITS = 19;
struct LocId
{
    int s, v;

    explicit LocId(int id) : s(id >> LOC_V_BITS), v(id % (1 << LOC_V_BITS)) {}
    LocId(int s_, int v_) : s(s_), v(v_) {}
    int id() const
    {
        return (s << LOC_V_BITS) | v;
    }
};


#endif