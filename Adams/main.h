#ifndef MAIN_H
#define MAIN_H
inline const char* PROGRAM = "Adams";
inline const char* VERSION = "Version:\n  3.1 (2023-05-07)";
inline const char* VERSION_NOTES = "Add t_max in version table. Change products table.";
inline constexpr int DB_ADAMS_VERSION = 2;

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