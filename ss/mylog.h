#ifndef MYLOG_H
#define MYLOG_H

#include "algebras/database.h"
#include "algebras/myexception.h"
#include "pigroebner.h"
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ostream.h>
#include <tuple>

enum class SSFlag : uint32_t;
enum class EnumReason : uint32_t
{
    manual,     /* By other source of knowledge */
    degree,     /* For degree reason */
    degree2,    /* For degree reason */
    deduce,     /* By deduction */
    deduce2,    /* Find ss source by deduction */
    dd_cof,     /* By deduction in cofseq */
    dd_cof2,    /* By deduction in cofseq */
    cofseq_in,  /* ss permanent cycle by deduction in cofseq */
    cofseq_out, /* ss Boundary from zero cofseq */
    nat,        /* By naturality */
    syn,     /* By naturality of synthetic */
    synext,     /* By naturality of synthetic */
    synext_p,   /* Permanent cycle by synthetic */
    deduce_xx,  /* d(x) -> d(x^2) */
    deduce_xy,  /* d(x)=? -> d(xy) */
    deduce_fx,  /* By deduction of d(xy) and d(x^2) */
    comm,       /* By commutativity */
    def,        /* By definition */
    try1,       /* Try dx=? */
    try2,       /* Try d?=y */
    migrate,    /* Migration */
    d2,         /* By Adams d2 computation */
};
constexpr std::array REASONS = {"manual", "degree",    "degree2",   "deduce",    "deduce2",   "deduce_cs", "deduce_cs2", "cs_in", "cs_out", "nat",     "syn_nat",
                                "syn_cs", "syn_cs_in", "deduce_xx", "deduce_xy", "deduce_fx", "comm",      "def",       "try",   "try2",   "migrate", "d2"};
constexpr std::array REASONS_DB = {"M", "G", "GI", "D", "DI", "D", "DI", "ToCs", "OutCsI", "N", "Syn", "SynCs", "SynCsIn", "XX", "XY", "FX", "CsCm", "Def", "T", "TI", "Mg", "d2"};
inline const char* INDENT = "          ";
static_assert(REASONS.size() == REASONS_DB.size());

class DbLog : public myio::Database
{
    using Statement = myio::Statement;

public:
    DbLog() = default;
    explicit DbLog(const std::string& filename) : Database(filename) {}
    virtual ~DbLog() {}

    void create_log() const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS log (id INTEGER PRIMARY KEY, depth TINYINT, reason TEXT, name TEXT, stem SMALLINT as (t-s), s SMALLINT, t SMALLINT, r SMALLINT, x TEXT, dx TEXT, info TEXT);");
        execute_cmd("CREATE TABLE IF NOT EXISTS exclusions (id INTEGER PRIMARY KEY);");
    }
    void drop_and_create_minimal() const
    {
        drop_table("minimal");
        execute_cmd("CREATE TABLE IF NOT EXISTS minimal (id INTEGER PRIMARY KEY);");
    }

    void InsertInfo(int depth, const std::string& info);
    void InsertExclusion(int id);
    void InsertDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, int r, const alg::int1d& x, const alg::int1d& dx, const std::string& info);
    void InsertDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r);
    void InsertNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, int r);
};

/* return (name, deg, r, x, dx, isCs, isDInv) */
inline auto columns_diff(myio::Statement& stmt)
{
    return std::make_tuple(stmt.column_int(0), stmt.column_str(1), alg::AdamsDeg(stmt.column_int(2), stmt.column_int(3)), stmt.column_int(4), myio::Deserialize<alg::int1d>(stmt.column_str(5)), myio::Deserialize<alg::int1d>(stmt.column_str(6)),
                           (bool)stmt.column_int(7), (bool)stmt.column_int(8));
}

/* There should be at least one global instance to close the files */
class Logger
{
private:
    static std::ofstream fout_main_;
    static std::string cmd_;
    static std::string line_;

public:
    static DbLog db_;

private:
    static std::string GetCmd(int argc, char** argv);

public:
    Logger() {}

    static void Reset();
    static void SetOutMain(const char* filename);
    static void SetOutDeduce(const char* filename);

    static void LogCmd(int argc, char** argv);
    static void LogSummary(const std::string& category, int count);
    static void LogTime(const std::string& time);

    static int GetCheckpoint();
    static void RollBackToCheckpoint(int checkpt);
    static void LogExclusion(int id);
    static void RollBackToExclusionsCheckpoint(int checkpt);

    static void LogDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, int r, const alg::int1d& x, const alg::int1d& dx, const std::string& proof, SSFlag flag);
    static void LogDiffInv(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_dx, int r, const alg::int1d& x, const alg::int1d& dx, const std::string& proof, SSFlag flag);
    static void LogGreyDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, int r, const alg::int1d& x, const alg::int1d& dx, SSFlag flag);
    static void LogNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, int r, const alg::int1d& x, SSFlag flag);
};

template <>
struct fmt::formatter<EnumReason>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx) const
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const EnumReason reason, FormatContext& ctx) const
    {

        return fmt::format_to(ctx.out(), "{}", REASONS[size_t(reason)]);
    }
};

#endif