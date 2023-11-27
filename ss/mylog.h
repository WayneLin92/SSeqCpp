#include "algebras/database.h"
#include "algebras/myexception.h"
#include "pigroebner.h"
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ostream.h>

enum class EnumReason : uint32_t
{
    manual,    /* By other source of knowledge */
    degree,    /* For degree reason */
    deduce,    /* By deduction */
    nat,       /* By naturality */
    deduce_xx, /* By deduction of d(xy) and d(x^2) */
    deduce_xy, /* By deduction of d(xy) and d(x^2) */
    deduce_fx, /* By deduction of d(xy) and d(x^2) */
    def,       /* By definition */
    cofseq_b,  /* ss Boundary from cofseq logic */
    try1,      /* Try dx=? */
    try2,      /* Try d?=y */
    migrate,   /* Migration */
};
constexpr std::array REASONS = {"manual", "degree", "deduce", "nat", "deduce_xx", "deduce_xy", "deduce_fx", "def", "cofseq_b", "try1", "try2", "migrate"};
inline const char* INDENT = "          ";

class DbLog : public myio::Database
{
    using Statement = myio::Statement;

public:
    DbLog() = default;
    explicit DbLog(const std::string& filename) : Database(filename) {}

    void create_log() const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS log (id INTEGER PRIMARY KEY, depth TINYINT, reason TEXT, name TEXT, stem SMALLINT as (t-s), s SMALLINT, t SMALLINT, r SMALLINT, x TEXT, dx TEXT, tag TEXT);");
    }

    void InsertTag(int depth, const std::string& tag);
    void InsertError(int depth, const std::string& name, alg::AdamsDeg deg_dx, const alg::int1d& dx, int r, const std::string& tag);
    void InsertDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r);
    void InsertNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, int r);
};

/* There should be at least one global instance to close the files */
class Logger
{
private:
    static std::ofstream fout_main_;
    static std::string cmd_;
    static std::string line_;
    static DbLog db_deduce_;

private:
    static std::string GetCmd(int argc, char** argv);

public:
    Logger() {}

    static void DeleteFromLog()
    {
        db_deduce_.execute_cmd("DELETE FROM log;");
    }

    static void SetOutMain(const char* filename);
    static void SetOutDeduce(const char* filename);

    static void LogCmd(int argc, char** argv);
    static void LogSummary(const std::string& category, int count);
    static void LogTime(const std::string& time);

    static void Flush()
    {
        fout_main_.flush();
    }

    static void LogDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r);
    static void LogNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, int r);
    static void LogDiffInv(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, alg::AdamsDeg deg_dx, const alg::int1d& x, const alg::int1d& dx, int r);
    static void LogDiffBoun(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_dx, const alg::int1d& dx);
    static void LogHtpyGen(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, size_t gen_id, const alg2::Poly& Einf);
    static void LogHtpyGen(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, size_t gen_id, const alg2::Mod& Einf);
    static void LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Poly& rel);
    static void LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Mod& rel);
    static void LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Poly& rel1, const algZ::Poly& rel2);
    static void LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Mod& rel1, const algZ::Mod& rel2);
    static void LogHtpyMap(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const std::string& f, size_t gen_id, const algZ::Poly& fx);
    static void LogHtpyMap(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const std::string& f, size_t gen_id, const algZ::Poly& fx1, const algZ::Poly& fx2);
    static void LogHtpyProd(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Poly& m, const algZ::Poly& hm1, const algZ::Poly& hm2);
    static void LogHtpyProd(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Mod& m, const algZ::Mod& hm1, const algZ::Mod& hm2);

    static void LogSSException(int depth, const std::string& name, alg::AdamsDeg deg_dx, const alg::int1d& dx, int r, unsigned code, alg::AdamsDeg deg_leibniz, const alg::int1d* a_leibniz);
};

template <>
struct fmt::formatter<EnumReason>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const EnumReason reason, FormatContext& ctx)
    {

        return fmt::format_to(ctx.out(), "{}", REASONS[size_t(reason)]);
    }
};

template <>
struct fmt::formatter<alg::AdamsDeg>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const alg::AdamsDeg& deg, FormatContext& ctx)
    {
        return fmt::format_to(ctx.out(), "({}, {})", deg.stem(), deg.s);
    }
};

template <>
struct fmt::formatter<alg2::Poly>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const alg2::Poly& x, FormatContext& ctx)
    {
        return fmt::format_to(ctx.out(), "{}", x.Str());
    }
};

template <>
struct fmt::formatter<algZ::Poly>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const algZ::Poly& x, FormatContext& ctx)
    {
        return fmt::format_to(ctx.out(), "{}", x.Str());
    }
};

template <>
struct fmt::formatter<alg2::Mod>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const alg2::Mod& x, FormatContext& ctx)
    {
        return fmt::format_to(ctx.out(), "{}", x.Str());
    }
};

template <>
struct fmt::formatter<algZ::Mod>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const algZ::Mod& x, FormatContext& ctx)
    {
        return fmt::format_to(ctx.out(), "{}", x.Str());
    }
};
