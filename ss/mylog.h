#include "algebras/myexception.h"
#include "pigroebner.h"
#include <fmt/color.h>
#include <fmt/core.h>
#include <fmt/ostream.h>

enum class enumReason : uint32_t
{
    htpy2ss,   /* Homotopy relations to ss boundaries */
    ss2htpy,   /* ss boundaries to Homotopy relations */
    degree,    /* For degree reason */
    nat,       /* By naturality */
    deduce,    /* By deduction */
    deduce_v2, /* By deduction of d(xy) */
    exact_hq,  /* By long exact sequence h*q */
    exact_ih,  /* By long exact sequence i*h */
    exact_qi,  /* By long exact sequence q*i */
    def,       /* By definition */
    try_,      /* Try */
    migrate,   /* Migration */
    manual,    /* By other source of knowledge */
};
constexpr std::array REASONS = {"htpy2ss", "ss2htpy", "degree", "nat", "deduce", "deduce_v2", "exact_hq", "exact_ih", "exact_qi", "def", "try", "migrate", "manual"};
inline const char* INDENT = "          ";

/* There should be at least one global instance to close the files */
class Logger
{
private:
    static std::ofstream fout_main_;
    static std::ofstream fout_deduce_;
    static std::string cmd_;
    static std::string line_;
    static std::string out_;

private:
    static std::string GetCmd(int argc, char** argv);

public:
    Logger() {}

    static void SetOutMain(const char* filename);
    static void SetOutDeduce(const char* filename);

    static void LogCmd(int argc, char** argv);
    static void LogSummary(std::string_view category, int count);
    static void LogTime(std::string_view time);

    /* Print the deductions in depth>=1 to fout_deduce_ */
    static void PrintDepth()
    {
        fmt::print(fout_deduce_, "{}", out_);
        out_.clear();
    }
    static void ClearDepth()
    {
        out_.clear();
    }

    static void Flush()
    {
        fout_main_.flush();
        fout_deduce_.flush();
    }

    static void LogDiff(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r);
    static void LogNullDiff(int depth, std::string_view name, alg::AdamsDeg deg_x, const alg::int1d& x, int r);
    static void LogDiffInv(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_dx, const alg::int1d& x, const alg::int1d& dx, int r);
    static void LogDiffPerm(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg);
    static void LogDiffBoun(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_dx, const alg::int1d& dx);
    static void LogHtpyGen(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg, size_t gen_id, const alg2::Poly& Einf);
    static void LogHtpyGen(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg, size_t gen_id, const alg2::Mod& Einf);
    static void LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Poly& rel);
    static void LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Mod& rel);
    static void LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Poly& rel1, const algZ::Poly& rel2);
    static void LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Mod& rel1, const algZ::Mod& rel2);
    static void LogHtpyMap(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, std::string_view f, size_t gen_id, const algZ::Poly& fx);
    static void LogHtpyMap(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, std::string_view f, size_t gen_id, const algZ::Poly& fx1, const algZ::Poly& fx2);
    static void LogHtpyProd(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Poly& m, const algZ::Poly& hm1, const algZ::Poly& hm2);
    static void LogHtpyProd(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Mod& m, const algZ::Mod& hm1, const algZ::Mod& hm2);

    template <typename... T>
    static void LogException(int depth, unsigned code, T&&... args)
    {
        std::string_view indent(INDENT, depth * 2);
        line_.clear();
        fmt::format_to(std::back_inserter(line_), "{}Error({:#x}) - ", indent, code);
        fmt::format_to(std::back_inserter(line_), args...);
        if (depth == 0) {
            fmt::print(fmt::fg(fmt::color::red), "{}", line_);
            fmt::print(fout_deduce_, "{}", line_);
            Flush();
        }
        else
            fmt::format_to(std::back_inserter(out_), "{}", line_);
    }
};

template <>
struct fmt::formatter<enumReason>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const enumReason reason, FormatContext& ctx)
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
