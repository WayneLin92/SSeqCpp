#include "mylog.h"
#include <fmt/format.h>
#include <fmt/ranges.h>

std::ofstream Logger::fout_main_;
std::ofstream Logger::fout_deduce_;
std::string Logger::cmd_;
std::string Logger::line_;
std::string Logger::out_;

void Logger::SetOutMain(const char* filename)
{
    if (!fout_main_.is_open()) {
        fout_main_.open(filename, std::ofstream::app);
        if (!fout_main_.is_open()) {
            Logger::LogException(0, 0xc396a312U, "Cannot open {}\n", filename);
            throw MyException(0xc396a312U, "Cannot open file.");
        }
    }
    else {
        Logger::LogException(0, 0x31ea3cefU, "File already set for fout_main_\n");
        throw MyException(0x31ea3cefU, "File already set for fout_main_");
    }
}

void Logger::SetOutDeduce(const char* filename)
{
    if (!fout_deduce_.is_open()) {
        fout_deduce_.open(filename, std::ofstream::app);
        if (!fout_deduce_.is_open()) {
            Logger::LogException(0, 0x464ac1ddU, "Cannot open {}\n", filename);
            throw MyException(0x464ac1ddU, "Cannot open file.");
        }
        fmt::print(fout_deduce_, "{}", cmd_);
        cmd_.clear();
    }
    else {
        Logger::LogException(0, 0xa9c25373U, "File already set for fout_deduce_\n");
        throw MyException(0xa9c25373U, "File already set for fout_deduce_");
    }
}

std::string Logger::GetCmd(int argc, char** argv)
{
    return fmt::format("\nLogging start at {}\ncmd: {}\n", ut::get_time(), fmt::join(argv, argv + argc, " "));
}

void Logger::LogCmd(int argc, char** argv)
{
    cmd_ = GetCmd(argc, argv);
    fmt::print(fout_main_, "{}", cmd_);
}

void Logger::LogSummary(std::string_view category, int count)
{
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "Summary - {}: {}          \n", category, count);
    fmt::print("{}", line_);
    fmt::print(fout_main_, "{}", line_);
    if (fout_deduce_)
        fmt::print(fout_deduce_, "{}\n", line_);
}

void Logger::LogTime(std::string_view time)
{
    fmt::print(fmt::fg(fmt::color::green), "{}\n", time);
    fmt::print(fout_main_, "{}\n", time);
    if (fout_deduce_)
        fmt::print(fout_deduce_, "{}\n", time);
}

std::string_view GetReason(enumReason reason)
{
    return REASONS[size_t(reason)];
}

void Logger::LogDiff(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} d_{}{}={}\n", indent, GetReason(reason), name, deg_x, r, x, dx);
    if (depth == 0) {
        fmt::print(fmt::fg(dx.empty() ? fmt::color::white_smoke : fmt::color::light_green), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), line_);
}

void Logger::LogNullDiff(int depth, std::string_view name, alg::AdamsDeg deg_x, const alg::int1d& x, int r)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} {} d_{}{}=?\n", indent, name, deg_x, r, x);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::gray), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), line_);
}

void Logger::LogDiffInv(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_dx, const alg::int1d& x, const alg::int1d& dx, int r)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}=d_{}{}\n", indent, GetReason(reason), name, deg_dx, dx, r, x);
    if (depth == 0) {
        fmt::print(fmt::fg(x.empty() ? fmt::color::white_smoke : fmt::color::light_green), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), line_);
}

void Logger::LogDiffPerm(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} needs one more permanent cycle\n", indent, GetReason(reason), name, deg);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::light_green), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogDiffBoun(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_dx, const alg::int1d& dx)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {} is a boundary\n", indent, GetReason(reason), name, deg_dx, dx);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::light_green), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyGen(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg, size_t gen_id, const alg2::Mod& Einf)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} x_{} detected by {}\n", indent, GetReason(reason), name, deg, gen_id, Einf);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyGen(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg, size_t gen_id, const alg2::Poly& Einf)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} x_{} detected by {}\n", indent, GetReason(reason), name, deg, gen_id, Einf);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Poly& rel)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}=0\n", indent, GetReason(reason), name, deg_rel, rel);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::white_smoke), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Mod& rel)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}=0\n", indent, GetReason(reason), name, deg_rel, rel);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::white_smoke), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Poly& rel1, const algZ::Poly& rel2)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {} --> {}\n", indent, GetReason(reason), name, deg_rel, rel1, rel2);
    if (depth == 0) {
        fmt::print(fmt::fg(rel1.data.size() == rel2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyRel(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_rel, const algZ::Mod& rel1, const algZ::Mod& rel2)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {} --> {}\n", indent, GetReason(reason), name, deg_rel, rel1, rel2);
    if (depth == 0) {
        fmt::print(fmt::fg(rel1.data.size() == rel2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyMap(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, std::string_view f, size_t gen_id, const algZ::Poly& fx)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}(v_{})={}\n", indent, GetReason(reason), name, deg_x, f, gen_id, fx);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::white_smoke), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyMap(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, std::string_view f, size_t gen_id, const algZ::Poly& fx1, const algZ::Poly& fx2)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}(v_{}) = {} --> {}\n", indent, GetReason(reason), name, deg_x, f, gen_id, fx1, fx2);
    if (depth == 0) {
        fmt::print(fmt::fg(fx1.data.size() == fx2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyProd(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Poly& m, const algZ::Poly& hm1, const algZ::Poly& hm2)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}*{} = {} --> {}\n", indent, GetReason(reason), name, deg_x, h, m, hm1, hm2);
    if (depth == 0) {
        fmt::print(fmt::fg(hm1.data.size() == hm2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}

void Logger::LogHtpyProd(int depth, enumReason reason, std::string_view name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Mod& m, const algZ::Mod& hm1, const algZ::Mod& hm2)
{
    std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}*{} = {} --> {}\n", indent, GetReason(reason), name, deg_x, h, m, hm1, hm2);
    if (depth == 0) {
        fmt::print(fmt::fg(hm1.data.size() == hm2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);
}
