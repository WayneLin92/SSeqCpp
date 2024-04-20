#include "mylog.h"
#include <fmt/format.h>
#include <fmt/ranges.h>

std::ofstream Logger::fout_main_;
std::string Logger::cmd_;
std::string Logger::line_;
DbLog Logger::db_deduce_;
int Logger::id_checkpoint_;

std::string_view GetReason(EnumReason reason)
{
    return REASONS[size_t(reason)];
}

/* -1: start time
 * -2: end time
 * -3:
 */
void DbLog::InsertTag(int depth, const std::string& tag)
{
    Statement stmt(*this, "INSERT INTO log (depth, tag) VALUES (?1, ?2);");
    stmt.bind_and_step(depth, std::string(tag));
}

void DbLog::InsertError(int depth, const std::string& name, alg::AdamsDeg deg_dx, const alg::int1d& dx, int r, const std::string& tag)
{
    Statement stmt(*this, "INSERT INTO log (depth, reason, name, s, t, r, dx, tag) VALUES (?1, \"Error\", ?2, ?3, ?4, ?5, ?6, ?7)");
    stmt.bind_and_step(depth, name, deg_dx.s, deg_dx.t, r, myio::Serialize(dx), tag);
}

void DbLog::InsertError(int depth, const std::string& tag)
{
    Statement stmt(*this, "INSERT INTO log (depth, reason, tag) VALUES (?1, \"Error\", ?2)");
    stmt.bind_and_step(depth, tag);
}

void DbLog::InsertDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r)
{
    Statement stmt(*this, "INSERT INTO log (depth, reason, name, s, t, r, x, dx) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)");
    stmt.bind_and_step(depth, std::string(GetReason(reason)), name, deg_x.s, deg_x.t, r, myio::Serialize(x), myio::Serialize(dx));
}

void DbLog::InsertNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, int r)
{
    Statement stmt(*this, "INSERT INTO log (depth, name, s, t, r, x) VALUES (?1, ?2, ?3, ?4, ?5, ?6)");
    stmt.bind_and_step(depth, name, deg_x.s, deg_x.t, r, myio::Serialize(x));
}

void Logger::DeleteFromLog()
{
    db_deduce_.execute_cmd("DELETE FROM log");
    db_deduce_.end_transaction();
    db_deduce_.execute_cmd("VACUUM");
}

void Logger::SetOutMain(const char* filename)
{
    if (!fout_main_.is_open()) {
        fout_main_.open(filename, std::ofstream::app);
        if (!fout_main_.is_open()) {
            fmt::print("Cannot open {}\n", filename);
            throw MyException(0xc396a312U, "Cannot open file.");
        }
    }
}

void Logger::SetOutDeduce(const char* filename)
{

    if (!db_deduce_.is_open()) {
        db_deduce_.open(filename);
        db_deduce_.create_log();
        db_deduce_.InsertTag(-1, cmd_);
        db_deduce_.begin_transaction();

        cmd_.clear();
    }
}

std::string Logger::GetCmd(int argc, char** argv)
{
    return fmt::format("cmd: {}\nLogging start at {}\n", fmt::join(argv, argv + argc, " "), ut::get_time());
}

void Logger::LogCmd(int argc, char** argv)
{
    cmd_ = GetCmd(argc, argv);
    fmt::print(fout_main_, "{}", cmd_);
}

void Logger::LogSummary(const std::string& category, int count)
{
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "Summary - {}: {}          \n", category, count);
    fmt::print("{}", line_);
    fmt::print(fout_main_, "{}", line_);
}

void Logger::LogTime(const std::string& time)
{
    fmt::print(fmt::fg(fmt::color::green), "{}\n", time);
    fmt::print(fout_main_, "{}\n", time);
    db_deduce_.InsertTag(-2, time);
}

void Logger::Checkpoint()
{
    id_checkpoint_ = db_deduce_.get_int("SELECT max(id) FROM log");
}

void Logger::RollBackToCheckpoint() {
    db_deduce_.execute_cmd(fmt::format("DELETE FROM log WHERE id>{}", id_checkpoint_));
}

void Logger::LogSSException(int depth, const std::string& name, alg::AdamsDeg deg_dx, const alg::int1d& dx, int r, unsigned code, alg::AdamsDeg deg_leibniz, const alg::int1d* a_leibniz)
{
    if (depth == 0)
        fmt::print("Error({:#x}) No source for the image. {} deg_dx={}, dx={}, r={}\n", code, name, deg_dx, dx, r);
    std::string tag;
    if (a_leibniz) {
        tag = fmt::format("{:#x} {} {}", code, deg_leibniz, *a_leibniz);
        a_leibniz = nullptr;
    }
    else
        tag = fmt::format("{:#x}", code);
    db_deduce_.InsertError(depth, name, deg_dx, dx, r, tag);
}

void Logger::LogSSSSException(int depth, unsigned code)
{
    std::string tag = fmt::format("{:#x}", code);
    db_deduce_.InsertError(depth, tag);
}

void Logger::LogDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0)
        fmt::print(fmt::fg(dx.empty() ? fmt::color::white_smoke : fmt::color::light_green), "{}{} - {} {} d_{}{}={}\n", indent, GetReason(reason), name, deg_x, r, x, dx);
    db_deduce_.InsertDiff(depth, reason, name, deg_x, x, dx, r);
}

void Logger::LogDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0)
        fmt::print(fmt::fg(fmt::color::gray), "{}{} {} d_{}{}={}\n", indent, name, deg_x, r, x, dx);
    db_deduce_.InsertNullDiff(depth, name, deg_x, x, r);
}

void Logger::LogNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, int r)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0)
        fmt::print(fmt::fg(fmt::color::gray), "{}{} {} d_{}{}=?\n", indent, name, deg_x, r, x);
    db_deduce_.InsertNullDiff(depth, name, deg_x, x, r);
}

void Logger::LogDiffInv(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, alg::AdamsDeg deg_dx, const alg::int1d& x, const alg::int1d& dx, int r)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0)
        fmt::print(fmt::fg(x.empty() ? fmt::color::white_smoke : fmt::color::light_green), "{}{} - {} {} {}=d_{}{}\n", indent, GetReason(reason), name, deg_dx, dx, r, x);
    db_deduce_.InsertDiff(depth, reason, name, deg_x, x, dx, r);
}

//void Logger::LogDiffBoun(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_dx, const alg::int1d& dx)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {} is a boundary\n", indent, GetReason(reason), name, deg_dx, dx);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::light_green), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyGen(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, size_t gen_id, const alg2::Mod& Einf)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} x_{} detected by {}\n", indent, GetReason(reason), name, deg, gen_id, Einf);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyGen(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, size_t gen_id, const alg2::Poly& Einf)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} x_{} detected by {}\n", indent, GetReason(reason), name, deg, gen_id, Einf);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Poly& rel)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}=0\n", indent, GetReason(reason), name, deg_rel, rel);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::white_smoke), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Mod& rel)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}=0\n", indent, GetReason(reason), name, deg_rel, rel);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::white_smoke), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Poly& rel1, const algZ::Poly& rel2)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {} --> {}\n", indent, GetReason(reason), name, deg_rel, rel1, rel2);
    if (depth == 0) {
        fmt::print(fmt::fg(rel1.data.size() == rel2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Mod& rel1, const algZ::Mod& rel2)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {} --> {}\n", indent, GetReason(reason), name, deg_rel, rel1, rel2);
    if (depth == 0) {
        fmt::print(fmt::fg(rel1.data.size() == rel2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyMap(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const std::string& f, size_t gen_id, const algZ::Poly& fx)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}(v_{})={}\n", indent, GetReason(reason), name, deg_x, f, gen_id, fx);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::white_smoke), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyMap(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const std::string& f, size_t gen_id, const algZ::Poly& fx1, const algZ::Poly& fx2)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}(v_{}) = {} --> {}\n", indent, GetReason(reason), name, deg_x, f, gen_id, fx1, fx2);
    if (depth == 0) {
        fmt::print(fmt::fg(fx1.data.size() == fx2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyProd(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Poly& m, const algZ::Poly& hm1, const algZ::Poly& hm2)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}*{} = {} --> {}\n", indent, GetReason(reason), name, deg_x, h, m, hm1, hm2);
    if (depth == 0) {
        fmt::print(fmt::fg(hm1.data.size() == hm2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}

//void Logger::LogHtpyProd(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Mod& m, const algZ::Mod& hm1, const algZ::Mod& hm2)
//{
    /*std::string_view indent(INDENT, depth * 2);
    line_.clear();
    fmt::format_to(std::back_inserter(line_), "{}{} - {} {} {}*{} = {} --> {}\n", indent, GetReason(reason), name, deg_x, h, m, hm1, hm2);
    if (depth == 0) {
        fmt::print(fmt::fg(hm1.data.size() == hm2.data.size() ? fmt::color::white_smoke : fmt::color::light_blue), "{}", line_);
        fmt::print(fout_deduce_, "{}", line_);
    }
    else
        fmt::format_to(std::back_inserter(out_), "{}", line_);*/
//}
