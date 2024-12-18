#include "mylog.h"
#include "main.h"
#include <fmt/format.h>
#include <fmt/ranges.h>

std::ofstream Logger::fout_main_;
std::string Logger::cmd_;
std::string Logger::line_;
DbLog Logger::db_;

std::string_view GetReason(EnumReason reason)
{
    return REASONS[size_t(reason)];
}

std::string_view GetReasonDb(EnumReason reason)
{
    return REASONS_DB[size_t(reason)];
}

/* -1: start time
 * -2: end time
 * -3:
 */
void DbLog::InsertInfo(int depth, const std::string& info)
{
    Statement stmt(*this, "INSERT INTO log (depth, info) VALUES (?1, ?2);");
    stmt.bind_and_step(depth, std::string(info));
}

void DbLog::InsertExclusion(int id)
{
    Statement stmt(*this, "INSERT INTO exclusions (id) VALUES (?1);");
    stmt.bind_and_step(id);
}

void DbLog::InsertDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, int r, const alg::int1d& x, const alg::int1d& dx, const std::string& info)
{
    if (info.empty()) {
        Statement stmt(*this, "INSERT INTO log (depth, reason, name, s, t, r, x, dx) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8)");
        stmt.bind_and_step(depth, std::string(GetReasonDb(reason)), name, deg.s, deg.t, r, myio::Serialize(x), myio::Serialize(dx));
    }
    else {
        Statement stmt(*this, "INSERT INTO log (depth, reason, name, s, t, r, x, dx, info) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)");
        stmt.bind_and_step(depth, std::string(GetReasonDb(reason)), name, deg.s, deg.t, r, myio::Serialize(x), myio::Serialize(dx), info);
    }
}

void DbLog::InsertDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, const alg::int1d& dx, int r)
{
    Statement stmt(*this, "INSERT INTO log (depth, name, s, t, r, x, dx) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)");
    stmt.bind_and_step(depth, name, deg_x.s, deg_x.t, r, myio::Serialize(x), myio::Serialize(dx));
}

void DbLog::InsertNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, const alg::int1d& x, int r)
{
    Statement stmt(*this, "INSERT INTO log (depth, name, s, t, r, x) VALUES (?1, ?2, ?3, ?4, ?5, ?6)");
    stmt.bind_and_step(depth, name, deg_x.s, deg_x.t, r, myio::Serialize(x));
}

void Logger::Reset()
{
    db_.drop_table("log");
    db_.drop_table("exclusions");
    db_.drop_table("minimal");
    db_.create_log();
    db_.end_transaction();
    db_.execute_cmd("VACUUM");
    db_.begin_transaction();
}

void Logger::SetOutMain(const char* filename)
{
    if (!fout_main_.is_open()) {
        fout_main_.open(filename, std::ofstream::app);
        if (!fout_main_.is_open()) {
            fmt::print("Cannot open {}\n", filename);
            throw ErrorIdMsg(0xc396a312U, "Cannot open file.");
        }
    }
}

void Logger::SetOutDeduce(const char* filename)
{

    if (!db_) {
        db_.open(filename);
        db_.create_log();
        db_.InsertInfo(-1, cmd_);
        db_.begin_transaction();

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
    if (db_)
        db_.InsertInfo(-2, time);
    else
        fmt::print("Warning: No database for logging.\n");
}

int Logger::GetCheckpoint()
{
    if (db_)
        return db_.get_int("SELECT COALESCE(MAX(id), -1) FROM log");
    return 0;
}

void Logger::RollBackToCheckpoint(int checkpt)
{
    if (db_)
        db_.execute_cmd(fmt::format("DELETE FROM log WHERE id>{}", checkpt));
}

void Logger::RollBackToExclusionsCheckpoint(int checkpt)
{
    if (db_)
        db_.execute_cmd(fmt::format("DELETE FROM exclusions WHERE id>{}", checkpt));
}

void Logger::LogExclusion(int id)
{
    if (db_)
        db_.InsertExclusion(id);
}

void Logger::LogDiff(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, int r, const alg::int1d& x, const alg::int1d& dx, const std::string& proof, SSFlag)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0)
        fmt::print(fmt::fg(dx.empty() ? fmt::color::white_smoke : fmt::color::light_green), "{}{} t={} - {} {} d_{}{}={}\n", indent, GetReason(reason), deg_x.t, name, deg_x, r, x, dx);

    if (depth > 0) { /* For depth>0 only log try */
        switch (reason) {
        case EnumReason::try1:
        case EnumReason::try2:
        case EnumReason::deduce:
        case EnumReason::deduce2:
            break;
        default:
            return;
        }
    }
    /*switch (reason) {
    case EnumReason::nat:
        if (!(flag & SSFlag::log_nat))
            return;
    case EnumReason::degree:
    case EnumReason::degree2:
    case EnumReason::deduce_xx:
        if (!(flag & SSFlag::log_deg))
            return;
    }*/
    if (db_)
        db_.InsertDiff(depth, reason, name, deg_x, r, x, dx, proof);
}

void Logger::LogDiffInv(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_dx, int r, const alg::int1d& x, const alg::int1d& dx, const std::string& proof, SSFlag)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0)
        fmt::print(fmt::fg(x.empty() ? fmt::color::white_smoke : fmt::color::light_green), "{}{} t={} - {} {} {}=d_{}{}\n", indent, GetReason(reason), deg_dx.t, name, deg_dx, dx, r, x);

    if (depth > 0) { /* For depth>0 only log try */
        switch (reason) {
        case EnumReason::try1:
        case EnumReason::try2:
            break;
        default:
            return;
        }
    }
   /* switch (reason) {
    case EnumReason::nat:
        if (!(flag & SSFlag::log_nat))
            return;
    case EnumReason::degree:
    case EnumReason::degree2:
        if (!(flag & SSFlag::log_deg))
            return;
    }*/
    if (db_)
        db_.InsertDiff(depth, reason, name, deg_dx, r, x, dx, proof);
}

void Logger::LogGreyDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, int r, const alg::int1d& x, const alg::int1d& dx, SSFlag)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::gray), "{}{} {} d_{}{}={}\n", indent, name, deg_x, r, x, dx);
        if (db_)
            db_.InsertDiff(depth, name, deg_x, x, dx, r); /* A hint of what is used */
    }
}

void Logger::LogNullDiff(int depth, const std::string& name, alg::AdamsDeg deg_x, int r, const alg::int1d& x, SSFlag)
{
    std::string_view indent(INDENT, depth * 2);
    if (depth == 0) {
        fmt::print(fmt::fg(fmt::color::gray), "{}{} {} d_{}{}=?\n", indent, name, deg_x, r, x);
        if (db_)
            db_.InsertNullDiff(depth, name, deg_x, x, r);
    }
}

// void Logger::LogDiffBoun(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_dx, const alg::int1d& dx)
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

// void Logger::LogHtpyGen(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, size_t gen_id, const alg2::Mod& Einf)
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

// void Logger::LogHtpyGen(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg, size_t gen_id, const alg2::Poly& Einf)
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

// void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Poly& rel)
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

// void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Mod& rel)
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

// void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Poly& rel1, const algZ::Poly& rel2)
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

// void Logger::LogHtpyRel(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_rel, const algZ::Mod& rel1, const algZ::Mod& rel2)
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

// void Logger::LogHtpyMap(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const std::string& f, size_t gen_id, const algZ::Poly& fx)
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

// void Logger::LogHtpyMap(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const std::string& f, size_t gen_id, const algZ::Poly& fx1, const algZ::Poly& fx2)
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

// void Logger::LogHtpyProd(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Poly& m, const algZ::Poly& hm1, const algZ::Poly& hm2)
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

// void Logger::LogHtpyProd(int depth, EnumReason reason, const std::string& name, alg::AdamsDeg deg_x, const algZ::Poly& h, const algZ::Mod& m, const algZ::Mod& hm1, const algZ::Mod& hm2)
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
