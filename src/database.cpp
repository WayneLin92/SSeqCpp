#include "database.h"
#include <fmt/format.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sqlite3.h>
#include <sstream>

static_assert(SQLITE_ROW == MYSQLITE_ROW);
static_assert(SQLITE3_TEXT == MYSQLITE3_TEXT);

namespace myio {

template <>
int1d Deserialize<int1d>(const std::string& str)
{
    int1d result;
    if (str.empty())
        return result;
    std::stringstream ss(str);
    while (ss.good()) {
        int i;
        ss >> i;
        result.push_back(i);
        if (ss.peek() == ',')
            ss.ignore();
    }
    return result;
}

Database::Database(std::string filename) : filename_(std::move(filename))
{
    if (myio::FileExists(filename_))
        newFile_ = false;
    if (sqlite3_open(filename_.c_str(), &conn_) != SQLITE_OK) {
        fmt::print("Cannot open database {}\n", filename_);
        throw ErrorIdMsg(0x8de81e80, "Cannot open database");
    }
}

Database::~Database()
{
    end_transaction();
    if (conn_) {
        sqlite3_close(conn_);
        conn_ = nullptr;
    }
}

void myio::Database::open(std::string filename)
{
    filename_ = std::move(filename);
    if (myio::FileExists(filename_))
        newFile_ = false;
    if (sqlite3_open(filename_.c_str(), &conn_) != SQLITE_OK)
        throw RunTimeError(fmt::format("Cannot open database: {}.", filename_));
}

void Database::disconnect()
{
    end_transaction();
    if (conn_) {
        sqlite3_close(conn_);
        conn_ = nullptr;
    }
}

void Database::sqlite3_prepare(const char* zSql, sqlite3_stmt** ppStmt) const
{
    int error_code = sqlite3_prepare_v2(conn_, zSql, int(strlen(zSql)) + 1, ppStmt, NULL);
    if (error_code != SQLITE_OK)
        throw RunTimeError(fmt::format("Sqlite3 failed to compile.\nDatabase: {}\nSQL: {}\nError: {}", filename_, zSql, sqlite3_errstr(error_code)));
}

void Database::sqlite3_prepare(const std::string& sql, sqlite3_stmt** ppStmt) const
{
    int error_code = sqlite3_prepare_v2(conn_, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL);
    if (error_code != SQLITE_OK)
        throw RunTimeError(fmt::format("Sqlite3 failed to compile.\nDatabase: {}\nSQL: {}\nError: {}", filename_, sql, sqlite3_errstr(error_code)));
}

void Database::execute_cmd(const std::string& sql) const
{
    char* zErrMsg = nullptr;
    int rc = sqlite3_exec(conn_, sql.c_str(), nullptr, nullptr, &zErrMsg);
    if (rc != SQLITE_OK) {
        std::string msg;
        if (zErrMsg) {
            msg = zErrMsg;
            sqlite3_free(zErrMsg);
        }
        throw RunTimeError(fmt::format("Sqlite3 failed to execute.\nDatabase: {}\nSQL: {}\nError: {}", filename_, sql, msg));
    }
}

int Database::get_int(const std::string& sql) const
{
    Statement stmt(*this, sql);
    if (stmt.step() == SQLITE_ROW)
        if (stmt.column_type(0) == SQLITE_INTEGER)
            return stmt.column_int(0);
    throw RunTimeError(sql);
}

int Database::get_int(const std::string& sql, int default_) const
{
    Statement stmt(*this, sql);
    if (stmt.step() == SQLITE_ROW) {
        if (stmt.column_type(0) == SQLITE_INTEGER)
            return stmt.column_int(0);
        else {
            throw RunTimeError(fmt::format("{}: Incorrect type", sql));
        }
    }
    return default_;
}

int Database::get_metadata_int(std::string_view key) const
{
    try {
        if (has_table("version")) {
            if (get_int(fmt::format("select count(*) from version where name=\"{}\"", key)) > 0)
                return get_int(fmt::format("select value from version where name=\"{}\"", key));
            else
                return -1000;
        }
    }
    catch (MyException&) { /* has_table("version") could throw an error when the database is locked or corupted */
        return -3000;
    }
    return -2000;
}

std::string Database::get_str(const std::string& sql) const
{
    Statement stmt(*this, sql);
    if (stmt.step() == SQLITE_ROW)
        if (stmt.column_type(0) == SQLITE_TEXT)
            return stmt.column_str(0);
    throw RunTimeError(sql);
}

std::string Database::get_metadata_str(std::string_view key) const
{
    try {
        if (has_table("version")) {
            if (get_int(fmt::format("select count(*) from version where name=\"{}\"", key)) > 0)
                return get_str(fmt::format("select value from version where name=\"{}\"", key));
            else
                return "Error: key not found";
        }
    }
    catch (MyException&) { /* has_table("version") could throw an error when the database is locked or corupted */
        return "Error: database is locked or corrupted";
    }
    return "Error: database has no version table";
}

std::vector<int> Database::get_column_int(const std::string& table_name, const std::string& column_name, const std::string& conditions) const
{
    std::vector<int> result;
    Statement stmt(*this, "SELECT " + column_name + " FROM " + table_name + ' ' + conditions + ';');
    while (stmt.step() == SQLITE_ROW)
        result.push_back(stmt.column_int(0));
    return result;
}

Statement::~Statement()
{
    sqlite3_finalize(stmt_);
}
Statement::Statement(const Database& db, const std::string& sql)
{
    db.sqlite3_prepare(sql, &stmt_);
    sql_ = sql;
}

void Statement::bind_int(int iCol, int i) const
{
    if (sqlite3_bind_int(stmt_, iCol, i) != SQLITE_OK)
        throw ErrorIdMsg(0xa61e05b2, "Sqlite bind_int fail"); // TODO: Use RunTimeError instead of ErrorIdMsg
}

void Statement::bind_int64(int iCol, int64_t i) const
{
    if (sqlite3_bind_int64(stmt_, iCol, i) != SQLITE_OK)
        throw ErrorIdMsg(0xb78985c, "Sqlite bind_int64 fail");
}

void Statement::bind_double(int iCol, double d) const
{
    if (sqlite3_bind_double(stmt_, iCol, d) != SQLITE_OK)
        throw ErrorIdMsg(0x6643edffU, "Sqlite bind_double fail");
}

void Statement::bind_null(int iCol) const
{
    if (sqlite3_bind_null(stmt_, iCol) != SQLITE_OK)
        throw ErrorIdMsg(0xe22b11c4, "Sqlite bind_null fail");
}

/* Warning: SQLITE_STATIC is used here. str should not change before step_and_reset */
void Statement::bind_str(int iCol, const std::string& str) const
{
    if (sqlite3_bind_text(stmt_, iCol, str.c_str(), -1, SQLITE_STATIC) != SQLITE_OK)
        throw ErrorIdMsg(0x29cc3b21, "Sqlite bind_str fail");
}

/* Warning: SQLITE_STATIC is used here. str should not change before step_and_reset */
void Statement::bind_blob(int iCol, const void* data, int nBytes) const
{
    if (sqlite3_bind_blob(stmt_, iCol, data, nBytes, SQLITE_STATIC) != SQLITE_OK)
        throw ErrorIdMsg(0x79d80302, "Sqlite bind_blob fail");
}

std::string Statement::column_str(int iCol) const
{
    return std::string(reinterpret_cast<const char*>(sqlite3_column_text(stmt_, iCol)));
}

int Statement::column_int(int iCol) const
{
    return sqlite3_column_int(stmt_, iCol);
}

int Statement::column_type(int iCol) const
{
    return sqlite3_column_type(stmt_, iCol);
}
const void* Statement::column_blob(int iCol) const
{
    return sqlite3_column_blob(stmt_, iCol);
}
int Statement::column_blob_size(int iCol) const
{
    return sqlite3_column_bytes(stmt_, iCol);
}

void Statement::step_and_reset() const
{
    if (int rc = step(); rc != SQLITE_DONE)
        throw RunTimeError(fmt::format("Sqlite step fail: {}. SQL: {}", sqlite3_errstr(rc), sql_));
    sqlite3_reset(stmt_);
}

int Statement::step() const
{
    return sqlite3_step(stmt_);
}

void Statement::reset() const
{
    if (int rc = sqlite3_reset(stmt_); rc != SQLITE_OK)
        throw RunTimeError(fmt::format("Sqlite reset fail: {}. SQL: {}", sqlite3_errstr(rc), sql_));
}

}  // namespace myio
