#include "database.h"
#include <sqlite3.h>
#include <cstring>
#include <iostream>
#include <sstream>

#if SQLITE_ROW != MYSQLITE_ROW
#error "MYSQLITE_ROW is not equal to SQLITE_ROW"
#endif
#if SQLITE3_TEXT != MYSQLITE3_TEXT
#error "MYSQLITE3_TEXT is not equal to SQLITE3_TEXT"
#endif

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

Database::Database(const std::string& filename)
{
    if (sqlite3_open(filename.c_str(), &conn_) != SQLITE_OK)
        throw MyException(0x8de81e80, std::string("database ") + filename + " not found");
}

Database::~Database()
{
    end_transaction();
    sqlite3_close(conn_);
}

void Database::sqlite3_prepare(const char* zSql, sqlite3_stmt** ppStmt) const
{
    int error_code = sqlite3_prepare_v2(conn_, zSql, int(strlen(zSql)) + 1, ppStmt, NULL);
    if (error_code != SQLITE_OK) {
        std::cout << "0xbce2dcfeU: Sqlite3 compiling " + std::string(zSql) + " :" + sqlite3_errstr(error_code) << '\n';
        throw MyException(0xbce2dcfeU, std::string("Sqlite3 compiling ") + zSql + " :" + sqlite3_errstr(error_code));
    }
}

void Database::sqlite3_prepare(const std::string& sql, sqlite3_stmt** ppStmt) const
{
    int error_code = sqlite3_prepare_v2(conn_, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL);
    if (error_code != SQLITE_OK) {
        std::cout << "0xda6ab7f6U: Sqlite3 compiling " + sql + " :" + sqlite3_errstr(error_code) << '\n';
        throw MyException(0xda6ab7f6U, "Sqlite3 compiling " + sql + " :" + sqlite3_errstr(error_code));
    }
}

void Database::execute_cmd(const std::string& sql) const
{
    sqlite3_stmt* stmt;
    sqlite3_prepare(sql, &stmt);
    sqlite3_step(stmt);
    sqlite3_finalize(stmt);
}

int Database::get_int(const std::string& sql) const
{
    Statement stmt(*this, sql);
    if (stmt.step() == SQLITE_ROW)
        if (stmt.column_type(0) == SQLITE_INTEGER)
            return stmt.column_int(0);
    throw MyException(0xeffbf28c, "Failed to get_int using: " + sql);
}

int Database::get_int(const std::string& sql, int default_) const
{
    Statement stmt(*this, sql);
    if (stmt.step() == SQLITE_ROW) {
        if (stmt.column_type(0) == SQLITE_INTEGER)
            return stmt.column_int(0);
        else
            throw MyException(0xa429d1aU, "Incorrect type using: " + sql);
    }
    return default_;
}

std::string Database::get_str(const std::string& sql) const
{
    Statement stmt(*this, sql);
    if (stmt.step() == SQLITE_ROW)
        if (stmt.column_type(0) == SQLITE_TEXT)
            return stmt.column_str(0);
    throw MyException(0xeffbf28c, "Failed to get_str using: " + sql); 
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
}

void Statement::bind_str(int iCol, const std::string& str) const
{
    if (sqlite3_bind_text(stmt_, iCol, str.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) /* SQLITE_TRANSIENT is important here */
        throw MyException(0x29cc3b21, "Sqlite bind_str fail");
}

void Statement::bind_int(int iCol, int i) const
{
    if (sqlite3_bind_int(stmt_, iCol, i) != SQLITE_OK)
        throw MyException(0xa61e05b2, "Sqlite bind_int fail");
}

void Statement::bind_int64(int iCol, int64_t i) const
{
    if (sqlite3_bind_int64(stmt_, iCol, i) != SQLITE_OK)
        throw MyException(0xb78985c, "Sqlite bind_int64 fail");
}

void Statement::bind_double(int iCol, double d) const
{
    if (sqlite3_bind_double(stmt_, iCol, d) != SQLITE_OK)
        throw MyException(0x6643edffU, "Sqlite bind_double fail");
}

void Statement::bind_null(int iCol) const
{
    if (sqlite3_bind_null(stmt_, iCol) != SQLITE_OK)
        throw MyException(0xe22b11c4, "Sqlite bind_null fail");
}

void Statement::bind_blob(int iCol, const void* data, int nBytes) const
{
    if (sqlite3_bind_blob(stmt_, iCol, data, nBytes, SQLITE_TRANSIENT) != SQLITE_OK)
        throw MyException(0x79d80302, "Sqlite bind_blob fail");
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

int Statement::step() const
{
    return sqlite3_step(stmt_);
}

int Statement::reset() const
{
    return sqlite3_reset(stmt_);
}

void Statement::step_and_reset() const
{
    step();
    reset();
}
}  // namespace myio

