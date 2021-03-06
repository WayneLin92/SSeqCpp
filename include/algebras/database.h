/** \file database.h
 * A module for interacting with a database.
 */
#ifndef DATABASE_H
#define DATABASE_H

#include "myexception.h"
#include "myio.h"
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

struct sqlite3;
struct sqlite3_stmt;
#define MYSQLITE_ROW 100
#define MYSQLITE3_TEXT 3

namespace ut {

template <typename>
struct is_vector : std::false_type
{
};

template <typename T>
struct is_vector<std::vector<T>> : std::true_type
{
    using value_type = T;
};

}  // namespace ut

/**
 * This namespace provides C++ wrappers for the sqlite3 library.
 */
namespace myio {

using int1d = std::vector<int>;

inline std::string Serialize(const int1d& arr)
{
    return StrCont("", ",", "", "", arr, [](int i) { return std::to_string(i); });
}

template <typename T>
T Deserialize(const std::string& str)
{
    throw MyException(0xb5c7695cU, "Must use a specialization");
}

template <>
int1d Deserialize<int1d>(const std::string& str);

class Database;

/**
 * Wrapper for `sqlite3_stmt*`
 */
class Statement
{
public:
    Statement() = default;
    Statement(const Database& db, const std::string& sql);
    ~Statement();

private:
    sqlite3_stmt* stmt_ = nullptr;

public:
    void bind_str(int iCol, const std::string& str) const;
    void bind_int(int iCol, int i) const;
    void bind_int64(int iCol, int64_t i) const;
    void bind_double(int iCol, double d) const;
    void bind_null(int iCol) const;
    void bind_blob(int iCol, const void* data, int nBytes) const;
    template <typename T>
    void bind_blob(int iCol, const std::vector<T>& data) const
    {
        static_assert(!ut::is_vector<T>::value, "T must not be a vector.");  // TODO: complete the bind_blob for vector of vectors
        if (data.data())
            bind_blob(iCol, data.data(), int(data.size() * sizeof(T)));
        else
            bind_blob(iCol, this, 0);
    }
    std::string column_str(int iCol) const;
    int column_int(int iCol) const;
    int column_type(int iCol) const;
    const void* column_blob(int iCol) const;
    int column_blob_size(int iCol) const;
    template <typename T>
    std::vector<T> column_blob_tpl(int iCol) const
    {
        std::vector<T> result;
        const void* data = column_blob(iCol);
        int bytes = column_blob_size(iCol);
        size_t data_size = (size_t)bytes / sizeof(T);
        result.resize(data_size);
        std::memcpy(result.data(), data, bytes);
        return result;
    }
    int step() const;
    int reset() const;
    void step_and_reset() const;
};

/**
 * Wrapper for `sqlite*`
 */
class Database
{
private:
    sqlite3* conn_ = nullptr;
    int numInTransaction_ = -1;

public:
    Database() = default;
    explicit Database(const std::string& filename);
    ~Database();

public:
    void sqlite3_prepare(const char* zSql, sqlite3_stmt** ppStmt) const;
    void sqlite3_prepare(const std::string& sql, sqlite3_stmt** ppStmt) const;
    void execute_cmd(const std::string& sql) const;
    void begin_transaction()
    {
        if (numInTransaction_ == -1) {
            execute_cmd("BEGIN TRANSACTION");
            numInTransaction_ = 0;
        }
    }
    void end_transaction(int min=0)
    {
        if (numInTransaction_ >= min) {
            execute_cmd("END TRANSACTION");
            numInTransaction_ = -1;
        }
    }
    void reg_transaction()
    {
        if (numInTransaction_ >= 0)
            ++numInTransaction_;
    }

public:
    void delete_from(const std::string& table_name) const
    {
        execute_cmd("DELETE FROM " + table_name);
    }

public:
    int get_int(const std::string& sql) const;
    std::string get_str(const std::string& sql) const;
    std::vector<int> get_column_int(const std::string& table_name, const std::string& column_name, const std::string& conditions) const;
    /*
     * This converts a column of strings to vector of type `T`.
     */
    template <typename T, typename FnMap>
    std::vector<T> get_column_from_str(const std::string& table_name, const std::string& column_name, const std::string& conditions, const FnMap& map) const
    {
        std::vector<T> result;
        Statement stmt(*this, "SELECT " + column_name + " FROM " + table_name + ' ' + conditions + ';');
        while (stmt.step() == MYSQLITE_ROW)
            result.push_back(map(stmt.column_str(0)));
        std::clog << column_name << "'s loaded from " << table_name << ", size=" << result.size() << '\n';
        return result;
    }
    /* `map` takes two arguments (void*, int)->T */
    template <typename T, typename FnMap>
    std::vector<T> get_column_from_blob(const std::string& table_name, const std::string& column_name, const std::string& conditions, const FnMap& map) const
    {
        std::vector<T> result;
        Statement stmt(*this, "SELECT " + column_name + " FROM " + table_name + ' ' + conditions + ';');
        while (stmt.step() == MYSQLITE_ROW)
            result.push_back(map(stmt.column_blob(0), stmt.column_blob_size(0)));
        std::clog << column_name << "'s loaded from " << table_name << ", size=" << result.size() << '\n';
        return result;
    }
    template <typename T, typename FnMap>
    std::vector<T> get_column_from_str_with_null(const std::string& table_name, const std::string& column_name, const T& null_value, const std::string& conditions, const FnMap& map) const
    {
        std::vector<T> result;
        Statement stmt(*this, "SELECT " + column_name + " FROM " + table_name + ' ' + conditions + ';');
        while (stmt.step() == MYSQLITE_ROW)
            if (stmt.column_type(0) == MYSQLITE3_TEXT)
                result.push_back(map(stmt.column_str(0)));
            else
                result.push_back(null_value);
        std::clog << column_name << "'s loaded from " << table_name << ", size=" << result.size() << '\n';
        return result;
    }
    std::vector<std::string> get_column_str(const std::string& table_name, const std::string& column_name, const std::string& conditions) const
    {
        return get_column_from_str<std::string>(table_name, column_name, conditions, [](std::string c) { return c; });
    }

public:
    template <typename T, typename FnMap>
    void update_str_column(const std::string& table_name, const std::string& column_name, const std::string& index_name, const std::vector<T>& column, const FnMap& map, size_t i_start) const  ///
    {
        Statement stmt(*this, "UPDATE " + table_name + " SET " + column_name + " = ?1 WHERE " + index_name + "= ?2;");
        for (size_t i = i_start; i < column.size(); ++i) {
            stmt.bind_str(1, map(column[i]));
            stmt.bind_int(2, (int)i);
            stmt.step_and_reset();
        }
    }
    /* `map` should return a pair of type (void*, int) */
    template <typename T, typename FnMap>
    void update_blob_column(const std::string& table_name, const std::string& column_name, const std::string& index_name, const std::vector<T>& column, const FnMap& map, size_t i_start) const
    {
        Statement stmt(*this, "UPDATE " + table_name + " SET " + column_name + " = ?1 WHERE " + index_name + "= ?2;");
        for (size_t i = i_start; i < column.size(); ++i) {
            stmt.bind_blob(1, map(column[i]).first, (int)map(column[i]).second);
            stmt.bind_int(2, (int)i);
            stmt.step_and_reset();
        }
    }
};
}  // namespace myio

#endif /* DATABASE_H */
