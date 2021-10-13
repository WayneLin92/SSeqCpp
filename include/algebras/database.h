/** \file database.h
 * A component for interaction with a database.
 * TODO: Template member functions.
 * TODO: Change load_gb/save_gb. Add version control system.
 */
#ifndef DATABASE_H
#define DATABASE_H

#define DATABASE_SAVE_LOGGING /* This is a switch to print what are saved to the database */

#include "algebras.h"
#include <map>
#include <string>
#include <vector>

struct sqlite3;
struct sqlite3_stmt;

/**
 * This is namespace provides C++ wrappers for the sqlite3 library.
 */
namespace myio {

enum struct SqlType
{
    SqlInt = 0,
    SqlStr = 1,
};

/**
 * Wrapper for `sqlite*`
 */
class Database
{
public:
    Database() = default;
    explicit Database(const char* filename);
    ~Database();

private:
    sqlite3* conn_ = nullptr;

public:
    void sqlite3_prepare(const char* zSql, sqlite3_stmt** ppStmt, bool bPrintError = false) const;
    void sqlite3_prepare(const std::string& sql, sqlite3_stmt** ppStmt, bool bPrintError = false) const;
    void execute_cmd(const std::string& sql) const;
    void begin_transaction() const
    {
        execute_cmd("BEGIN TRANSACTION");
    }
    void end_transaction() const
    {
        execute_cmd("END TRANSACTION");
    }

public:
    int get_int(const std::string& sql) const;
    std::vector<int> get_column_int(const std::string& table_name, const std::string& column_name, const std::string& conditions) const;
};

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
    void bind_null(int iCol) const;
    const char* column_str(int iCol) const;
    int column_int(int iCol) const;
    int column_type(int iCol) const;
    int step() const;
    int reset() const;
    void step_and_reset() const;
};
}  // namespace myio

namespace alg {

constexpr auto kLevelMax = 10000;

struct BasisComplex
{
    alg::array2d boundaries;
    alg::array2d cycles;
};

struct Staircase
{
    alg::array2d basis_ind;
    alg::array2d diffs_ind; /* ={-1} means null */
    alg::array levels;
};

class Database : public myio::Database
{
public:
    Database() = default;
    explicit Database(const char* filename) : myio::Database(filename) {}

public:
    std::vector<array> get_column_array(const std::string& table_name, const std::string& column_name, const std::string& conditions) const;
    Poly1d get_column_Poly(const std::string& table_name, const std::string& column_name, const std::string& conditions) const;

public:
    std::vector<Deg> load_gen_degs(const std::string& table_name) const;
    std::vector<std::string> load_gen_names(const std::string& table_name) const;
    Poly1d load_gen_diffs(const std::string& table_name) const;
    Poly1d load_gen_reprs(const std::string& table_name) const;
    Poly1d load_gen_images(const std::string& table_name, const std::string& column_name, int t_max) const;
    Mon2d load_leading_terms(const std::string& table_name, int t_max) const;
    Poly1d load_gb(const std::string& table_name, int t_max) const;
    Poly1d load_gb_s(const std::string& table_name, int s_max) const;
    std::map<Deg, int> load_indices(const std::string& table_name, int t_max) const;
    std::map<Deg, Mon1d> load_basis(const std::string& table_name, int t_max) const;
    std::map<Deg, array2d> load_mon_diffs_ind(const std::string& table_name, int t_max, int withnull = false) const;
    std::map<Deg, Poly1d> load_mon_diffs(const std::string& table_name, const std::map<Deg, Mon1d>& basis, int r, int t_max) const;
    std::map<Deg, BasisComplex> load_basis_ss(const std::string& table_name_ss, int r, int t_max) const; /* load E_r-cycles and E_r-boundaries */
    std::map<Deg, Staircase> load_basis_ss(const std::string& table_name_ss, int t_max) const;

public:
    void save_generators(const std::string& table_name, const std::vector<std::string>& gen_names, const std::vector<Deg>& gen_degs) const;
    void save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs, size_t i_start) const;
    void save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs) const
    {
        save_generators(table_name, gen_degs, gen_reprs, 0);
    }
    void save_generators(const std::string& table_name, const std::vector<std::string>& gen_names, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs) const;
    void save_gen_images(const std::string& table_name, const std::string& column_name, const Poly1d& gen_images) const;
    void save_gb(const std::string& table_name, const Poly1d& gb, const std::vector<Deg>& gen_degs, size_t i_start) const;
    void save_gb(const std::string& table_name, const Poly1d& gb, const std::vector<Deg>& gen_degs) const
    {
        save_gb(table_name, gb, gen_degs, 0);
    }
    void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis) const;
    void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, array2d>& mon_reprs) const;
    void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, Poly1d>& mon_reprs) const;
    void save_basis_ss(const std::string& table_name, const std::map<Deg, Staircase>& basis_ss) const;
    void update_basis_ss(const std::string& table_name, const std::map<Deg, Staircase>& basis_ss) const;
};

array DbStr2Array(const char* str_mon);
Mon DbStr2Mon(const char* str_mon);
Poly DbStr2Poly(const char* str_poly);
std::string Array2DbStr(array::const_iterator pbegin, array::const_iterator pend);
inline std::string Array2DbStr(const array& a)
{
    return Array2DbStr(a.begin(), a.end());
}
std::string Mon2DbStr(MonInd pbegin, MonInd pend);
inline std::string Mon2DbStr(const Mon& mon)
{
    return Mon2DbStr(mon.begin(), mon.end());
}
std::string Poly_to_str(Poly::const_iterator pbegin, Poly::const_iterator pend);
inline std::string Poly_to_str(const Poly& poly)
{
    return Poly_to_str(poly.begin(), poly.end());
}

Poly indices_to_Poly(const array& indices, const Mon1d& basis);
array Poly_to_indices(const Poly& poly, const Mon1d& basis);

}  // namespace alg

#endif /* DATABASE_H */
