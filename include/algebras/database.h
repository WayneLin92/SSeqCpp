/** \file database.h
 * A component for interaction with a database.
 * TODO: Template member functions.
 * TODO: Change load_db/save_db.
 */
#ifndef DATABASE_H
#define DATABASE_H

#define DATABASE_SAVE_LOGGING /* This is a switch to print what are saved to the database */

#include "algebras.h"
#include <map>
#include <string>
#include <vector>

/********** STRUCTS AND CLASSES **********/
/* Declarations */
constexpr auto kLevelMax = 10000;
struct sqlite3;
struct sqlite3_stmt;

struct BasisComplex
{
    alg::array2d boundaries;
    alg::array2d cycles;
};

struct Staircase
{
    alg::array2d basis_ind;
    alg::array2d diffs_ind; /* element {-1} means null */
    alg::array levels;
};

/* The wrapper for sqlite* */
class Database
{
public:
    Database() = default;
    explicit Database(const char* filename);
    ~Database();

public:
    void execute_cmd(const std::string& sql) const;
    int get_int(const std::string& sql) const;
    alg::array get_ints(const std::string& table_name, const std::string& column_name, const std::string& conditions = "") const;
    void sqlite3_prepare_v100(const char* zSql, sqlite3_stmt** ppStmt, bool bPrintError = false) const;
    void sqlite3_prepare_v100(const std::string& sql, sqlite3_stmt** ppStmt, bool bPrintError = false) const;

public:
    void begin_transaction() const
    {
        execute_cmd("BEGIN TRANSACTION");
    }
    void end_transaction() const
    {
        execute_cmd("END TRANSACTION");
    }

public:
    std::vector<alg::Deg> load_gen_degs(const std::string& table_name) const;
    std::vector<std::string> load_gen_names(const std::string& table_name) const;
    alg::Poly1d load_gen_diffs(const std::string& table_name) const;
    alg::Poly1d load_gen_reprs(const std::string& table_name) const;
    alg::Poly1d load_gen_images(const std::string& table_name, const std::string& column_name, int t_max) const;
    alg::Mon2d load_leading_terms(const std::string& table_name, int t_max) const;
    alg::Poly1d load_gb(const std::string& table_name, int t_max) const;
    std::map<alg::Deg, int> load_indices(const std::string& table_name, int t_max) const;
    std::map<alg::Deg, alg::Mon1d> load_basis(const std::string& table_name, int t_max) const;
    std::map<alg::Deg, alg::array2d> load_mon_diffs_ind(const std::string& table_name, int t_max, int withnull = false) const;
    std::map<alg::Deg, alg::Poly1d> load_mon_diffs(const std::string& table_name, const std::map<alg::Deg, alg::Mon1d>& basis, int r, int t_max) const;
    std::map<alg::Deg, BasisComplex> load_basis_ss(const std::string& table_name_ss, int r, int t_max) const; /* load E_r-cycles and E_r-boundaries */
    std::map<alg::Deg, Staircase> load_basis_ss(const std::string& table_name_ss, int t_max) const;

public:
    void save_generators(const std::string& table_name, const std::vector<std::string>& gen_names, const std::vector<alg::Deg>& gen_degs) const;
    void save_generators(const std::string& table_name, const std::vector<alg::Deg>& gen_degs, const alg::Poly1d& gen_reprs, size_t i_start) const;
    void save_generators(const std::string& table_name, const std::vector<alg::Deg>& gen_degs, const alg::Poly1d& gen_reprs) const
    {
        save_generators(table_name, gen_degs, gen_reprs, 0);
    }
    void save_generators(const std::string& table_name, const std::vector<std::string>& gen_names, const std::vector<alg::Deg>& gen_degs, const alg::Poly1d& gen_reprs) const;
    void save_gen_images(const std::string& table_name, const std::string& column_name, const alg::Poly1d& gen_images) const;
    void save_gb(const std::string& table_name, const alg::Poly1d& gb, const std::vector<alg::Deg>& gen_degs, size_t i_start) const;
    void save_gb(const std::string& table_name, const alg::Poly1d& gb, const std::vector<alg::Deg>& gen_degs) const
    {
        save_gb(table_name, gb, gen_degs, 0);
    }
    void save_basis(const std::string& table_name, const std::map<alg::Deg, alg::Mon1d>& basis) const;
    void save_basis(const std::string& table_name, const std::map<alg::Deg, alg::Mon1d>& basis, const std::map<alg::Deg, alg::array2d>& mon_reprs) const;
    void save_basis(const std::string& table_name, const std::map<alg::Deg, alg::Mon1d>& basis, const std::map<alg::Deg, alg::Poly1d>& mon_reprs) const;
    void save_basis_ss(const std::string& table_name, const std::map<alg::Deg, Staircase>& basis_ss) const;
    void update_basis_ss(const std::string& table_name, const std::map<alg::Deg, Staircase>& basis_ss) const;

private:
    sqlite3* conn_ = nullptr;
};

/* The wrapper for sqlite3_stmt* */
class Statement
{
public:
    Statement() = default;
    explicit Statement(const Database& db, const std::string& sql);
    ~Statement();

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

private:
    sqlite3_stmt* stmt_ = nullptr;
};

alg::array str_to_array(const char* str_mon);
alg::Mon str_to_Mon(const char* str_mon);
alg::Poly str_to_Poly(const char* str_poly);
std::string array_to_str(alg::array::const_iterator pbegin, alg::array::const_iterator pend);
inline std::string array_to_str(const alg::array& a)
{
    return array_to_str(a.begin(), a.end());
}
std::string Mon_to_str(alg::MonInd pbegin, alg::MonInd pend);
inline std::string Mon_to_str(const alg::Mon& mon)
{
    return Mon_to_str(mon.begin(), mon.end());
}
std::string Poly_to_str(alg::Poly::const_iterator pbegin, alg::Poly::const_iterator pend);
inline std::string Poly_to_str(const alg::Poly& poly)
{
    return Poly_to_str(poly.begin(), poly.end());
}

alg::Poly indices_to_Poly(const alg::array& indices, const alg::Mon1d& basis);
alg::array Poly_to_indices(const alg::Poly& poly, const alg::Mon1d& basis);

#endif /* DATABASE_H */
