#ifndef MAIN_H
#define MAIN_H

#include "algebras/benchmark.h"
#include "algebras/dbAdamsSS.h"
#include "algebras/groebner.h"

inline const char* DB_DEFAULT = "S0_AdamsSS_t255.db";
inline const char* TABLE_DEFAULT = "S0_AdamsE2";
inline const char* VERSION = "Version:\n  2.2 (2022-9-27)";

using namespace alg;

constexpr int kLevelMax = 10000;
constexpr int kLevelMin = 2;
constexpr int kRPC = 200;
constexpr int kLevelPC = kLevelMax - kRPC; /* Level of Permanant cycles */

struct Staircase
{
    int2d basis_ind;
    int2d diffs_ind; /* element {-1} means null */
    int1d levels;
};

using Staircases = std::map<AdamsDeg, Staircase>;
using Staircases1d = std::vector<Staircases>;

/* A custom Exception class */
class SSException : public MyException
{
public:
    SSException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
    SSException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

/* Never to be thrown. A placeholder for some exception to catch */
class NoException : public MyException
{
public:
    NoException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
    NoException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

class Timer : private bench::Timer
{
private:
    double stop_time_;

public:
    Timer(double stop_time) : stop_time_(stop_time)
    {
        SuppressPrint();
    }

    bool timeout() const
    {
        return Elapsed() > stop_time_;
    }
};

struct NullDiff
{
    AdamsDeg deg;
    unsigned index;
    int first, count;
};

using NullDiff1d = std::vector<NullDiff>;
using NullDiff2d = std::vector<NullDiff1d>;

struct SS
{
    int t_max = -1;
    Groebner gb;
    std::map<AdamsDeg, Mon1d> basis;
    Staircases1d basis_ss;
    NullDiff2d nd;
};

struct SSMod
{
    int t_max = -1;
    Poly1d f_top_cell;
    AdamsDeg deg_f_top_cell;
    GroebnerMod gb;
    std::map<AdamsDeg, MMod1d> basis;
    Staircases1d basis_ss;
    NullDiff2d nd;
};

using SSMod1d = std::vector<SSMod>;

class Diagram
{
protected:
    SS ssS0_;
    SSMod1d ssCofs_;

    std::vector<Staircases1d*> all_basis_ss_;
    std::vector<NullDiff2d*> all_nd_;
    int1d all_t_max_;

public:
    Diagram(const std::vector<std::string>& dbnames);

public:
    /* Return the newest version of the staircase in history */
    const Staircase& GetRecentStaircase(const Staircases1d& basis_ss, AdamsDeg deg) const;
    Staircase& GetRecentStaircase(Staircases1d& basis_ss, AdamsDeg deg)
    {
        return const_cast<Staircase&>(const_cast<const Diagram*>(this)->GetRecentStaircase(basis_ss, deg));
    }

    auto& GetS0() const
    {
        return ssS0_;
    }

    auto& GetCofs() const
    {
        return ssCofs_;
    }

    auto& GetAllBasisSs() const
    {
        return all_basis_ss_;
    }

    auto& GetAllNd() const
    {
        return all_nd_;
    }

    auto& GetAllTMax() const
    {
        return all_t_max_;
    }

    /* Return it is a possible target */
    bool IsPossTgt(const Staircases1d& basis_ss, AdamsDeg deg, int r) const;

    /* Return it is a possible source */
    bool IsPossSrc(const Staircases1d& basis_ss, int t_max, AdamsDeg deg, int r) const;

    /* This is used for plotting Er pages. The actual result might differ by a linear combination.
     * Return a level such that all levels above will not decrease further.
     */
    int Diagram::GetFirstFixedLevelForPlot(const Staircases1d& basis_ss, AdamsDeg deg) const;

    /* Return the first index (>=`level`) such that all levels above are already fixed */
    size_t GetFirstIndexOfFixedLevels(const Staircases1d& basis_ss, AdamsDeg deg, int level) const;

    /* Count the number of all possible d_r targets. Return (count, index). */
    std::pair<int, int> CountPossDrTgt(const Staircases1d& basis_ss, int t_max, const AdamsDeg& deg_tgt, int r) const;

    /* Count the number of all possible d_r sources. Return (count, index). */
    std::pair<int, int> CountPossDrSrc(const Staircases1d& basis_ss, const AdamsDeg& deg_src, int r) const;

    /*
     * Count the number of all possible d_r1 targets where r<=r1<=r_max.
     * Return (count, index).
     * count>=2 are all treated as the same.
     */
    std::tuple<AdamsDeg, int, int> CountPossTgt(const Staircases1d& basis_ss, int t_max, const AdamsDeg& deg, int r, int r_max) const;

    /*
     * Count the number of all possible d_r1 sources where r1<=r.
     * Return (count, index).
     * count>=2 are all treated as the same.
     */
    std::tuple<AdamsDeg, int, int> CountPossSrc(const Staircases1d& basis_ss, AdamsDeg deg, int level) const;

    /*
     * Return the smallest r1>=r such that d_{r1} has a possible target
     *
     * Range >= t_max is deem unknown and thus possiple.
     * Return -1 if not found
     */
    int NextRTgt(const Staircases1d& basis_ss, int t_max, AdamsDeg deg, int r) const;

    /*
     * Return the largest r1<=r such that d_{r1} has a possible source
     *
     * Return -1 if not found
     */
    int NextRSrc(const Staircases1d& basis_ss, AdamsDeg deg, int r) const;

    int1d GetDiff(const Staircases1d& basis_ss, AdamsDeg deg_x, int1d x, int r) const;
    bool IsNewDiff(const Staircases1d& basis_ss, AdamsDeg deg_x, int1d x, int1d dx, int r) const;

    auto& GetChanges(size_t iBasisSS) const
    {
        return (*all_basis_ss_[iBasisSS])[1];
    }

protected:
    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    void SetDiff(Staircases1d& basis_ss, AdamsDeg deg_x, int1d x, int1d dx, int r);

    /* Add an image. dx must be nonempty. */
    void SetImage(Staircases1d& basis_ss, AdamsDeg deg_dx, int1d dx, int1d x, int r);

public:
    /* Combine the changes to basis_ss_[index] */
    void ApplyChanges(size_t index);

    /* Combine the last two records in basis_ss */
    void ApplyRecentChanges();

    /* Add a node */
    void AddNode()
    {
        for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
            all_basis_ss_[k]->push_back({});
            all_nd_[k]->push_back({});
        }
    }

    /* Pop the lastest node */
    void PopNode()
    {
        for (size_t k = 0; k < all_basis_ss_.size(); ++k) {
            all_basis_ss_[k]->pop_back();
            all_nd_[k]->pop_back();
        }
    }
    /* Apply the change of the staircase to the current history based on the most recent history */
    void UpdateStaircase(Staircases1d& basis_ss, AdamsDeg deg, const Staircase& sc_i, size_t i_insert, int1d x, int1d dx, int level, int1d& image, int& level_image);

    /* Cache null diffs to the most recent node. */
    void CacheNullDiffs(int maxPoss);

    /**
     * Add d_r(x)=dx;
     * Add d_s(xy)=d_s(x)y+xd_s(y) for y with level=kLevelMax-s and s>=r_min.
     * Return the number of changed degrees.
     *
     * dx should not be null.
     */
    int SetS0DiffLeibniz(AdamsDeg deg_x, int1d x, int1d dx, int r, int r_min);
    int SetCofDiffLeibniz(size_t iCof, AdamsDeg deg_x, int1d x, int1d dx, int r, int r_min);

    /**
     * Check first if it is a new differential before adding it.
     * Do some deductions by degree.
     *
     * Return the number of changed degrees.
     */
    int SetS0DiffLeibnizV2(AdamsDeg deg_x, int1d x, int1d dx, int r);
    int SetCofDiffLeibnizV2(size_t iCof, AdamsDeg deg_x, int1d x, int1d dx, int r);
    int SetDiffLeibnizV2(size_t index, AdamsDeg deg_x, int1d x, int1d dx, int r);

    /* Add d_r(?)=x;
     * Add d_r(?)=xy for d_ry=0 (y in level < kLevelMax - r);
     */
    int SetS0ImageLeibniz(AdamsDeg deg_x, int1d x, int r);
    int SetCofImageLeibniz(size_t iCof, AdamsDeg deg_x, int1d x, int r);

public:
    int DeduceZeroDiffs();
    int DeduceDiffs(int r_max, int maxPoss, int top_depth, int depth, Timer& timer);
    int DeduceImageJ();
};

class DBSS : public myio::DbAdamsSS
{
    using Statement = myio::Statement;

public:
    DBSS() = default;
    explicit DBSS(const std::string& filename) : DbAdamsSS(filename) {}

    void create_basis_ss(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss (id INTEGER PRIMARY KEY, base TEXT, diff TEXT, level SMALLINT, s SMALLINT, t SMALLINT);");
    }

    void drop_and_create_basis_ss(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss");
        create_basis_ss(table_prefix);
    }

    void save_basis_ss(const std::string& table_prefix, const Staircases& basis_ss) const;
    /* load the minimum id in every degree */
    std::map<AdamsDeg, int> load_indices(const std::string& table_prefix) const;
    void update_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& basis_ss) const;
    Staircases load_basis_ss(const std::string& table_prefix) const;
};

std::ostream& operator<<(std::ostream& sout, const int1d& arr);
std::string GetTablePrefix(const std::string& db);
int GetTopCellT(const std::string& db);

template <typename T>
AdamsDeg1d OrderDegsV2(const T& cont)
{
    AdamsDeg1d result;
    for (auto& [d, _] : cont) {
        result.push_back(d);
    }
    std::sort(result.begin(), result.end(), [](const AdamsDeg& d1, const AdamsDeg& d2) { return d1.t < d2.t || (d1.t == d2.t && d1.s > d2.s); });
    return result;
}

int main_basis_prod(int argc, char** argv, int index);
int main_plot(int argc, char** argv, int index);
int main_generate_ss(int argc, char** argv, int index);
int main_add_diff(int argc, char** argv, int index);
int main_try_add_diff(int argc, char** argv, int index);

int main_deduce(int argc, char** argv, int index);
int main_deduce_migrate(int argc, char** argv, int index);

int main_mod(int argc, char** argv, int index);

#endif
