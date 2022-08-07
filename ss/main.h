#ifndef MAIN_H
#define MAIN_H

#include "algebras/benchmark.h"
#include "algebras/dbAdamsSS.h"
#include "algebras/groebner.h"

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

struct NullDiff
{
    AdamsDeg deg;
    unsigned index;
    int first, count;
};

using NullDiff1d = std::vector<NullDiff>;
using NullDiff2d = std::vector<NullDiff1d>;

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

class SS
{
private:
    int t_max_;
    Groebner gb_;
    std::map<AdamsDeg, Mon1d> basis_;
    Staircases1d basis_ss_;
    NullDiff2d nd_;

public:
    SS(Groebner gb, std::map<AdamsDeg, Mon1d> basis, Staircases basis_ss) : gb_(std::move(gb)), basis_(std::move(basis)), basis_ss_({{std::move(basis_ss)}, {}}), nd_({{}, {}})
    {
        t_max_ = basis_.rbegin()->first.t;
    }

public:
    /* Return the newest version of the staircase in history */
    const Staircase& GetRecentStaircase(AdamsDeg deg) const;
    Staircase& GetRecentStaircase(AdamsDeg deg)
    {
        return const_cast<Staircase&>(const_cast<const SS*>(this)->GetRecentStaircase(deg));
    }

    auto& GetBasis() const
    {
        return basis_;
    }

    auto& GetND() const
    {
        return nd_;
    }

    /* Return it is a possible target */
    bool IsPossTgt(AdamsDeg deg, int r) const;

    /* Return it is a possible source */
    bool IsPossSrc(AdamsDeg deg, int r) const;

    /* This is used for plotting Er pages. The actual result might differ by a linear combination. */
    int SS::GetFirstFixedLevelForPlot(AdamsDeg deg) const;

    /* Return the first index (>=`level`) such that all levels above are already fixed */
    size_t GetFirstIndexOfFixedLevels(AdamsDeg deg, int level) const;

    /* Count the number of all possible d_r targets. Return (count, index). */
    std::pair<int, int> CountPossDrTgt(const AdamsDeg& deg_tgt, int r) const;

    /* Count the number of all possible d_r sources. Return (count, index). */
    std::pair<int, int> CountPossDrSrc(const AdamsDeg& deg_src, int r) const;

    /*
     * Count the number of all possible d_r1 targets where r<=r1<=r_max.
     * Return (count, index).
     * count>=2 are all treated as the same.
     */
    std::tuple<AdamsDeg, int, int> CountPossTgt(const AdamsDeg& deg, int r, int r_max) const;

    /*
     * Count the number of all possible d_r1 sources where r1<=r.
     * Return (count, index).
     * count>=2 are all treated as the same.
     */
    std::tuple<AdamsDeg, int, int> CountPossSrc(const AdamsDeg& deg, int level) const;

    /*
     * Return the smallest r1>=r such that d_{r1} has a possible target
     *
     * Range >= t_max is deem unknown and thus possiple.
     * Return -1 if not found
     */
    int NextRTgt(AdamsDeg deg, int r) const;

    /*
     * Return the largest r1<=r such that d_{r1} has a possible source
     *
     * Return -1 if not found
     */
    int NextRSrc(AdamsDeg deg, int r) const;

    int1d GetDiff(AdamsDeg deg_x, int1d x, int r) const;
    bool IsNewDiff(AdamsDeg deg_x, int1d x, int1d dx, int r) const;

    auto& GetChanges() const
    {
        return basis_ss_[1];
    }

private:
    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    void SetDiff(AdamsDeg deg_x, int1d x, int1d dx, int r);

    /* Add an image. dx must be nonempty. */
    void SetImage(AdamsDeg deg_dx, int1d dx, int1d x, int r);

public:
    /* Combine the changes to basis_ss_[index] */
    void ApplyChanges(size_t index);

    /* Combine the last two records in basis_ss */
    void ApplyRecentChanges();

    /* Add a node */
    void AddNode()
    {
        basis_ss_.push_back({});
        nd_.push_back({});
    }

    /* Pop the lastest node */
    void PopNode()
    {
        basis_ss_.pop_back();
        nd_.pop_back();
    }
    /* Apply the change of the staircase to the current history based on the most recent history */
    void UpdateStaircase(AdamsDeg deg, const Staircase& sc_i, size_t i_insert, int1d x, int1d dx, int level, int1d& image, int& level_image);

    /* Cache null diffs to the most recent node. */
    void CacheNullDiffs(int maxPoss);

    /**
     * Add d_r(x)=dx;
     * Add d_s(xy)=d_s(x)y+xd_s(y) for y with level <= kLevelMax - s and s>=r_min.
     * Return the number of changed degrees.
     *
     * dx should not be null.
     */
    int SetDiffLeibniz(AdamsDeg deg_x, int1d x, int1d dx, int r, int r_min);

    /**
     * Check first if it is a new differential before adding it.
     * Do some deductions by degree.
     * 
     * Return the number of changed degrees.
     */
    int SetDiffLeibnizV2(AdamsDeg deg_x, int1d x, int1d dx, int r);

    /* Add d_r(?)=x;
     * Add d_r(?)=xy for d_ry=0 (y in level < kLevelMax - r);
     */
    int SetImageLeibniz(AdamsDeg deg_x, int1d x, int r);

public:
    int DeduceZeroDiffs();
    int DeduceImageJ();
    int DeduceDiffs(int r_max, int maxPoss, int top_depth, int depth, Timer& timer);
};

class SSDB : public myio::DbAdamsSS
{
    using Statement = myio::Statement;

public:
    SSDB() = default;
    explicit SSDB(const std::string& filename) : DbAdamsSS(filename) {}

    void create_basis_ss(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss (id INTEGER PRIMARY KEY, base TEXT, diff TEXT, level SMALLINT, s SMALLINT, t SMALLINT);");
    }

    void create_basis_ss_and_delete(const std::string& table_prefix) const
    {
        create_basis_ss(table_prefix);
        delete_from(table_prefix + "_ss");
    }

    void save_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& basis_ss) const;
    /* load the minimum id in every degree */
    std::map<AdamsDeg, int> load_indices(const std::string& table_prefix) const;
    void update_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& basis_ss) const;
    Staircases load_basis_ss(const std::string& table_prefix) const;

    SS load_ss(const std::string& table_prefix) const;
};

std::ostream& operator<<(std::ostream& sout, const int1d& arr);

int main_basis_prod(int argc, char** argv, int index);
int main_plot(int argc, char** argv, int index);
int main_generate_ss(int argc, char** argv, int index);
int main_add_diff(int argc, char** argv, int index);
int main_try_add_diff(int argc, char** argv, int index);
int main_deduce(int argc, char** argv, int index);

#endif
