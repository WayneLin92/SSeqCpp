#ifndef MAIN_H
#define MAIN_H

#include "algebras/benchmark.h"
#include "algebras/dbAdamsSS.h"
#include "algebras/groebnerZ.h"
#include <set>

inline const char* VERSION = "Version:\n  3.0 (2022-11-14)";
using namespace alg2;

constexpr int kLevelMax = 10000;
constexpr int kLevelMin = 2;
constexpr int kRPC = 200;
constexpr int kLevelPC = kLevelMax - kRPC; /* Level of Permanant cycles */

enum class DeduceFlag : uint32_t
{
    no_op = 0,
    set_diff = 1,
    homotopy = 2,
    fast_try_diff = 4,   /* SetDiffLeibniz will update only partially in the try node */
    check_exactness = 8, /* Check exactness of htpy in the try node */
    all_x = 16,          /* Deduce dx for all x including linear combinations */
};

inline DeduceFlag operator|(DeduceFlag lhs, DeduceFlag rhs)
{
    return DeduceFlag(uint32_t(lhs) | uint32_t(rhs));
}

inline bool operator&(DeduceFlag lhs, DeduceFlag rhs)
{
    return uint32_t(lhs) & uint32_t(rhs);
}

struct Staircase
{
    int2d basis_ind;
    int2d diffs_ind; /* element {-1} means null */
    int1d levels;
};

using Staircases = std::map<AdamsDeg, Staircase>;
using Staircases1d = std::vector<Staircases>;

struct PiBase
{
    algZ::Mon1d pi_basis;
    int2d Einf;
};

using PiBasis = std::map<AdamsDeg, PiBase>;
using PiBasis1d = std::vector<PiBasis>;

struct PiBaseMod
{
    algZ::MMod1d pi_basis;
    int2d Einf;
};

using PiBasisMod = std::map<AdamsDeg, PiBaseMod>;
using PiBasisMod1d = std::vector<PiBasisMod>;

/* A custom Exception class */
class SSException : public MyException
{
public:
    SSException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
    SSException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

/* A custom Exception class */
class SSPiException : public SSException
{
public:
    SSPiException(unsigned int error_id, const char* message) : SSException(error_id, message) {}
    SSPiException(unsigned int error_id, const std::string& message) : SSException(error_id, message) {}
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
    int1d x;
    int r;
    int first, count;
};

using NullDiff1d = std::vector<NullDiff>;

struct SS
{
    std::string name; /* Constant after initialization */
    int t_max = -1;
    Groebner gb;
    std::map<AdamsDeg, Mon1d> basis;
    AdamsDeg1d degs_basis_order_by_stem; /* Constant after initialization */
    Staircases1d basis_ss;

    algZ::Groebner pi_gb;
    Poly1d pi_gen_Einf = {Poly::Gen(0)}; /* Projection onto the E_infty page */
    PiBasis1d pi_basis = {{{AdamsDeg(0, 0), {{algZ::Mon()}, {{0}}}}}};
    std::vector<size_t> pi_nodes_gen;
    std::vector<size_t> pi_nodes_rel;
    int2d pi_nodes_gen_2tor_degs;
};

struct SSMod
{
    std::string name; /* Constant after initialization */
    int t_max = -1;
    AdamsDeg deg_qt;
    GroebnerMod gb;
    std::map<AdamsDeg, MMod1d> basis;
    AdamsDeg1d degs_basis_order_by_stem; /* Constant after initialization */
    Poly1d qt;
    Staircases1d basis_ss;

    algZ::GroebnerMod pi_gb;
    Mod1d pi_gen_Einf; /* Projection onto the E_infty page */
    PiBasisMod1d pi_basis = {{}};
    algZ::Poly2d pi_qt;
    std::vector<size_t> pi_nodes_gen;
    std::vector<size_t> pi_nodes_rel;
};

using SSMod1d = std::vector<SSMod>;

constexpr size_t MAX_NUM_NODES = 5;

class Diagram
{

public: /* Settings */
    int stem_max_exactness_ = 100;

protected:
    SS ssS0_;
    SSMod1d ssCofs_;

    std::vector<std::string> all_names_;
    std::vector<Staircases1d*> all_basis_ss_;
    int1d all_t_max_;

public:
    Diagram(const std::vector<std::string>& dbnames);

public:
    /* Return the newest version of the staircase in history */
    static auto GetRecentStaircase(const Staircases1d& basis_ss, AdamsDeg deg) -> const Staircase&;
    static auto GetRecentStaircase(Staircases1d& basis_ss, AdamsDeg deg) -> Staircase&
    {
        return const_cast<Staircase&>(GetRecentStaircase((const Staircases1d&)(basis_ss), deg));
    }
    static auto GetRecentPiBasis(const PiBasis1d& pi_basis, AdamsDeg deg) -> const PiBase*;
    static auto GetRecentPiBasis(const PiBasisMod1d& pi_basis, AdamsDeg deg) -> const PiBaseMod*;

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

    auto& GetAllTMax() const
    {
        return all_t_max_;
    }

    /* Return if it is possibly a new dr target for r<=r_max */
    static bool IsPossTgt(const Staircases1d& basis_ss, AdamsDeg deg, int r_max);

    /* Return if it is possibly a new dr source for r>=r_min */
    static bool IsPossSrc(const Staircases1d& basis_ss, int t_max, AdamsDeg deg, int r_min);

    /* This is used for plotting Er pages. The actual result might differ by a linear combination.
     * Return a level such that all levels above will not decrease further.
     */
    static int GetFirstFixedLevelForPlot(const Staircases1d& basis_ss, AdamsDeg deg);

    /* Return the first index (with level >=`level_min`) such that all levels above are already fixed */
    static size_t GetFirstIndexOfFixedLevels(const Staircases1d& basis_ss, AdamsDeg deg, int level_min);

    /* Count the number of all possible d_r targets. Return (count, index). */
    auto CountPossDrTgt(const Staircases1d& basis_ss, int t_max, const AdamsDeg& deg_tgt, int r) const -> std::pair<int, int>;

    /* Count the number of all possible d_r sources. Return (count, index). */
    std::pair<int, int> CountPossDrSrc(const Staircases1d& basis_ss, const AdamsDeg& deg_src, int r) const;

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

    int1d GetDiff(const Staircases1d& basis_ss, AdamsDeg deg_x, const int1d& x, int r) const;
    bool IsNewDiff(const Staircases1d& basis_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const;

    auto& GetBasisSSChanges(size_t iBasisSS) const
    {
        if ((*all_basis_ss_[iBasisSS]).size() != 2)
            throw MyException(0xc2fa755cU, "Not on the change node");
        return (*all_basis_ss_[iBasisSS])[1];
    }

protected:
    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    void SetDiff(Staircases1d& basis_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r);

    /* Add an image. dx must be nonempty. */
    void SetImage(Staircases1d& basis_ss, AdamsDeg deg_dx, const int1d& dx, const int1d& x, int r);

public:
    /* Add a node */
    void AddNode();

    /* Pop the lastest node */
    void PopNode();

    /* Apply the change of the staircase to the current history */
    void UpdateStaircase(Staircases1d& basis_ss, AdamsDeg deg, const Staircase& sc_i, size_t i_insert, const int1d& x, const int1d& dx, int level, int1d& image, int& level_image);

    /* Cache null diffs to the most recent node. */
    void CacheNullDiffs(size_t iSS, AdamsDeg deg, DeduceFlag flag, NullDiff1d& nds);

    /**
     * Add d_r(x)=dx;
     * Add d_s(xy)=d_s(x)y+xd_s(y) for y with level=kLevelMax-s and s>=r_min.
     * Return the number of changed degrees.
     *
     * dx should not be null.
     */
    int SetS0DiffLeibniz(AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, bool bFastTry = false);
    int SetCofDiffLeibniz(size_t iCof, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, bool bFastTry = false);

    /**
     * Check first if it is a new differential before adding it.
     * Do some deductions by degree.
     *
     * Return the number of changed degrees.
     */
    int SetS0DiffLeibnizV2(AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry = false);
    int SetCofDiffLeibnizV2(size_t iCof, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry = false);
    int SetDiffLeibnizV2(size_t iSS, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry = false);

    /* Add d_r(?)=x;
     * Add d_r(?)=xy for d_ry=0 (y in level < kLevelMax - r);
     */
    int SetS0ImageLeibniz(AdamsDeg deg_x, const int1d& x, int r);
    int SetCofImageLeibniz(size_t iCof, AdamsDeg deg_x, const int1d& x, int r);

public:
    int DeduceTrivialDiffs();
    /* Return 0 if there is no exception */
    int TryDiff(size_t iSS, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int depth, DeduceFlag flag);
    int DeduceDiffs(size_t iSS, AdamsDeg deg, int depth, DeduceFlag flag);
    int DeduceDiffs(int stem_min, int stem_max, int depth, DeduceFlag flag);

public:
    /* Return if Einf at deg is possibly nontrivial */
    static int PossEinf(const Staircases1d& basis_ss, AdamsDeg deg);
    /* Return if Einf at deg could have more nontrivial elements */
    static int PossMoreEinf(const Staircases1d& basis_ss, AdamsDeg deg);
    int1d PossMoreEinfFirstS_S0() const;
    int1d PossMoreEinfFirstS_Cof(size_t iCof) const;

    /*
     * Return the smallest s1>=s such that extension has a possible target
     *
     * Range >= t_max is deem unknown and thus possiple.
     * Return FIL_MAX if not found
     */
    int ExtendRelS0(int stem, const algZ::Poly& rel, algZ::Poly& rel_extended) const;
    int ExtendRelCof(size_t iCof, int stem, const algZ::Mod& rel, algZ::Mod& rel_extended) const;
    void ExtendRelS0(int stem, algZ::Poly& rel) const;
    void ExtendRelCof(size_t iCof, int stem, algZ::Mod& rel) const;
    void ExtendRelS0V2(int stem, algZ::Poly& rel, std::unordered_map<int, int>& num_leads) const;
    void ExtendRelCofV2(size_t iCof, int stem, algZ::Mod& rel, std::unordered_map<int, int>& num_leads) const;
    int2d GetS0GbEinf(AdamsDeg deg) const;
    std::map<AdamsDeg, int2d> GetS0GbEinf() const;
    int2d GetCofGbEinf(int iCof, AdamsDeg deg) const;
    std::map<AdamsDeg, int2d> GetCofGbEinf(int iCof) const;

public: /* homotopy groups */
    void AddPiRelsS0(algZ::Poly1d rels);
    void AddPiRelsCof(size_t iCof, algZ::Mod1d rels);
    void AddPiRelsCof2S0(size_t iCof);
    void SimplifyPiRels()
    {
        ssS0_.pi_gb.SimplifyRels();
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
            ssCofs_[iCof].pi_gb.SimplifyRels();
    }
    static algZ::Mon1d GenBasis(const algZ::Groebner& gb, AdamsDeg deg, const PiBasis1d& pi_basis);
    static algZ::MMod1d GenBasis(const algZ::GroebnerMod& gb, AdamsDeg deg, const PiBasis1d& basis);
    void SyncS0Homotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopyy, int depth);
    void SyncCofHomotopy(int iCof, AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth);
    void SyncHomotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth)
    {
        SyncS0Homotopy(deg_min, count_ss, count_homotopy, depth);
        for (size_t iCof = 0; iCof < ssCofs_.size(); ++iCof)
            SyncCofHomotopy((int)iCof, deg_min, count_ss, count_homotopy, depth);
    }
    int DeduceTrivialExtensions(int depth);
    int DeduceExtensionsByExactness(int stem_min, int stem_max, int depth);

    unsigned TryExtS0(algZ::Poly rel, AdamsDeg deg_change, int depth, DeduceFlag flag);
    unsigned TryExtCof(size_t iCof, algZ::Mod rel, AdamsDeg deg_change, int depth, DeduceFlag flag);
    unsigned TryExtQ(size_t iCof, size_t gen_id, algZ::Poly q, AdamsDeg deg_change, int depth, DeduceFlag flag);
    void DeduceExtensions(int stem_min, int stem_max, int& count_ss, int& count_homotopy, int depth, DeduceFlag flag);
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

    void create_pi_generators_mod(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, Einf TEXT, to_S0 TEXT, s SMALLINT, t SMALLINT);");
    }

    void drop_and_create_basis_ss(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss");
        create_basis_ss(table_prefix);
    }

    void drop_and_create_pi_generators_mod(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_pi_generators");
        create_pi_generators_mod(table_prefix);
    }

    void save_pi_generators_mod(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Mod1d& gen_Einf, const algZ::Poly1d& to_S0) const;
    void save_basis_ss(const std::string& table_prefix, const Staircases& basis_ss) const;
    void save_pi_basis(const std::string& table_prefix, const PiBasis& basis) const;
    void save_pi_basis_mod(const std::string& table_prefix, const PiBasisMod& basis) const;
    /* load the minimum id in every degree */
    std::map<AdamsDeg, int> load_basis_indices(const std::string& table_prefix) const;
    void update_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& basis_ss) const;
    Staircases load_basis_ss(const std::string& table_prefix) const;
};

std::ostream& operator<<(std::ostream& sout, const int1d& arr);
std::string GetComplexName(const std::string& db);
std::string GetE2TablePrefix(const std::string& db);
int GetTopCellT(const std::string& db);

/* Order by (t, -s) */
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

/* Order by (t, -s) */
template <typename T>
AdamsDeg1d OrderDegsByStem(const T& cont)
{
    AdamsDeg1d result;
    for (auto& [d, _] : cont) {
        result.push_back(d);
    }
    std::sort(result.begin(), result.end(), [](const AdamsDeg& d1, const AdamsDeg& d2) { return d1.stem() < d2.stem() || (d1.stem() == d2.stem() && d1.s < d2.s); });
    return result;
}

/* If n = 2^k1 + ... + 2^kn,
 * return the array k1, ..., kn. */
inline int1d two_expansion(unsigned n)
{
    int1d result;
    int k = 0;
    while (n > 0) {
        if (n & 1)
            result.push_back(k);
        n >>= 1;
        ++k;
    }
    return result;
}

inline bool BelowS0VanishingLine(AdamsDeg deg)
{
    return 3 * deg.s <= deg.t + 3;
}

size_t GetFirstIndexOnLevel(const Staircase& sc, int level);

std::vector<std::string> GetDbNames(const std::string& selector);

int main_basis_prod(int argc, char** argv, int index);
int main_plot(int argc, char** argv, int index);
int main_plotpi(int argc, char** argv, int index);
int main_reset(int argc, char** argv, int index);
int main_add_diff(int argc, char** argv, int index);
int main_add_ext(int argc, char** argv, int index);

int main_deduce(int argc, char** argv, int index);
int main_deduce_ext(int argc, char** argv, int index);
int main_deduce_migrate(int argc, char** argv, int index);

int main_mod(int argc, char** argv, int index);

int main_resetpi(int argc, char** argv, int index);
int main_truncate(int argc, char** argv, int index);

#endif
