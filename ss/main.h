#ifndef MAIN_H
#define MAIN_H

#include "algebras/benchmark.h"
#include "algebras/dbAdamsSS.h"
#include "pigroebner.h"
#include <set>
#include <variant>

inline const char* VERSION = "Version:\n  3.1 (2023-01-11)";
using namespace alg2;

constexpr int LEVEL_MAX = 10000;
constexpr int LEVEL_MIN = 2;
constexpr int R_PERM = 200;
constexpr int LEVEL_PERM = LEVEL_MAX - R_PERM; /* Level of Permanant cycles */

constexpr size_t MAX_DEPTH = 3; /* Maximum deduction depth */

inline const auto NULL_DIFF = int1d{-1};
inline const algZ::Mod MOD_V0 = algZ::MMod(algZ::Mon(), 0, 0);

enum class DeduceFlag : uint32_t
{
    no_op = 0,
    set_diff = 1,
    fast_try_diff = 2, /* SetDiffLeibniz will update only partially in the try node */
    all_x = 4,         /* Deduce dx for all x including linear combinations */
    homotopy = 8,
    homotopy_exact = 16, /* Check exactness of htpy in the try node */
    homotopy_def = 32,   /* Check exactness of htpy in the try node */
};

inline DeduceFlag operator|(DeduceFlag lhs, DeduceFlag rhs)
{
    return DeduceFlag(uint32_t(lhs) | uint32_t(rhs));
}

inline bool operator&(DeduceFlag lhs, DeduceFlag rhs)
{
    return uint32_t(lhs) & uint32_t(rhs);
}

enum class EnumDef : int
{
    no_def = 0,
    dec = 1,         /* Decomposable */
    constraints = 2, /* The indeterminancies have some constraints */
};

struct Staircase
{
    int2d basis;
    int2d diffs; /* element {-1} means null */
    int1d levels;
};

using Staircases = std::map<AdamsDeg, Staircase>;
using Staircases1d = std::vector<Staircases>;

struct PiBase
{
    algZ::Mon1d nodes_pi_basis;
    int2d Einf;
};

using PiBasis = std::map<AdamsDeg, PiBase>;
using PiBasis1d = std::vector<PiBasis>;

struct PiBaseMod
{
    algZ::MMod1d nodes_pi_basis;
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

/* Never to be thrown. A dummy Exception to be filled in the catch statement */
class NoException : public MyException
{
public:
    NoException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
    NoException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

class TerminationException : public MyException
{
public:
    TerminationException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
    TerminationException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

struct NullDiff
{
    int1d x;
    int r;
    int first, count;
};

using NullDiff1d = std::vector<NullDiff>;

/* This constrains the choice of a generator g by requiring that Unknownfil(g*m)=O. */
struct GenConstraint
{
    int map_index;
    algZ::Mon m;
    int O;
};

struct RingSp
{
    /* #metadata */
    std::string name;
    int t_max = -1;  ////TODO: remove
    std::vector<size_t> ind_mods, ind_maps;

    /* #ss */
    Groebner gb;
    std::map<AdamsDeg, Mon1d> basis;
    AdamsDeg1d degs_basis_order_by_stem;      /* Constant after initialization */
    Staircases1d nodes_ss;                    /* size = depth + 2 */
    ut::map_seq2d<int, 0> basis_ss_possEinf;  ////TODO: change to int2d

    /* #pi */
    algZ::Groebner pi_gb;
    Poly1d pi_gen_Einf = {Poly::Gen(0)};
    PiBasis1d nodes_pi_basis = {{{AdamsDeg(0, 0), {{algZ::Mon()}, {{0}}}}}}; /* size = depth + 1 */

    std::vector<EnumDef> pi_gen_defs;
    std::vector<std::vector<GenConstraint>> pi_gen_def_mons;
};
using RingSp1d = std::vector<RingSp>;

struct ModSp
{
    /* #metadata */
    std::string name;
    int t_max = -1;  ////TODO: remove
    AdamsDeg deg_qt;
    size_t iRing;
    std::vector<size_t> ind_maps;

    /* #ss */
    GroebnerMod gb;
    std::map<AdamsDeg, MMod1d> basis;
    AdamsDeg1d degs_basis_order_by_stem;      /* Constant after initialization */
    Staircases1d nodes_ss;                    /* size = depth + 2 */
    ut::map_seq2d<int, 0> basis_ss_possEinf;  ////TODO: change to int2d

    /* #pi */
    algZ::GroebnerMod pi_gb;
    Mod1d pi_gen_Einf;
    PiBasisMod1d nodes_pi_basis = {{}}; /* size = depth + 1 */

    std::vector<EnumDef> pi_gen_defs;
    std::vector<std::vector<GenConstraint>> pi_gen_def_mons;
};
using ModSp1d = std::vector<ModSp>;

struct MapRing2Ring  ////TODO: Add pi_images
{
    size_t from, to;
    std::vector<Poly> images;
    int1d map(const int1d& x, AdamsDeg deg_x, const RingSp1d& rings);
};

struct MapMod2Mod
{
    size_t from, to;
    int sus;
    std::vector<Mod> images;
    int1d map(const int1d& x, AdamsDeg deg_x, const ModSp1d& mods);
};

struct MapMod2Ring
{
    size_t from, to;
    int sus;
    std::vector<Poly> images;
    int1d map(const int1d& x, AdamsDeg deg_x, const ModSp1d& mods, const RingSp1d& rings);
};

struct MapMulRing
{
    size_t index;
    int sus;
    Poly factor;
};

struct MapRing2Mod
{
    size_t from, to;
    int sus;
    Mod factor;
};

struct Map
{
    std::string name;
    std::variant<MapRing2Ring, MapMod2Mod, MapMod2Ring, MapMulRing, MapRing2Mod> map;
};

using Map1d = std::vector<Map>;

class Diagram
{

public: /* Settings */
protected:
    RingSp1d rings_;
    ModSp1d modules_;
    Map1d maps_;

public:
    Diagram(std::string diagram_name, DeduceFlag flag, bool log = true);
    void VersionConvertReorderRels()
    {
        for (auto& ring : rings_)
            ring.pi_gb.SimplifyRelsReorder(ring.basis_ss_possEinf);
        for (auto& mod : modules_)
            mod.pi_gb.SimplifyRelsReorder(mod.basis_ss_possEinf);
    }
    void save(std::string diagram_name, DeduceFlag flag);

public: /* Getters */
    /* Return the newest version of the staircase in history */
    static auto GetRecentSc(const Staircases1d& nodes_ss, AdamsDeg deg) -> const Staircase&;
    static auto GetRecentSc(Staircases1d& nodes_ss, AdamsDeg deg) -> Staircase&
    {
        return const_cast<Staircase&>(GetRecentSc((const Staircases1d&)(nodes_ss), deg));
    }
    static auto GetRecentPiBasis(const PiBasis1d& nodes_pi_basis, AdamsDeg deg) -> const PiBase*;
    static auto GetRecentPiBasis(const PiBasisMod1d& nodes_pi_basis, AdamsDeg deg) -> const PiBaseMod*;

    auto& GetRings() const
    {
        return rings_;
    }

    auto& GetRings()
    {
        return rings_;
    }

    auto& GetModules() const
    {
        return modules_;
    }

    auto& GetModules()
    {
        return modules_;
    }

    int GetRingIndexByName(const std::string& name) const
    {
        for (size_t i = 0; i < rings_.size(); ++i)
            if (rings_[i].name == name)
                return (int)i;
        return -1;
    }

    int GetModuleIndexByName(const std::string& name) const
    {
        for (size_t i = 0; i < modules_.size(); ++i)
            if (modules_[i].name == name)
                return (int)i;
        return -1;
    }

    auto& GetRingByName(const std::string& name) const
    {
        for (auto& ring : rings_)
            if (ring.name == name)
                return ring;
        throw MyException(0xdb4dc253, "Incorrect name");
    }

    auto& GetModuleByName(const std::string& name) const
    {
        for (auto& mod : modules_)
            if (mod.name == name)
                return mod;
        throw MyException(0x8d0358d4, "Incorrect name");
    }

    /* This is used for plotting Er pages. The actual result might differ by a linear combination.
     * Return a level such that all levels above will not decrease further.
     */
    static int GetFirstFixedLevelForPlot(const Staircases1d& nodes_ss, AdamsDeg deg);

private: /* Staircase */
    static size_t GetFirstIndexOfNullOnLevel(const Staircase& sc, int level);
    static int GetMaxLevelWithNull(const Staircase& sc);
    static bool IsZeroOnLevel(const Staircase& sc, const int1d& x, int level);

private: /* ss */
    /* Return if it is possibly a new dr target for r<=r_max
     * Warning: This does check if ss[deg] is trivial
     */
    static bool IsPossTgt(const Staircases1d& nodes_ss, AdamsDeg deg, int r_max);

    static bool IsPossTgt(const Staircases1d& nodes_ss, AdamsDeg deg)
    {
        return IsPossTgt(nodes_ss, deg, R_PERM);
    }

    /* Return if it is possibly a new dr source for r>=r_min */
    static bool IsPossSrc(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, int r_min);

    /* Return the first index (with level >=`level_min`) such that all levels above are already fixed */
    static size_t GetFirstIndexOfFixedLevels(const Staircases1d& nodes_ss, AdamsDeg deg, int level_min);

    /* Count the number of all possible d_r targets. Return (count, index). */
    auto CountPossDrTgt(const Staircases1d& nodes_ss, int t_max, const AdamsDeg& deg_tgt, int r) const -> std::pair<int, int>;

    /* Count the number of all possible d_r sources. Return (count, index). */
    std::pair<int, int> CountPossDrSrc(const Staircases1d& nodes_ss, const AdamsDeg& deg_src, int r) const;

    /*
     * Return the smallest r1>=r such that d_{r1} has a possible target
     *
     * Range >= t_max is deem unknown and thus possiple.
     * Return R_PERM if not found
     */
    int NextRTgt(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, int r) const;

    /*
     * Return the largest r1<=r such that d_{r1} has a possible source
     *
     * Return -1 if not found
     */
    int NextRSrc(const Staircases1d& nodes_ss, AdamsDeg deg, int r) const;

    int1d GetDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, int r) const;
    bool IsNewDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const;

private:
    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    void SetDiffSc(std::string_view name, Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r);

    /* Add an image. dx must be nonempty. */
    void SetImageSc(std::string_view name, Staircases1d& nodes_ss, AdamsDeg deg_dx, const int1d& dx, const int1d& x, int r);

    /**
     * Add d_r(x)=dx;
     * Add d_s(xy)=d_s(x)y+xd_s(y) for y with level=LEVEL_MAX-s and s>=r_min.
     * Return the number of changed degrees.
     *
     * dx should not be null.
     */
    int SetRingDiffLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, bool bFastTry = false);
    int SetModuleDiffLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, bool bFastTry = false);

public:
    /* Add a node */
    void AddNode(DeduceFlag flag);

    /* Pop the lastest node */
    void PopNode(DeduceFlag flag);

    /* Apply the change of the staircase to the current history */
    void UpdateStaircase(Staircases1d& nodes_ss, AdamsDeg deg, const Staircase& sc_i, size_t i_insert, const int1d& x, const int1d& dx, int level, int1d& image, int& level_image);

    /* Cache null diffs to the most recent node. */
    void CacheNullDiffs(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, DeduceFlag flag, NullDiff1d& nds);

    /**
     * Check first if it is a new differential before adding it.
     * Do some deductions by degree.
     *
     * Return the number of changed degrees.
     */
    int SetRingDiffGlobal(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry = false);
    int SetModuleDiffGlobal(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry = false);
    int SetCwDiffGlobal(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool bFastTry = false);

    /* Add d_r(?)=x;
     * Add d_r(?)=xy for d_r(y)=0 (y on level < LEVEL_MAX - r);
     */
    int SetRingBoundaryLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, int r);
    int SetModuleBoundaryLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, int r);

public: /* Differentials */
    int DeduceTrivialDiffs();
    /* Return 0 if there is no exception */
    int TryDiff(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int depth, DeduceFlag flag);
    int DeduceDiffs(size_t iCw, AdamsDeg deg, int depth, DeduceFlag flag);
    int DeduceDiffs(int stem_min, int stem_max, int depth, DeduceFlag flag);

public:
    /* Return if Einf at deg is possibly nontrivial */
    static int PossEinf(const Staircases1d& nodes_ss, AdamsDeg deg);
    static void UpdatePossEinf(const Staircases1d& nodes_ss, ut::map_seq2d<int, 0>& basis_ss_possEinf);
    void UpdateAllPossEinf()
    {
        for (auto& ring : rings_)
            UpdatePossEinf(ring.nodes_ss, ring.basis_ss_possEinf);
        for (auto& mod : modules_)
            UpdatePossEinf(mod.nodes_ss, mod.basis_ss_possEinf);
    }
    /* Return if Einf at deg could have more nontrivial elements */
    static int PossMoreEinf(const Staircases1d& nodes_ss, AdamsDeg deg);
    void PossMoreEinfFirstS_Ring(size_t iRing, int1d& O1s, int1d& O2s, int1d& isSingle) const;
    void PossMoreEinfFirstS_Mod(size_t iMod, int1d& O1s, int1d& O2s, int1d& isSingle) const;

    /*
     * Return the smallest s1>=s such that extension has a possible target
     *
     * Range >= t_max is deem unknown and thus possiple.
     * Return FIL_MAX if not found
     */
    int ExtendRelRing(size_t iRing, int stem, const algZ::Poly& rel, algZ::Poly& rel_extended) const;
    int ExtendRelMod(size_t iCof, int stem, const algZ::Mod& rel, algZ::Mod& rel_extended) const;
    int ExtendRelRing(size_t iRing, int stem, algZ::Poly& rel) const;
    int ExtendRelMod(size_t iCof, int stem, algZ::Mod& rel) const;
    int ExtendRelRingV2(size_t iRing, int stem, algZ::Poly& rel, ut::map_seq<int, 0>& num_leads) const;
    int ExtendRelCofV2(size_t iCof, int stem, algZ::Mod& rel, ut::map_seq<int, 0>& num_leads) const;
    int2d GetRingGbEinf(size_t iRing, AdamsDeg deg) const;
    std::map<AdamsDeg, int2d> GetRingGbEinf(size_t iRing) const;
    int2d GetModuleGbEinf(size_t iMod, AdamsDeg deg) const;
    std::map<AdamsDeg, int2d> GetModuleGbEinf(size_t iMod) const;

public: /* homotopy groups */
    void SetPermanentCycle(int depth, size_t iCof, AdamsDeg deg_x);
    void AddPiRelsRing(size_t iRing, algZ::Poly1d rels);
    void AddPiRelsCof(size_t iMod, algZ::Mod1d rels);
    void AddPiRelsByNat(size_t iMod);
    void SimplifyPiRels()
    {
        for (auto& ring : rings_)
            ring.pi_gb.SimplifyRels(ring.basis_ss_possEinf);
        for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
            modules_[iCof].pi_gb.SimplifyRels(modules_[iCof].basis_ss_possEinf);
    }
    static algZ::Mon1d GenBasis(const algZ::Groebner& gb, AdamsDeg deg, const PiBasis1d& nodes_pi_basis);
    static algZ::MMod1d GenBasis(const algZ::GroebnerMod& gb, AdamsDeg deg, const PiBasis1d& basis);
    void SyncS0Homotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopyy, int depth);
    void SyncCofHomotopy(int iCof, AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth);
    void SyncHomotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth)
    {
        SyncS0Homotopy(deg_min, count_ss, count_homotopy, depth);
        for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
            SyncCofHomotopy((int)iCof, deg_min, count_ss, count_homotopy, depth);
    }
    int DeduceTrivialExtensions(int depth);
    int DeduceExtensions2tor();
    int DeduceExtensionsByExactness(int stem_min, int stem_max, int depth);

    unsigned TryExtS0(algZ::Poly rel, AdamsDeg deg_change, int depth, DeduceFlag flag);
    unsigned TryExtCof(size_t iCof, algZ::Mod rel, AdamsDeg deg_change, int depth, DeduceFlag flag);
    unsigned TryExtQ(size_t iCof, size_t gen_id, algZ::Poly q, AdamsDeg deg_change, int depth, DeduceFlag flag);
    void DeduceExtensions(int stem_min, int stem_max, int& count_ss, int& count_homotopy, int depth, DeduceFlag flag);
    int DefineDependenceInExtensions(int stem_min, int stem_max, int depth);
    int DefineDependenceInExtensionsV2(int stem_min, int stem_max, int stem_max_mult, int depth);
};

class DBSS : public myio::DbAdamsSS
{
    using Statement = myio::Statement;

public:
    DBSS() = default;
    explicit DBSS(const std::string& filename) : DbAdamsSS(filename) {}

    void create_basis_ss(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss (id INTEGER PRIMARY KEY, base TEXT, diff TEXT, level SMALLINT, s SMALLINT, t SMALLINT)");
    }

    void create_pi_generators_mod(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, Einf TEXT, to_S0 TEXT, s SMALLINT, t SMALLINT)");
    }

    void create_pi_def(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_generators_def (id INTEGER PRIMARY KEY, def TINYINT, map_ids TEXT, multipliers TEXT, mult_name TEXT, fils TEXT)");
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

    void drop_and_create_pi_definitions(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_pi_generators_def");
        create_pi_def(table_prefix);
    }

    void save_pi_generators_mod(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Mod1d& gen_Einf) const;
    void save_basis_ss(const std::string& table_prefix, const Staircases& nodes_ss) const;
    void save_pi_basis(const std::string& table_prefix, const PiBasis& basis) const;
    void save_pi_basis_mod(const std::string& table_prefix, const PiBasisMod& basis) const;
    void save_pi_def(const std::string& table_prefix, const std::vector<EnumDef>& pi_gen_defs, const std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const;
    /* load the minimum id in every degree */
    std::map<AdamsDeg, int> load_basis_indices(const std::string& table_prefix) const;
    void update_basis_ss(const std::string& table_prefix, const std::map<AdamsDeg, Staircase>& nodes_ss) const;
    Staircases load_basis_ss(const std::string& table_prefix) const;
    void load_pi_def(const std::string& table_prefix, std::vector<EnumDef>& pi_gen_defs, std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const;
};

std::ostream& operator<<(std::ostream& sout, const int1d& arr);

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
void GetAllDbNames(const std::string& diagram_name, std::vector<std::string>& names, std::vector<std::string>& paths, std::vector<int>& isRing, bool log=false);

int main_deduce(int argc, char** argv, int index);
int main_deduce_ext(int argc, char** argv, int index);
int main_deduce_ext_def(int argc, char** argv, int index);
int main_deduce_ext_def2(int argc, char** argv, int index);
int main_deduce_ext_2tor(int argc, char** argv, int index);

#endif
