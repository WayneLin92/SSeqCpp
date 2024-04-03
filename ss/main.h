#ifndef MAIN_H
#define MAIN_H

#include "algebras/benchmark.h"
#include "algebras/dbAdamsSS.h"
#include "json.h"
#include "pigroebner.h"
#include <memory>
#include <set>
#include <variant>

inline const char* PROGRAM = "ss";
inline const char* VERSION = "Version:\n  2.1 (2024-03-25)";
using namespace alg2;

constexpr int LEVEL_MAX = 10000;
constexpr int LEVEL_MIN = 2;
constexpr int R_PERM = 1000;
constexpr int LEVEL_PERM = LEVEL_MAX - R_PERM; /* Level of Permanant cycles */

constexpr size_t MAX_DEPTH = 3; /* Maximum deduction depth */

inline const auto NULL_DIFF = int1d{-1};
inline const algZ::Mod MOD_V0 = algZ::MMod(algZ::Mon(), 0, 0);
using size_t1d = std::vector<size_t>;

enum class SSFlag : uint32_t
{
    no_op = 0,
    all_x = 1,                /* Deduce dx for all x including linear combinations */
    xy = 2,                   /* Deduce d(xy) for even if dx is uncertain */
    cofseq = 4,               /* Consider cofiber sequence */
    pi = 8,                   /* Consider homotopy */
    pi_def = 16,              /* Define generators in pi */
    no_save = 32,             /* Do not save the database */
    depth_ss_cofseq = 4 + 64, /* Deduce cof inside ss */
    depth_ss_ss = 128,        /* Deduce cof inside ss */
    naming = 256,             /* naming mode */
    try_all = 512,            /* Try all possible differentials */
};

inline SSFlag operator|(SSFlag lhs, SSFlag rhs)
{
    return SSFlag(uint32_t(lhs) | uint32_t(rhs));
}

inline bool operator&(SSFlag lhs, SSFlag rhs)
{
    return uint32_t(lhs) & uint32_t(rhs);
}

enum class EnumDef : int
{
    no_def = 0,
    dec = 1,         /* Decomposable */
    constraints = 2, /* The indeterminancies have some constraints */
};

/*
 * [0]=d2[0]
 * [1]=d2[?]
 * [2]=d3[1]
 * [3]=d3[?]
 * ...
 * d2[7]=[7]
 * d2[8]=[?]
 */
struct Staircase
{
    int2d basis;
    int2d diffs; /* element {-1} means null */
    int1d levels;
};
using Staircases = std::map<AdamsDeg, Staircase>;
using Staircases1d = std::vector<Staircases>;

template <>
struct fmt::formatter<Staircase>
{
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const Staircase& sc, FormatContext& ctx)
    {
        for (size_t i = 0; i < sc.levels.size(); ++i) {
            if (sc.levels[i] < LEVEL_MAX / 2) {
                if (sc.diffs[i] == NULL_DIFF)
                    fmt::format_to(ctx.out(), "{}=d_{}[?]\n", sc.basis[i], sc.levels[i]);
                else
                    fmt::format_to(ctx.out(), "{}=d_{}{}\n", sc.basis[i], sc.levels[i], sc.diffs[i]);
            }
            else {
                if (sc.diffs[i] == NULL_DIFF)
                    fmt::format_to(ctx.out(), "d_{}{}=[?]\n", LEVEL_MAX - sc.levels[i], sc.basis[i]);
                else
                    fmt::format_to(ctx.out(), "d_{}{}={}\n", LEVEL_MAX - sc.levels[i], sc.basis[i], sc.diffs[i]);
            }
        }
        return ctx.out();
    }
};

/* cw1 --map1--> cw2 --map2--> cw3 --map3--> cw1 */
struct CofSeq
{
    std::string name;
    std::array<size_t, 3> indexMap;
    std::array<AdamsDeg, 3> degMap;
    std::array<bool, 3> isRing;
    std::array<size_t, 3> indexCw;
    std::array<std::string, 3> nameCw;
    std::array<int, 3> t_max;
    std::array<Staircases1d*, 3> nodes_ss;
    std::array<Staircases1d, 3> nodes_cofseq; /* size = depth + 2 */
};
using CofSeq1d = std::vector<CofSeq>;

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

class InteruptAndSaveException : public MyException
{
public:
    InteruptAndSaveException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
    InteruptAndSaveException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

struct NullDiff
{
    int1d x;
    int r;
    int direction = 0;
    int first, count;
};
using NullDiff1d = std::vector<NullDiff>;

struct NullDiffCofseq
{
    int1d x;
    int r;
    int direction = 0;
    int first, count;
    int first_ss, count_ss;
};
using NullDiffCofseq1d = std::vector<NullDiffCofseq>;

/* This constrains the choice of a generator g by requiring that Unknownfil(g*m)=O. */
struct GenConstraint
{
    int map_index;
    algZ::Mon m;
    int O;
};

struct IndexCof
{
    int iCof = -1, iCs = -1;
};

struct RingSp
{
    /* #metadata */
    std::string name;
    int t_max = -1;
    std::vector<size_t> ind_mods, ind_maps;
    std::vector<IndexCof> ind_cofs;

    /* #ss */
    Groebner gb;
    std::map<AdamsDeg, Mon1d> basis;
    AdamsDeg1d degs_basis_order_by_stem;      /* Constant after initialization */
    Staircases1d nodes_ss;                    /* size = depth + 2 */
    ut::map_seq2d<int, 0> basis_ss_possEinf;  //// TODO: change to int2d

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
    int t_max = -1;
    size_t iRing;
    std::vector<size_t> ind_maps;
    std::vector<IndexCof> ind_cofs;

    /* #ss */
    GroebnerMod gb;
    std::map<AdamsDeg, MMod1d> basis;
    AdamsDeg1d degs_basis_order_by_stem;      /* Constant after initialization */
    Staircases1d nodes_ss;                    /* size = depth + 2 */
    ut::map_seq2d<int, 0> basis_ss_possEinf;  //// TODO: change to int2d

    /* #pi */
    algZ::GroebnerMod pi_gb;
    Mod1d pi_gen_Einf;
    PiBasisMod1d nodes_pi_basis = {{}}; /* size = depth + 1 */

    std::vector<EnumDef> pi_gen_defs;
    std::vector<std::vector<GenConstraint>> pi_gen_def_mons;
};
using ModSp1d = std::vector<ModSp>;

class Diagram;

class Map
{
public:
    std::string name, display;
    int t_max = -1;
    AdamsDeg deg;
    IndexCof ind_cof;

public:
    Map() {}
    Map(std::string name, std::string display, int t_max, AdamsDeg deg) : name(std::move(name)), display(std::move(display)), t_max(t_max), deg(deg) {}
    virtual ~Map() {}
    virtual bool IsFromRing(size_t& from) const = 0;
    virtual bool IsToRing(size_t& to) const = 0;
    virtual bool IsMul() /* in maps_v2 */
    {
        return false;
    }
    virtual int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const = 0;
};
using PMap1d = std::vector<std::unique_ptr<Map>>;

class MapRing2Ring : public Map  ////TODO: Add pi_images
{
public:
    size_t from, to;
    std::vector<Poly> images;
    MapRing2Ring(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, std::vector<Poly> images) : Map(std::move(name), std::move(display), t_max, deg), from(from), to(to), images(std::move(images)) {}
    bool IsFromRing(size_t& from_) const
    {
        from_ = from;
        return true;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = to;
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
};

class MapMod2Mod : public Map
{
public:
    size_t from, to;
    std::vector<Mod> images;
    MapMod2Mod(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, std::vector<Mod> images) : Map(std::move(name), std::move(display), t_max, deg), from(from), to(to), images(std::move(images)) {}
    bool IsFromRing(size_t& from_) const
    {
        from_ = from;
        return false;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = to;
        return false;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
    void Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs);
};

class MapMod2ModV2 : public Map
{
public:
    size_t from, to, over;
    std::vector<Mod> images;
    MapMod2ModV2(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, size_t over, std::vector<Mod> images)
        : Map(std::move(name), std::move(display), t_max, deg), from(from), to(to), over(over), images(std::move(images))
    {
    }
    bool IsFromRing(size_t& from_) const
    {
        from_ = from;
        return false;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = to;
        return false;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
    void Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs);
};

class MapMod2Ring : public Map
{
public:
    size_t from, to;
    std::vector<Poly> images;
    MapMod2Ring(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, std::vector<Poly> images) : Map(std::move(name), std::move(display), t_max, deg), from(from), to(to), images(std::move(images)) {}
    bool IsFromRing(size_t& from_) const
    {
        from_ = from;
        return false;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = to;
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
    void Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs);
};

class MapMod2RingV2 : public Map
{
public:
    size_t from, to, over;
    std::vector<Poly> images;
    MapMod2RingV2(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, size_t over, std::vector<Poly> images)
        : Map(std::move(name), std::move(display), t_max, deg), from(from), to(to), over(over), images(std::move(images))
    {
    }
    bool IsFromRing(size_t& from_) const
    {
        from_ = from;
        return false;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = to;
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
    void Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs);
};

class MapMulRing2Ring : public Map
{
public:
    size_t index;
    Poly factor;
    MapMulRing2Ring(std::string name, std::string display, int t_max, AdamsDeg deg, size_t index, Poly factor) : Map(std::move(name), std::move(display), t_max, deg), index(index), factor(std::move(factor)) {}
    bool IsFromRing(size_t& from_) const
    {
        from_ = index;
        return true;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = index;
        return true;
    }
    bool IsMul()
    {
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
};

class MapMulRing2Mod : public Map
{
public:
    size_t from, to;
    Mod factor;
    MapMulRing2Mod(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, Mod factor) : Map(std::move(name), std::move(display), t_max, deg), from(from), to(to), factor(std::move(factor)) {}
    bool IsFromRing(size_t& from_) const
    {
        from_ = from;
        return true;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = to;
        return false;
    }
    bool IsMul()
    {
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
};

class MapMulMod2Mod : public Map
{
public:
    size_t index;
    Poly factor;
    MapMulMod2Mod(std::string name, std::string display, int t_max, AdamsDeg deg, size_t index, Poly factor) : Map(std::move(name), std::move(display), t_max, deg), index(index), factor(std::move(factor)) {}
    bool IsFromRing(size_t& from_) const
    {
        from_ = index;
        return false;
    }
    bool IsToRing(size_t& to_) const
    {
        to_ = index;
        return false;
    }
    bool IsMul()
    {
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const;
};

struct Commutativity
{
    std::string name;
    int f0, f1, g0, g1;
};
using Commutativity1d = std::vector<Commutativity>;

constexpr size_t FLAG_MOD = 1 << 16;

class Diagram
{
protected:
    RingSp1d rings_;
    ModSp1d modules_;
    PMap1d maps_;
    CofSeq1d cofseqs_;
    Commutativity1d comms_;

protected:
    nlohmann::json js_;
    std::vector<size_t> deduce_list_spectra_;
    std::vector<size_t> deduce_list_cofseq_;
    int depth_ = 0;
    int deduce_count_max_ = 10;
    AdamsDeg deg_leibniz_;             /* For logging */
    const int1d* a_leibniz_ = nullptr; /* For logging */

public:
    Diagram(std::string diagram_name, SSFlag flag, bool log = true, bool loadD2 = false);
    void VersionConvertReorderRels()
    {
        for (auto& ring : rings_)
            ring.pi_gb.SimplifyRelsReorder(ring.basis_ss_possEinf);
        for (auto& mod : modules_)
            mod.pi_gb.SimplifyRelsReorder(mod.basis_ss_possEinf);
    }
    void save(std::string diagram_name, SSFlag flag);

public: /* Getters */
    static auto GetRecentPiBasis(const PiBasis1d& nodes_pi_basis, AdamsDeg deg) -> const PiBase*;
    static auto GetRecentPiBasis(const PiBasisMod1d& nodes_pi_basis, AdamsDeg deg) -> const PiBaseMod*;

    auto& GetJs() const
    {
        return js_;
    }

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

    auto& GetMaps()
    {
        return maps_;
    }

    auto& GetMaps() const
    {
        return maps_;
    }

    auto& GetCofSeqs()
    {
        return cofseqs_;
    }

    auto& GetCofSeqs() const
    {
        return cofseqs_;
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

    int GetMapIndexByName(const std::string& name) const
    {
        for (size_t i = 0; i < maps_.size(); ++i)
            if (maps_[i]->name == name)
                return (int)i;
        return -1;
    }

    int GetCofSeqIndexByName(const std::string& name) const
    {
        for (size_t i = 0; i < cofseqs_.size(); ++i)
            if (cofseqs_[i].name == name)
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
     * Return a level such that nontrivial diffs in >= level will not be longer.
     */
    static int GetFirstFixedLevelForPlot(const Staircases1d& nodes_ss, AdamsDeg deg);
    static int GetFirstFixedLevelForPlotCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg);

public:
    void SetDeduceList(const std::vector<std::string>& cws);

private: /* Staircase */
    static size_t GetFirstIndexOfNullOnLevel(const Staircase& sc, int level);
    static int GetMaxLevelWithND(const Staircase& sc);
    static bool IsZeroOnLevel(const Staircase& sc, const int1d& x, int level);

private: /* ss */
    /* Warning: The following IsPoss functions do not check if ss[deg] exists */
    /* Check if deg can be hit by dr for r<=r_max */
    static bool IsPossTgt(const Staircases1d& nodes_ss, AdamsDeg deg, int r_max);
    static bool IsPossTgtCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r_max);
    static bool IsPossTgt(const Staircases1d& nodes_ss, AdamsDeg deg)
    {
        return IsPossTgt(nodes_ss, deg, R_PERM);
    }

    /* Return the first index with level >=`level_min` such that all levels above are already fixed */
    static size_t GetFirstIndexOfFixedLevels(const Staircases1d& nodes_ss, AdamsDeg deg, int level_min);
    static size_t GetFirstIndexOfFixedLevelsCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int level_min);

    /* Count the number of all possible d_r targets. Return (index, count). */
    auto CountPossDrTgt(const Staircases1d& nodes_ss, int t_max, const AdamsDeg& deg_tgt, int r) const -> std::pair<int, int>;
    /* Warning: this assumes that there shall not be more Einf elements */
    auto CountPossDrTgtCofseq(const CofSeq& cofseq, size_t iCs, const AdamsDeg& deg_tgt, int r) const -> std::pair<int, int>;

    /* Count the number of all possible d_r sources. Return (index, count). */
    std::pair<int, int> CountPossDrSrc(const Staircases1d& nodes_ss, const AdamsDeg& deg_src, int r) const;
    /* Count the number of all potential permant cycles. Return (index, count). */
    auto CountPossMorePerm(const Staircases1d& nodes_ss, const AdamsDeg& deg) const -> std::pair<int, int>
    {
        return CountPossDrSrc(nodes_ss, deg, R_PERM - 1);
    };
    /* Warning: this assumes that there shall not be more Einf elements */
    std::pair<int, int> CountPossDrSrcCofseq(const CofSeq& cofseq, size_t iCs, const AdamsDeg& deg_src, int r) const;

    /*
     * Return the smallest r1>=r such that d_{r1} has a possible target
     *
     * Range >= t_max is unknown and thus always possible.
     * Return R_PERM if not found
     */
    int NextRTgt(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, int r) const;
    int NextRTgtCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r) const;

    /**
     * Return the largest r1<=r such that d_{r1} has a possible source
     *
     * Return -1 if not found
     */
    int NextRSrc(const Staircases1d& nodes_ss, AdamsDeg deg, int r) const;
    int NextRSrcCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r) const;

    int1d GetDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, int r) const;
    int1d GetLevelAndDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, int1d x, int& level) const;
    /* Return the minimal length of the cross differentials */
    int GetCofseqCrossR(const Staircases1d& nodes_cofseq, const Staircases1d& nodes_ss, AdamsDeg deg, int t_max, int r_min) const;

private:
    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    void SetDiffSc(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag);
    /* Add an image. dx must be nonempty and not null. x must be nonempty. */
    void SetImageSc(size_t iCw, AdamsDeg deg_dx, const int1d& dx, const int1d& x, int r, SSFlag flag);

    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    void SetDiffScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x_, const int1d& dx, int r, SSFlag flag);
    /* Add an image. dx must be nonempty and not null. x must be nonempty. */
    void SetImageScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_dx, const int1d& dx_, const int1d& x, int r, SSFlag flag);
    /* Retriangulate when ss changes */
    void ReSetScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg, SSFlag flag);

    /**
     * Add d_r(x)=dx;
     * Add d_s(xy)=d_s(x)y+xd_s(y) for y with level=LEVEL_MAX-s and s>=r_min.
     * Return the number of changed degrees.
     *
     * dx should not be null.
     */
    int SetRingDiffLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag);
    /** This version deduces d(xy) even if dx is uncertain */
    int SetRingDiffLeibnizV2(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);
    int SetModuleDiffLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag);
    /** This version deduces d(xy) even if dx is uncertain */
    int SetModuleDiffLeibnizV2(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);

public:
    /* Add a node */
    void AddNode(SSFlag flag);

    /* Pop the lastest node */
    void PopNode(SSFlag flag);

    /* Apply the change of the staircase to the current history */
    void UpdateStaircase(Staircases1d& nodes_ss, AdamsDeg deg, const Staircase& sc_i, size_t i_insert, const int1d& x, const int1d& dx, int level, int1d& image, int& level_image);

    /* Cache null diffs to the most recent node. */
    void CacheNullDiffs(const Staircases1d& nodes_ss, int t_max, AdamsDeg deg, SSFlag flag, NullDiff1d& nds) const;
    void CacheNullDiffsCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, SSFlag flag, NullDiffCofseq1d& nds) const;

    /**
     * Check first if it is a new differential before adding it.
     * Do some deductions by degree.
     *
     * Return the number of changed degrees.
     */
    int SetRingDiffGlobal(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);
    int SetModuleDiffGlobal(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);
    int SetCwDiffGlobal(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);

    /* Add d_r(?)=x;
     * Add d_r(?)=xy for d_r(y)=0 (y on level < LEVEL_MAX - r);
     */
    int SetRingBoundaryLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);
    int SetModuleBoundaryLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);

    int SetDiffLeibnizCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag);
    int SetDiffGlobalCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);

public: /* Differentials */
    bool IsNewDiff(const Staircases1d& nodes_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const;
    bool IsNewDiffCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const;
    int DeduceTrivialDiffs(SSFlag flag);
    int DeduceTrivialDiffsCofseq(SSFlag flag);
    int DeduceManual(SSFlag flag);
    /* Return 0 if there is no exception */
    int TryDiff(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int depth, SSFlag flag, bool tryY);
    int TryDiffCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, AdamsDeg deg_dx, const int1d& x, const int1d& dx, const int1d& perm, int r, int depth, SSFlag flag, bool tryY);
    int DeduceDiffs(size_t iCw, AdamsDeg deg, int depth, SSFlag flag);
    int DeduceDiffs(const size_t1d& cws, int stem_min, int stem_max, int depth, SSFlag flag);
    int DeduceDiffs(int stem_min, int stem_max, int depth, SSFlag flag);
    /* Deduce d(xy) no matter what dx is */
    int DeduceDiffsV2();

    int DeduceDiffsCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg, int depth, SSFlag flag);
    int DeduceDiffsCofseq(int stem_min, int stem_max, int depth, SSFlag flag);
    int DeduceDiffsNbhdCofseq(CofSeq& cofseq, size_t iCs, int stem, int depth, SSFlag flag);
    int CommuteCofseq(SSFlag flag);
    void SyncCofseq(SSFlag flag);

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
    //void SetPermanentCycle(int depth, size_t iCof, AdamsDeg deg_x);
    void AddPiRelsRing(size_t iRing, algZ::Poly1d rels);
    void AddPiRelsCof(size_t iMod, algZ::Mod1d rels);
    //void AddPiRelsByNat(size_t iMod);
    void SimplifyPiRels()
    {
        for (auto& ring : rings_)
            ring.pi_gb.SimplifyRels();
        for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
            modules_[iCof].pi_gb.SimplifyRels();
    }
    static algZ::Mon1d GenBasis(const algZ::Groebner& gb, AdamsDeg deg, const PiBasis1d& nodes_pi_basis);
    static algZ::MMod1d GenBasis(const algZ::GroebnerMod& gb, AdamsDeg deg, const PiBasis1d& basis);
    //void SyncS0Homotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopyy, int depth);
    //void SyncCofHomotopy(int iCof, AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth);
    /*void SyncHomotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth)
    {
        SyncS0Homotopy(deg_min, count_ss, count_homotopy, depth);
        for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
            SyncCofHomotopy((int)iCof, deg_min, count_ss, count_homotopy, depth);
    }*/
    //int DeduceTrivialExtensions(int depth);
    //int DeduceExtensions2tor();
    //int DeduceExtensionsByExactness(int stem_min, int stem_max, int depth);

    //unsigned TryExtS0(algZ::Poly rel, AdamsDeg deg_change, int depth, SSFlag flag);
    //unsigned TryExtCof(size_t iCof, algZ::Mod rel, AdamsDeg deg_change, int depth, SSFlag flag);
    //unsigned TryExtQ(size_t iCof, size_t gen_id, algZ::Poly q, AdamsDeg deg_change, int depth, SSFlag flag);
    //void DeduceExtensions(int stem_min, int stem_max, int& count_ss, int& count_homotopy, int depth, SSFlag flag);
    //int DefineDependenceInExtensions(int stem_min, int stem_max, int depth);  ////;
    //int DefineDependenceInExtensionsV2(int stem_min, int stem_max, int stem_max_mult, int depth);
};

class DBSS : public myio::DbAdamsSS
{
    using Statement = myio::Statement;

public:
    DBSS() = default;
    explicit DBSS(const std::string& filename) : DbAdamsSS(filename) {}

    void create_ss(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_ss (id INTEGER PRIMARY KEY, s SMALLINT, t SMALLINT, base TEXT, diff TEXT, level SMALLINT)");
    }

    void create_cofseq(const std::string& table) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table + " (iC SMALLINT, s SMALLINT, t SMALLINT, base TEXT, diff TEXT, level SMALLINT)");
    }

    void create_pi_generators_mod(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, Einf TEXT, to_S0 TEXT, s SMALLINT, t SMALLINT)");
    }

    void create_pi_def(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_pi_generators_def (id INTEGER PRIMARY KEY, def TINYINT, map_ids TEXT, multipliers TEXT, mult_name TEXT, fils TEXT)");
    }

    void drop_and_create_ss(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss");
        create_ss(table_prefix);
    }

    void drop_and_create_cofseq_ss(const std::string& table_prefix) const
    {
        drop_table(table_prefix + "_ss");
        create_cofseq(table_prefix);
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

public:
    void update_ss(const std::string& table_prefix, const Staircases& nodes_ss) const;
    void save_ss(const std::string& table_prefix, const Staircases& nodes_ss) const;
    void save_cofseq(const std::string& table, const CofSeq& cofseq) const;
    void save_pi_generators_mod(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Mod1d& gen_Einf) const;
    void save_pi_basis(const std::string& table_prefix, const PiBasis& basis) const;
    void save_pi_basis_mod(const std::string& table_prefix, const PiBasisMod& basis) const;
    void save_pi_def(const std::string& table_prefix, const std::vector<EnumDef>& pi_gen_defs, const std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const;

public:
    /* load the minimum id in every degree */
    std::map<AdamsDeg, int> load_basis_indices(const std::string& table_prefix) const;
    Staircases load_ss(const std::string& table_prefix) const;
    std::array<Staircases, 3> load_cofseq(const std::string& table) const;
    void load_pi_def(const std::string& table_prefix, std::vector<EnumDef>& pi_gen_defs, std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const;
};

inline std::optional<std::string> SerializeDiff(const int1d& dx)
{
    std::optional<std::string> result;
    if (dx == NULL_DIFF)
        return result;
    result = myio::Serialize(dx);
    return result;
}
std::ostream& operator<<(std::ostream& sout, const int1d& arr);

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

inline bool BelowS0VanishingLine(AdamsDeg deg)
{
    return 3 * deg.s <= deg.t + 3;
}

/* Strictly above the vanishing line */
inline bool AboveS0Vanishing(AdamsDeg deg)
{
    return 3 * (deg.s - 1) > deg.t;
}

/* Strictly above the vanishing line */
inline bool AboveJ(AdamsDeg deg)
{
    return 3 * (deg.s + 1) >= deg.t;
}

/* Strictly above the vanishing line */
inline bool BelowCokerJ(AdamsDeg deg)
{
    return 5 * deg.s <= deg.stem() && deg.stem() >= 10;
}

size_t GetFirstIndexOnLevel(const Staircase& sc, int level);

/* Compute x mod (level-1) */
int1d Residue(int1d x, const Staircases1d& nodes_ss, AdamsDeg deg, int level);

void GetAllDbNames(const std::string& diagram_name, std::vector<std::string>& names, std::vector<std::string>& paths, std::vector<int>& isRing, bool log = false);

#endif
