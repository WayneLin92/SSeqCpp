#ifndef MAIN_H
#define MAIN_H

#include "algebras/benchmark.h"
#include "algebras/dbAdamsSS.h"
#include "algebras/linalg.h"
#include "json.h"
#include "pigroebner.h"
#include <memory>
#include <ranges>
#include <set>
#include <variant>

inline const char* PROGRAM = "ss";
inline const char* VERSION = "1.0.1";
using namespace alg2;
using size1d = std::vector<size_t>;
using json = nlohmann::json;

inline constexpr int LEVEL_MAX = 10000;
inline constexpr int LEVEL_MIN = 2;
inline constexpr int R_PERM = 1000;
inline constexpr int LEVEL_PERM = LEVEL_MAX - R_PERM; /* Level of Permanant cycles */

inline const auto NULL_DIFF = int1d{-1};
inline const algZ::Mod MOD_V0 = algZ::MMod(algZ::Mon(), 0, 0);
inline constexpr size_t MAX_DEPTH = 3; /* Maximum deduction depth */

enum class SSFlag : uint32_t
{
    no_op = 0,
    cofseq = 1,              /* Consider cofiber sequences */
    pi = 2,                  /* Consider homotopy */
    pi_def = 4,              /* Define generators in pi */
    no_save = 8,             /* Do not save the database */
    deduce_4_all_x = 16,     /* Deduce dx for all x including linear combinations */
    deduce_dxy = 32,         /* Deduce d_r(xy) even if dx is uncertain */
    deduce_zero = 64,        /* Expand zero d_r differentials */
    deduce_pullback = 128,   /* Consider pullbacks upon assumtion */
    naming = 256,            /* Naming mode */
    try_all = 512,           /* Try all possible differentials */
    synthetic = 1024,        /* Sync method */
    no_exclusions = 1 << 11, /* Do not save the database */
    stacked = 1 << 12,       /* In a SetDiffGlobal function called in a SetDiffGlobal function */
    log_proof = 1 << 13,     /* Log proof */
    log_deg = 1 << 14,       /* Log for degree reason */
    log_nat = 1 << 15,       /* Log naturality */
};

inline SSFlag operator|(SSFlag lhs, SSFlag rhs)
{
    return SSFlag(uint32_t(lhs) | uint32_t(rhs));
}

inline bool operator&(SSFlag lhs, SSFlag rhs)
{
    return uint32_t(lhs) & uint32_t(rhs);
}

enum class CrossType
{
    no_cross,
    no_strict_cross,
    all_cases
};

enum class EnumDef : int
{
    no_def = 0,
    dec = 1,         /* Decomposable */
    constraints = 2, /* The indeterminancies have some constraints */
};

/* A return type for the logical state */
class [[nodiscard]] SSRet
{
public:
    unsigned code = 0;
    std::string err_msg;

public:
    SSRet() = default;
    explicit SSRet(int code) : code(code) {}
    static SSRet NUL()
    {
        return SSRet(0);
    }
    static SSRet CHANGE()
    {
        return SSRet(1);
    }
    static SSRet CHANGE_AT_X()
    {
        return SSRet(2);
    }
    static SSRet FAIL()
    {
        return SSRet(10);
    }
    static SSRet FAIL_SS()
    {
        return SSRet(11);
    }
    static SSRet FAIL_SS_SS()
    {
        return SSRet(100);
    }
    /* Return if this is a failure code */
    explicit operator bool() const
    {
        return code >= FAIL().code;
    }
    bool operator!() const
    {
        return !bool(*this);
    }
    SSRet& operator+=(SSRet rhs)
    {
        code = std::max(code, rhs.code);
        if (rhs.err_msg.size()) {
            if (err_msg.empty())
                std::swap(err_msg, rhs.err_msg);
            else
                err_msg += rhs.err_msg;
        }
        return *this;
    }

    bool IsChanged() const
    {
        return code >= CHANGE().code;
    }

    bool IsChangedAtX() const
    {
        return code >= CHANGE_AT_X().code;
    }
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

    bool non_empty() const
    {
        return levels.size();
    }

    std::string Str() const;
};

inline const auto EmptyStaircase = Staircase{};

class SS
{
protected:
    std::vector<std::vector<Staircase>> data_;
    /* The first bit marks the existance of data.
     * The second bit marks that the data is synced.
     * The third bit marks that the data is old
     * A new data is always marked with 1.
     */
    std::vector<std::vector<uint32_t>> valid_;

protected:
    /* ordered by (stem, s) */
    auto ijs() const
    {
        return ut::Iter2d(&valid_) | std::views::filter([this](ut::IJ ij) { return valid_[ij.first][ij.second]; });
    }

public:
    size_t size() const
    {
        return valid_.size();
    }
    void resize(size_t size)
    {
        data_.resize(size);
        valid_.resize(size);
    }
    auto& at(AdamsDeg deg) const
    {
        return data_[deg.stem()][deg.s];
    }
    auto& operator[](AdamsDeg deg)
    {
        ut::get(ut::get(valid_, deg.stem()), deg.s) = 1;
        return ut::get(ut::get(data_, deg.stem()), deg.s);
    }
    bool has(AdamsDeg deg) const
    {
        int n = deg.stem();
        return n < data_.size() && deg.s < data_[n].size() && valid_[n][deg.s];
    }
    bool Sync(AdamsDeg deg)
    {
        bool result = valid_[deg.stem()][deg.s] & 2;
        valid_[deg.stem()][deg.s] |= 2;
        return result;
    }
    void MarkModified(AdamsDeg deg)
    {
        valid_[deg.stem()][deg.s] = 1;
    }
    bool empty() const
    {
        return data_.empty();
    }
    int t_max() const
    {
        int t_max = -1;
        for (size_t n = 0; n < data_.size(); ++n) {
            int t = (int)(n + data_[n].size() - 1);
            t_max = std::max(t_max, t);
        }
        return t_max;
    }
    /* ordered by (stem, s) */
    auto degs() const
    {
        return ijs() | std::views::transform([](ut::IJ ij) { return AdamsDeg((int)ij.second, int(ij.first + ij.second)); });
    }
    /* ordered by (stem, s) */
    auto unsynced_degs()
    {
        return ut::Iter2d(&valid_) | std::views::filter([this](ut::IJ ij) {
                   if (valid_[ij.first][ij.second]) {
                       if (valid_[ij.first][ij.second] & 2)
                           return false;
                       valid_[ij.first][ij.second] |= 2;
                       return true;
                   }
                   return false;
               })
               | std::views::transform([](ut::IJ ij) { return AdamsDeg((int)ij.second, int(ij.first + ij.second)); });
    }
    /* ordered by (stem, s) */
    auto items() const
    {
        return ijs() | std::views::transform([this](ut::IJ ij) {
                   auto d = AdamsDeg((int)ij.second, int(ij.first + ij.second));
                   return std::pair<AdamsDeg, const Staircase&>(d, at(d));
               });
    }

    /* ordered by (t, s) */
    AdamsDeg1d arr_degs() const
    {
        AdamsDeg1d result;
        for (AdamsDeg d : degs())
            result.push_back(d);
        std::sort(result.begin(), result.end());
        return result;
    }
};

/* incorrect differentials that can be ruled out */
struct Exclusion
{
    int1d x;
    int r = 10086;
    int2d dxs;
};
using Exclusion1d = std::vector<Exclusion>;
class Exclusions
{
protected:
    std::vector<std::vector<Exclusion1d>> data_;
    std::vector<std::vector<uint32_t>> valid_;

public:
    auto& get(AdamsDeg deg, const int1d& x, int r)
    {
        ut::get(ut::get(valid_, deg.stem()), deg.s) = 1;
        auto& e1d = ut::get(ut::get(data_, deg.stem()), deg.s);
        for (auto& e : e1d)
            if (e.x == x && e.r == r)
                return e;
        e1d.push_back(Exclusion{x, r, {}}); /* Add a new exclusion */
        return e1d.back();
    }
    Exclusion* has(AdamsDeg deg, const int1d& x, int r)
    {
        int n = deg.stem();
        if (n < data_.size() && deg.s < data_[n].size() && valid_[n][deg.s]) {
            auto& e1d = data_[n][deg.s];
            for (auto& e : e1d)
                if (e.x == x && e.r == r)
                    return &e;
        }
        return nullptr;
    }
    void Sort()
    {
        for (auto p : ut::Iter2d(&data_))
            if (valid_[p.first][p.second])
                for (auto& e : data_[p.first][p.second])
                    std::sort(e.dxs.begin(), e.dxs.end());
    }
};

class SSNodes;
Staircase GetSc4Display(const SSNodes& nodes_this, const SSNodes& nodes_src, AdamsDeg deg, int stem_map);

class SSNodes : public ut::vector<SS, MAX_DEPTH + 2>
{
public:
    Exclusions exclusions;

public:
    const auto& GetRecentSc(AdamsDeg deg) const
    {
        for (auto p = rbegin(); p != rend(); ++p)
            if (p->has(deg))
                return p->at(deg);
        throw RunTimeError(fmt::format("Recent Value not found. deg={}.\n", deg));
    }

    /* ordered by (t, s) */
    auto unsynced_degs()
    {
        return back().unsynced_degs();
    }
    Staircase GetSc4Display(AdamsDeg deg) const
    {
        return ::GetSc4Display(*this, *this, deg, -1);
    }
    auto& changes() const
    {
        return data_[1];
    }
    void AddNode()
    {
        push_back({});
        back().resize(front().size());
    }
};

struct IndexUniv
{
    uint32_t type = uint32_t(-1), iTri = uint32_t(-1);
    size_t index = size_t(-1);

    IndexUniv() = default;
    IndexUniv(uint32_t type, size_t index) : type(type), index(index) {}
    IndexUniv(uint32_t type, uint32_t iTri, size_t index) : type(type), iTri(iTri), index(index) {}
    explicit operator bool() const
    {
        return type != uint32_t(-1);
    }
    bool isRing() const
    {
        return type == 0;
    }
    bool isModule() const
    {
        return type == 1;
    }
    bool isCofseq() const
    {
        return type == 2;
    }
    bool operator==(const IndexUniv& rhs) const = default;
    std::string Str() const
    {
        if (type == 0)
            return fmt::format("R{}", index);
        else if (type == 1)
            return fmt::format("M{}", index);
        else
            return fmt::format("Cs{}:{}", index, iTri);
    }
    IndexUniv next() const
    {
        return IndexUniv(type, (iTri + 1) % 3, index);
    }
    IndexUniv prev() const
    {
        return IndexUniv(type, (iTri + 2) % 3, index);
    }
};
inline IndexUniv IndexRing(size_t iRing)
{
    return IndexUniv(0, iRing);
}
inline IndexUniv IndexMod(size_t iMod)
{
    return IndexUniv(1, iMod);
}
inline IndexUniv IndexCof(size_t index, uint32_t iTri)
{
    return IndexUniv(2, iTri, index);
}
inline IndexUniv IndexCof(size_t index, size_t iTri)
{
    return IndexUniv(2, (uint32_t)iTri, index);
}

/* cw1 --map1--> cw2 --map2--> cw3 --map3--> cw1 */
struct CofSeq
{
    std::string name;
    std::array<size_t, 3> indexMap;
    std::array<AdamsDeg, 3> degMap;
    std::array<IndexUniv, 3> indexCw;
    std::array<std::string, 3> nameCw;
    std::array<int, 3> t_max;
    std::array<SSNodes*, 3> nodes_ss;
    std::array<SSNodes, 3> nodes_cofseq;
};
using CofSeq1d = std::vector<CofSeq>;

struct PiBase
{
    algZ::Mon1d nodes_pi_basis;
    int2d Einf;
};

// TODO: use vector instead of map
using PiBasis = std::map<AdamsDeg, PiBase>;
using PiBasis1d = std::vector<PiBasis>;

struct PiBaseMod
{
    algZ::MMod1d nodes_pi_basis;
    int2d Einf;
};

using PiBasisMod = std::map<AdamsDeg, PiBaseMod>;
using PiBasisMod1d = std::vector<PiBasisMod>;

/* Never to be thrown. A dummy Exception to be filled in the catch statement */
class NoException : public ErrorIdMsg
{
public:
    NoException(unsigned int error_id, const char* message) : ErrorIdMsg(error_id, message) {}
    NoException(unsigned int error_id, const std::string& message) : ErrorIdMsg(error_id, message) {}
};

class InteruptAndSaveException : public RunTimeError
{
public:
    InteruptAndSaveException(const char* message, const std::source_location location = std::source_location::current()) : RunTimeError(message, location)
    {
        fmt::print("Interupted\n");
    }
    InteruptAndSaveException(std::string message, const std::source_location location = std::source_location::current()) : RunTimeError(message, location)
    {
        fmt::print("Interupted\n");
    }
};

struct NullDiff
{
    int1d x;
    int r;
    int first, count;
};
using NullDiff1d = std::vector<NullDiff>;

struct NullDiffCofseq
{
    int1d x;
    int r;
    int first, count;
    int first_ss, count_ss;
    std::string Str() const;
};
using NullDiffCofseq1d = std::vector<NullDiffCofseq>;

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
    int t_max = -1;
    std::vector<size_t> ind_mods, ind_maps, ind_maps_prev;
    std::vector<IndexUniv> ind_cofs;

    /* #ss */
    AdamsDeg1d gen_degs;
    std::vector<std::string> gen_names;
    Groebner gb;
    BasisMon basis;
    AdamsDeg1d degs_ss;                       /* Constant after initialization */
    SSNodes nodes_ss;                         /* size = depth + 2 */
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
    std::vector<size_t> ind_maps, ind_maps_prev;
    std::vector<IndexUniv> ind_cofs;

    /* #ss */
    AdamsDeg1d v_degs;
    std::vector<std::string> v_names;
    GroebnerMod gb;
    BasisMMod basis;
    AdamsDeg1d degs_ss;                       /* Constant after initialization */
    SSNodes nodes_ss;                         /* size = depth + 2 */
    ut::map_seq2d<int, 0> basis_ss_possEinf;  //// TODO: change to int2d

    /* #pi */
    algZ::GroebnerMod pi_gb;
    Mod1d pi_gen_Einf;
    PiBasisMod1d nodes_pi_basis = {{}}; /* size = depth + 1 */

    std::vector<EnumDef> pi_gen_defs;
    std::vector<std::vector<GenConstraint>> pi_gen_def_mons;
};
using ModSp1d = std::vector<ModSp>;

class Category;

class Map
{
public:
    std::string name, display;
    int t_max = -1;
    AdamsDeg deg;
    IndexUniv from, to;
    IndexUniv iCof;

public:
    Map() {}
    Map(std::string name, std::string display, int t_max, AdamsDeg deg, IndexUniv from, IndexUniv to) : name(std::move(name)), display(std::move(display)), t_max(t_max), deg(deg), from(from), to(to) {}
    virtual ~Map() {}
    virtual bool IsMul() /* in maps_v2 */
    {
        return false;
    }
    virtual int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const = 0;
};
using PMap1d = std::vector<std::unique_ptr<Map>>;

class MapRing2Ring : public Map  ////TODO: Add pi_images
{
public:
    std::vector<Poly> images;
    MapRing2Ring(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, std::vector<Poly> images) : Map(std::move(name), std::move(display), t_max, deg, IndexRing(from), IndexRing(to)), images(std::move(images)) {}
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
};

class MapMod2Mod : public Map
{
public:
    std::vector<Mod> images;
    MapMod2Mod(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, std::vector<Mod> images) : Map(std::move(name), std::move(display), t_max, deg, IndexMod(from), IndexMod(to)), images(std::move(images)) {}
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
    void Verify(const Category& category, const AdamsDeg2d& ring_gen_degs);
};

class MapMod2ModV2 : public Map
{
public:
    size_t over;
    std::vector<Mod> images;
    MapMod2ModV2(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, size_t over, std::vector<Mod> images)
        : Map(std::move(name), std::move(display), t_max, deg, IndexMod(from), IndexMod(to)), over(over), images(std::move(images))
    {
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
    void Verify(const Category& category, const AdamsDeg2d& ring_gen_degs);
};

class MapMod2Ring : public Map
{
public:
    std::vector<Poly> images;
    MapMod2Ring(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, std::vector<Poly> images) : Map(std::move(name), std::move(display), t_max, deg, IndexMod(from), IndexRing(to)), images(std::move(images)) {}
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
    void Verify(const Category& category, const AdamsDeg2d& ring_gen_degs);
};

class MapMod2RingV2 : public Map
{
public:
    size_t over;
    std::vector<Poly> images;
    MapMod2RingV2(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, size_t over, std::vector<Poly> images)
        : Map(std::move(name), std::move(display), t_max, deg, IndexMod(from), IndexRing(to)), over(over), images(std::move(images))
    {
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
    void Verify(const Category& category, const AdamsDeg2d& ring_gen_degs);
};

class MapMulRing2Ring : public Map
{
public:
    Poly factor;
    MapMulRing2Ring(std::string name, std::string display, int t_max, AdamsDeg deg, size_t index, Poly factor) : Map(std::move(name), std::move(display), t_max, deg, IndexRing(index), IndexRing(index)), factor(std::move(factor)) {}
    bool IsMul()
    {
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
};

class MapMulRing2Mod : public Map
{
public:
    Mod factor;
    MapMulRing2Mod(std::string name, std::string display, int t_max, AdamsDeg deg, size_t from, size_t to, Mod factor) : Map(std::move(name), std::move(display), t_max, deg, IndexRing(from), IndexMod(to)), factor(std::move(factor)) {}
    bool IsMul()
    {
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
};

class MapMulMod2Mod : public Map
{
public:
    Poly factor;
    MapMulMod2Mod(std::string name, std::string display, int t_max, AdamsDeg deg, size_t index, Poly factor) : Map(std::move(name), std::move(display), t_max, deg, IndexMod(index), IndexMod(index)), factor(std::move(factor)) {}
    bool IsMul()
    {
        return true;
    }
    int1d map(const int1d& x, AdamsDeg deg_x, const Category& category) const;
};

struct Commutativity
{
    std::string name;
    int f0, f1, g0, g1;
};
using Commutativity1d = std::vector<Commutativity>;

struct ALeibnize
{
    std::string* pt_name = nullptr;
    AdamsDeg deg;
    int1d a;
};

class Category
{
protected:
    RingSp1d rings_;
    ModSp1d modules_;
    PMap1d maps_;
    CofSeq1d cofseqs_;
    Commutativity1d comms_;

protected:
    nlohmann::json js_;
    std::vector<IndexUniv> deduce_list_spectra_;
    std::vector<size_t> deduce_list_cofseq_;
    int depth_ = 0;
    int deduce_count_max_ = 10;

public:
    Category(const std::string& cat_root, const std::string& ckpt_name, SSFlag flag, bool log = true, bool loadD2 = false);
    void VersionConvertReorderRels()
    {
        for (auto& ring : rings_)
            ring.pi_gb.SimplifyRelsReorder(ring.basis_ss_possEinf);
        for (auto& mod : modules_)
            mod.pi_gb.SimplifyRelsReorder(mod.basis_ss_possEinf);
    }
    void LoadNodes(const std::string& cat_root, const std::string& ckpt_name, SSFlag flag);
    void SaveNodes(const std::string& cat_root, const std::string& ckpt_name, bool ss_incremental, SSFlag flag);
    void PrintSummary() const;
    void LoadExclusions(int id_start, SSFlag flag);

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

    IndexUniv GetIndexCw(size_t iCw) const
    {
        return iCw < rings_.size() ? IndexRing(iCw) : IndexMod(iCw - rings_.size());
    }

    auto& GetCwName(IndexUniv iCw) const
    {
        return iCw.isRing() ? rings_[iCw.index].name : modules_[iCw.index].name;
    }

    std::string GetName(IndexUniv iUniv) const
    {
        switch (iUniv.type) {
        case 0:
            return rings_[iUniv.index].name;
        case 1:
            return modules_[iUniv.index].name;
        case 2:
            return fmt::format("{}:{}", cofseqs_[iUniv.index].name, iUniv.iTri);
        default:
            throw RunTimeError(fmt::format("Incorrect type {}", iUniv.type));
        }
    }

    auto GetMapName(IndexUniv iCs) const -> std::tuple<std::string, std::string, std::string, std::string>;

    auto& GetIndexCofs(IndexUniv iCw) const
    {
        return iCw.isRing() ? rings_[iCw.index].ind_cofs : modules_[iCw.index].ind_cofs;
    }

    auto& GetIndexCofs(IndexUniv iCw)
    {
        return iCw.isRing() ? rings_[iCw.index].ind_cofs : modules_[iCw.index].ind_cofs;
    }

    auto& GetTMax(IndexUniv iCw) const
    {
        return iCw.isRing() ? rings_[iCw.index].t_max : modules_[iCw.index].t_max;
    }

    auto& GetNodesSS(IndexUniv iCw) const
    {
        return iCw.isRing() ? rings_[iCw.index].nodes_ss : modules_[iCw.index].nodes_ss;
    }

    auto& GetNodesSS(IndexUniv iUniv)
    {
        return iUniv.isRing() ? rings_[iUniv.index].nodes_ss : (iUniv.isModule() ? modules_[iUniv.index].nodes_ss : cofseqs_[iUniv.index].nodes_cofseq[iUniv.iTri]);
    }

    auto& GetSSDegs(IndexUniv iCw) const
    {
        return iCw.isRing() ? rings_[iCw.index].degs_ss : modules_[iCw.index].degs_ss;
    }

    IndexUniv GetIndexCwByName(const std::string& name) const
    {
        for (size_t i = 0; i < rings_.size(); ++i)
            if (rings_[i].name == name)
                return IndexRing(i);
        for (size_t i = 0; i < modules_.size(); ++i)
            if (modules_[i].name == name)
                return IndexMod(i);
        return IndexUniv(uint32_t(-1), size_t(-1));
    }
    IndexUniv GetIndexCofByName(const std::string& name) const
    {
        for (size_t i = 0; i < cofseqs_.size(); ++i)
            if (name.starts_with(cofseqs_[i].name))
                return IndexCof(i, uint32_t(name.back() - '0'));
        throw RunTimeError(fmt::format("Incorrect Cs name {}", name));
    }
    AdamsDeg GetMapDeg(IndexUniv iCs) const
    {
        return cofseqs_[iCs.index].degMap[iCs.iTri];
    }
    AdamsDeg GetPrevMapDeg(IndexUniv iCs) const
    {
        return cofseqs_[iCs.index].degMap[size_t(iCs.iTri + 2) % 3];
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
        throw RunTimeError("Incorrect name");
    }

    auto& GetModuleByName(const std::string& name) const
    {
        for (auto& mod : modules_)
            if (mod.name == name)
                return mod;
        throw RunTimeError("Incorrect name");
    }

    /* This is used for plotting Er pages. The actual result might differ by a linear combination.
     * Return a level such that nontrivial diffs in >= level will not be longer.
     */
    static int GetFirstFixedLevelForPlot(const SSNodes& nodes_ss, AdamsDeg deg);
    static int GetFirstFixedLevelForPlotCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg);

public:
    void SetDeduceList(const std::vector<std::string>& cws);
    void SetCofDeduceList(const std::vector<std::string>& deduce_list_cofseq_);
    int1d Multiply(IndexUniv iCw_x, AdamsDeg deg_x, const int1d& x, IndexUniv iCw_y, AdamsDeg deg_y, const int1d& y) const;
    void FillGeneratorNames();
    std::string GetName(IndexUniv iCw_x, AdamsDeg deg_x, const int1d& x) const;
    json Sc2Json(IndexUniv iCw_x, AdamsDeg deg_x) const;
    json Sc2JsonCofseq(IndexUniv iCs, AdamsDeg deg_x) const;
    /* return B_{level} = span{?} */
    std::string ScB(IndexUniv iCw_x, AdamsDeg deg_x, int level_max) const;
    std::string ScBCofseq(IndexUniv iCs, AdamsDeg deg, int level_max) const;
    json Local2Json(IndexUniv iCw_x, AdamsDeg deg_x, int r) const;
    json LocalCofseqToJson(IndexUniv iCs_x, AdamsDeg deg_x, AdamsDeg deg_dx, int r) const;
    std::string Local2Table(IndexUniv iCw_x, int stem, int s_min, int s_max) const;

private:
    /* Count the number of all possible d_r targets. Return (index, count). */
    auto CountPossDrTgt(const SSNodes& nodes_ss, int t_max, const AdamsDeg& deg_tgt, int r) const -> std::pair<int, int>;
    /* Warning: this assumes that there shall not be more Einf elements */
    auto CountPossDrTgtCofseq(const CofSeq& cofseq, size_t iCs, const AdamsDeg& deg_tgt, int r, bool possMorePerm) const -> std::pair<int, int>;

    /* Count the number of all possible d_r sources. Return (index, count). */
    std::pair<int, int> CountPossDrSrc(const SSNodes& nodes_ss, const AdamsDeg& deg_src, int r) const;
    /* Count the number of all potential permant cycles. Return (index, count). */
    auto CountPossMorePerm(const SSNodes& nodes_ss, const AdamsDeg& deg) const -> std::pair<int, int>
    {
        return CountPossDrSrc(nodes_ss, deg, R_PERM - 1);
    };
    /* Warning: this assumes that there shall not be more Einf elements */
    std::pair<int, int> CountPossDrSrcCofseq(const CofSeq& cofseq, size_t iCs, const AdamsDeg& deg_src, int r, bool possMorePerm) const;

    /*
     * Return the smallest r1>=r such that d_{r1} has a possible target
     *
     * Range >= t_max is unknown and thus always possible.
     * Return R_PERM if not found
     */
    int NextRTgt(const SSNodes& nodes_ss, int t_max, AdamsDeg deg, int r) const;
    int NextRTgtCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r) const;

    /**
     * Return the largest r1<=r such that d_{r1} has a possible source
     *
     * Return -1 if not found
     */
    int NextRSrc(const SSNodes& nodes_ss, AdamsDeg deg, int r) const;
    int NextRSrcCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r) const;

private:
    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    SSRet SetDiffSc(IndexUniv iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag);
    /* Add an image. dx must be nonempty and not null. x must be nonempty. */
    SSRet SetImageSc(IndexUniv iCw, AdamsDeg deg_dx, const int1d& dx, const int1d& x, int r, SSFlag flag);

    /* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
    SSRet SetDiffScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x_, const int1d& dx, int r, SSFlag flag);
    /* Add an image. dx must be nonempty and not null. x must be nonempty. */
    SSRet SetImageScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_dx, const int1d& dx_, const int1d& x, int r, SSFlag flag);
    /* Retriangulate when ss changes */
    SSRet ReSetScCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg, SSFlag flag);

public:  //// TOOD: remove this public
    SSRet SyncCofseq(SSFlag flag);
    /* Nullify Adams boundaries in cofseq */
    SSRet SyncToCofseq(SSFlag flag);

    /**
     * Add d_r(x)=dx;
     * Add d_s(xy)=d_s(x)y+xd_s(y) for y with level=LEVEL_MAX-s and s>=r_min.
     * Return the number of changed degrees.
     *
     * dx should not be null.
     */
    SSRet SetRingDiffLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag);
    /** This version deduces d(xy) even if dx is uncertain */
    SSRet SetRingDiffLeibnizV2(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);
    SSRet SetModuleDiffLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, int r_min, SSFlag flag);
    /** This version deduces d(xy) even if dx is uncertain */
    SSRet SetModuleDiffLeibnizV2(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);

public:
    /* Add a node */
    void AddNode(SSFlag flag);

    /* Pop the lastest node */
    void PopNode(SSFlag flag);

    /* Apply the change of the staircase to the current history */
    static void UpdateStaircase(SSNodes& nodes_ss, AdamsDeg deg, const Staircase& sc_i, size_t i_insert, const int1d& x, const int1d& dx, int level, int1d& image, int& level_image);

    /* Cache null diffs to the most recent node. */
    void CacheNullDiffs(const SSNodes& nodes_ss, int t_max, AdamsDeg deg, SSFlag flag, NullDiff1d& nds) const;
    void CacheNullCofDiffs(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, SSFlag flag, NullDiffCofseq1d& nds) const;

    /**
     * Check first if it is a new differential before adding it.
     * Do some deductions by degree.
     *
     * Return the number of changed degrees.
     */
    SSRet SetRingDiffGlobal(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);
    SSRet SetModuleDiffGlobal(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);
    SSRet SetCwDiffGlobal(IndexUniv iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);

    [[nodiscard]] int GetSynImage(IndexUniv iCof, AdamsDeg deg_x, const int1d& x, int level_x, AdamsDeg& deg_fx, int1d& fx, int s_f_dinv_x, int cross_min);
    SSRet SetCwDiffSynthetic(IndexUniv iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool hasCross, SSFlag flag);

    /* Add d_r(?)=x;
     * Add d_r(?)=xy for d_r(y)=0 (y on level < LEVEL_MAX - r);
     */
    SSRet SetRingBoundaryLeibniz(size_t iRing, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);
    SSRet SetModuleBoundaryLeibniz(size_t iMod, AdamsDeg deg_x, const int1d& x, int r, SSFlag flag);

    SSRet SetDiffLeibnizCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag);
    SSRet SetDiffLeibnizCofseq(IndexUniv iRing, AdamsDeg deg_a, const int1d& a, SSFlag flag);
    SSRet SetDiffGlobalCofseq(CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag);
    void SetUnivDiffGlobal(std::string& name, AdamsDeg deg, int r, const int1d& x, const int1d& dx, bool isCs, bool isDInv, SSFlag flag);

public: /* Differentials */
    bool IsNewDiff(const SSNodes& nodes_ss, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const;
    bool IsNewDiffCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r) const;
    SSRet DeduceTrivialCwDiffs(IndexUniv iCw, SSFlag flag);
    SSRet DeduceTrivialCwDiffs(SSFlag flag);
    SSRet DeduceTrivialCofDiffs(size_t iCs_index, SSFlag flag);
    SSRet DeduceTrivialCofDiffs(SSFlag flag);
    SSRet DeduceTrivialDiffs(SSFlag flag);
    SSRet DeduceManual(SSFlag flag);
    SSRet DeduceDiffBySynthetic(SSFlag flag);
    SSRet DeduceDiffBySyntheticCofseq(SSFlag flag);

    SSRet TryDiff(IndexUniv iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, SSFlag flag, bool tryY);
    SSRet DeduceDiff4Nd(IndexUniv iCw, AdamsDeg deg, const NullDiff& nd, SSFlag flag);
    SSRet DeduceDiff4X(IndexUniv iCw, AdamsDeg deg, int1d x, int level, SSFlag flag);
    SSRet DeduceDiff4XDepth(IndexUniv iCw, AdamsDeg deg, int1d x, int original_level, SSFlag flag); /* For SSFlag::deduce_zero */
    SSRet DeduceDiffs4Deg(IndexUniv iCw, AdamsDeg deg, SSFlag flag);
    SSRet DeduceDiffs(IndexUniv& iCw, int stem_min, int stem_max, SSFlag flag);
    SSRet DeduceDiffs(int stem_min, int stem_max, int T, int id_thread, SSFlag flag);
    /* Deduce d(xy) no matter what dx is */
    SSRet DeduceDiffsV2(SSFlag flag);

    SSRet TryCofDiff(IndexUniv iCof, AdamsDeg deg_x, AdamsDeg deg_dx, const int1d& x, const int1d& dx, const int1d& perm, int r, SSFlag flag, bool tryY);
    SSRet DeduceCofDiff4Nd(IndexUniv iCof, AdamsDeg deg, const NullDiffCofseq& nd, SSFlag flag);
    SSRet DeduceCofDiff4X(IndexUniv iCof, AdamsDeg deg, int1d x, int level, SSFlag flag);
    SSRet DeduceCofDiff4XDepth(IndexUniv iCof, AdamsDeg deg, int1d x, int original_level, SSFlag flag);
    SSRet DeduceCofDiffs4Deg(IndexUniv iCof, AdamsDeg deg, SSFlag flag);
    SSRet DeduceCofDiffs(int stem_min, int stem_max, int T, int id_thread, SSFlag flag);
    SSRet DeduceCofDiffsNbhd(IndexUniv iCof, int stem, SSFlag flag);

    SSRet DeduceDiff4Logic(IndexUniv iUniv, AdamsDeg deg, int1d x, int level, SSFlag flag);

    SSRet CommuteCofseq(size_t iComm, SSFlag flag);
    SSRet CommuteCofseq(SSFlag flag);

public:
    /* Return if Einf at deg is possibly nontrivial */
    static void UpdatePossEinf(const SSNodes& nodes_ss, ut::map_seq2d<int, 0>& basis_ss_possEinf);
    void UpdateAllPossEinf()
    {
        for (auto& ring : rings_)
            UpdatePossEinf(ring.nodes_ss, ring.basis_ss_possEinf);
        for (auto& mod : modules_)
            UpdatePossEinf(mod.nodes_ss, mod.basis_ss_possEinf);
    }
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
    // void SetPermanentCycle(int depth, size_t iCof, AdamsDeg deg_x);
    void AddPiRelsRing(size_t iRing, algZ::Poly1d rels);
    void AddPiRelsCof(size_t iMod, algZ::Mod1d rels);
    // void AddPiRelsByNat(size_t iMod);
    void SimplifyPiRels()
    {
        for (auto& ring : rings_)
            ring.pi_gb.SimplifyRels();
        for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
            modules_[iCof].pi_gb.SimplifyRels();
    }
    static algZ::Mon1d GenBasis(const algZ::Groebner& gb, AdamsDeg deg, const PiBasis1d& nodes_pi_basis);
    static algZ::MMod1d GenBasis(const algZ::GroebnerMod& gb, AdamsDeg deg, const PiBasis1d& basis);
    // void SyncS0Homotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopyy, int depth);
    // void SyncCofHomotopy(int iCof, AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth);
    /*void SyncHomotopy(AdamsDeg deg_min, int& count_ss, int& count_homotopy, int depth)
    {
        SyncS0Homotopy(deg_min, count_ss, count_homotopy, depth);
        for (size_t iCof = 0; iCof < modules_.size(); ++iCof)
            SyncCofHomotopy((int)iCof, deg_min, count_ss, count_homotopy, depth);
    }*/
    // int DeduceTrivialExtensions(int depth);
    // int DeduceExtensions2tor();
    // int DeduceExtensionsByExactness(int stem_min, int stem_max, int depth);

    // unsigned TryExtS0(algZ::Poly rel, AdamsDeg deg_change, int depth, SSFlag flag);
    // unsigned TryExtCof(size_t iCof, algZ::Mod rel, AdamsDeg deg_change, int depth, SSFlag flag);
    // unsigned TryExtQ(size_t iCof, size_t gen_id, algZ::Poly q, AdamsDeg deg_change, int depth, SSFlag flag);
    // void DeduceExtensions(int stem_min, int stem_max, int& count_ss, int& count_homotopy, int depth, SSFlag flag);
    // int DefineDependenceInExtensions(int stem_min, int stem_max, int depth);  ////;
    // int DefineDependenceInExtensionsV2(int stem_min, int stem_max, int stem_max_mult, int depth);
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

    void drop_and_create_cofseq(const std::string& table) const
    {
        drop_table(table);
        create_cofseq(table);
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
    void update_ss(const std::string& table_prefix, const SSNodes& nodes_ss) const;
    void save_ss(const std::string& table_prefix, const SSNodes& nodes_ss) const;
    void save_cofseq(const std::string& table, const CofSeq& cofseq, bool incremental) const;
    void save_pi_generators_mod(const std::string& table_prefix, const AdamsDeg1d& gen_degs, const Mod1d& gen_Einf) const;
    void save_pi_basis(const std::string& table_prefix, const PiBasis& basis) const;
    void save_pi_basis_mod(const std::string& table_prefix, const PiBasisMod& basis) const;
    void save_pi_def(const std::string& table_prefix, const std::vector<EnumDef>& pi_gen_defs, const std::vector<std::vector<GenConstraint>>& pi_gen_def_mons) const;

public:
    /* load the minimum id in every degree */
    std::map<AdamsDeg, int> load_ss_indices(const std::string& table_prefix) const;  // TODO: change
    SS load_ss(const std::string& table_prefix) const;
    std::array<SS, 3> load_cofseq(const std::string& table) const;
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

/*--------------------------------------------------------------------------------------------
----------------------------------------   Staircase  ----------------------------------------
---------------------------------------------------------------------------------------------*/

size_t GetFirstIndexOnLevel(const Staircase& sc, int level);
size_t GetFirstIndexOfNullOnLevel(const Staircase& sc, int level);
int GetMaxLevelWithND(const Staircase& sc);
bool IsZeroOnLevel(const Staircase& sc, const int1d& x, int level);

/*--------------------------------------------------------------------------------------------
------------------------------------------    ss    ------------------------------------------
---------------------------------------------------------------------------------------------*/

/* Warning: The following IsPoss functions do not check if ss[deg] exists */
/* Check if deg can be hit by dr for r<=r_max */
bool IsPossTgt(const SSNodes& nodes_ss, AdamsDeg deg, int r_max);
inline bool IsPossTgt(const SSNodes& nodes_ss, AdamsDeg deg)
{
    return IsPossTgt(nodes_ss, deg, R_PERM);
}

/* Return if Einf at deg could have nontrivial elements */
int PossEinf(const SSNodes& nodes_ss, AdamsDeg deg);
/* Return if Einf at deg could have more nontrivial elements */
int PossMoreEinf(const SSNodes& nodes_ss, AdamsDeg deg);

/* Return the level of x or d^{-1}x */
int1d GetDiff(const SSNodes& nodes_ss, AdamsDeg deg_x, const int1d& x, int r);
int GetLevel(const SSNodes& nodes_ss, AdamsDeg deg_x, int1d x);
std::pair<int1d, int> GetDiffAndLevel(const SSNodes& nodes_ss, AdamsDeg deg_x, int1d x);
/* When level of x is smaller than 5000 we set r=R_PERM-1 and diff={} */
void GetRAndDiff(const SSNodes& nodes_ss, AdamsDeg deg_x, int1d x, int& r, int1d& diff);

/* Return the minimal length of the crossing differentials */
int GetCrossR(const SSNodes& nodes_ss, AdamsDeg deg, int t_max, int Er);

/* Return the first index with level >=`level_min` such that all levels above are already fixed
 * We do not support level_min < LEVEL_PERM yet
 */
size_t GetFirstIndexOfFixedLevels(const SSNodes& nodes_ss, AdamsDeg deg, int level_min);

/*--------------------------------------------------------------------------------------------
------------------------------------------    cofseq    --------------------------------------
---------------------------------------------------------------------------------------------*/
inline size_t NextiTri(auto iTri)
{
    return size_t(iTri + 1) % 3;
}
inline size_t PreviTri(auto iTri)
{
    return size_t(iTri + 2) % 3;
}
bool IsPossTgtCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int r_max);
size_t GetFirstIndexOfFixedLevelsCofseq(const CofSeq& cofseq, size_t iCs, AdamsDeg deg, int level_min);
int GetCofseqCrossR(const SSNodes& nodes_cofseq, const SSNodes& nodes_ss, AdamsDeg deg, int t_max, int r_min, int result_min);

/*--------------------------------------------------------------------------------------------
------------------------------------------    Utilities    -----------------------------------
---------------------------------------------------------------------------------------------*/

/* Compute x mod (level-1) */
inline int1d& ResidueInplace(int1d& x, const Staircase& sc, int level)
{
    size_t first_l = GetFirstIndexOnLevel(sc, level);
    lina::ResidueInplace(sc.basis.begin(), sc.basis.begin() + first_l, x);
    return x;
}

/* Compute x mod (level-1) */
inline int1d Residue(int1d x, const SSNodes& nodes_ss, AdamsDeg deg, int level)
{
    if (x.empty())
        return x;
    ResidueInplace(x, nodes_ss.GetRecentSc(deg), level);
    return x;
}

void GetAllDbNames(const std::string& cat_name, std::vector<std::string>& names, std::vector<std::string>& paths, std::vector<int>& isRing, bool log = false);
std::pair<std::string, std::string> ParseCatName(const std::string& cat_name);

std::string Str(const Mon& m, const std::vector<std::string>& gen_names);
std::string Str(const Poly& p, const std::vector<std::string>& gen_names);
std::string Str(const MMod& m, const std::vector<std::string>& gen_names, const std::vector<std::string>& v_names);
std::string Str(const Mod& p, const std::vector<std::string>& gen_names, const std::vector<std::string>& v_names);

#endif
