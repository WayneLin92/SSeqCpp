/** \file groebner_steenrod_cpt.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_STEENROD_CONST_H
#define GROEBNER_STEENROD_CONST_H

#include "algebras/steenrod.h"
#include <map>

using namespace steenrod;

struct AdamsDegV2
{
    int s, t;
    constexpr AdamsDegV2() : s(0), t(0) {}
    constexpr AdamsDegV2(int s_, int t_) : s(s_), t(t_) {}

    /* Order by (t, -s) */
    bool operator<(const AdamsDegV2& rhs) const
    {
        if (t < rhs.t)
            return true;
        if (t > rhs.t)
            return false;
        if (s > rhs.s)
            return true;
        return false;
    };
    bool operator==(const AdamsDegV2& rhs) const
    {
        return s == rhs.s && t == rhs.t;
    };
    int stem() const
    {
        return t - s;
    }
};

int1d HomToK(const Mod& x);

struct DataMResConst
{
    Mod x1, x2;
    DataMResConst() {}
    DataMResConst(Mod x1_, Mod x2_) : x1(std::move(x1_)), x2(std::move(x2_)) {}
    DataMResConst& operator+=(const DataMResConst& rhs)
    {
        x1 += rhs.x1;
        x2 += rhs.x2;
        return *this;
    }
};

using DataMResConst1d = std::vector<DataMResConst>;
using DataMResConst2d = std::vector<DataMResConst1d>;

inline void Reduce(Mod& x, const DataMResConst1d& y, Mod& tmp)
{
    for (size_t i = 0; i < y.size(); ++i)
        if (std::binary_search(x.data.begin(), x.data.end(), y[i].x1.GetLead()))
            x.iaddP(y[i].x1, tmp);
}

class DbAdamsResLoader : public myio::Database
{
    using Statement = myio::Statement;

public:
    explicit DbAdamsResLoader(const std::string& filename) : Database(filename) {}

public:
    /**
     * loc2glo[s][i] = id
     * glo2loc[id] = (s, i)
     *
     * From SteenrodMRes
     */
    void load_id_converter(const std::string& table_in, int2d& loc2glo, std::vector<std::pair<int, int>>& glo2loc) const;
    int2d load_basis_degrees(const std::string& table_prefix, int t_trunc) const;
    void load_generators(const std::string& table_prefix, std::map<int, AdamsDegV2>& id_st, int2d& vid_num, std::map<AdamsDegV2, Mod1d>& diffs, int t_trunc) const;
    DataMResConst2d load_data(const std::string& table_prefix, int t_trunc) const;
};

class AdamsResConst
{
    using TypeIndices = std::vector<std::unordered_map<uint64_t, int1d>>;

private:

    DataMResConst2d gb_;
    MMod2d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int2d basis_degrees_; /* `basis_degrees[s][i]` is the degree of v_{s,i} */

public:
    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    AdamsResConst(DataMResConst2d data, int2d basis_degrees);
    static AdamsResConst load(const DbAdamsResLoader& db, int t_trunc);

public:
    const int1d& basis_degrees(size_t s) const
    {
        return basis_degrees_[s];
    }

    int dim_Ext() const
    {
        size_t dim = 0;
        for (size_t s = 0; s < basis_degrees_.size(); ++s)
            dim += basis_degrees_[s].size();
        return (int)dim;
    }

    int dim_Gb() const
    {
        size_t dim = 0;
        for (size_t s = 0; s < gb_.size(); ++s)
            dim += gb_[s].size();
        return (int)dim;
    }

    void push_back(DataMResConst g, size_t s)
    {
        MMod m = g.x1.GetLead();

        leads_[s].push_back(m);
        indices_[s][m.v_raw()].push_back((int)gb_[s].size());
        gb_[s].push_back(std::move(g));
    }

    const auto& data() const
    {
        return gb_;
    }

    Mod DiffInv(Mod x, size_t s) const;
    void DiffInvBatch(Mod1d xs, Mod1d& xs_reduced, size_t s) const;
};

struct GenMRes
{
    int id, t;
    Mod diff;
};

using GenMRes1d = std::vector<GenMRes>;
using GenMRes2d = std::vector<GenMRes1d>;

#endif