/** \file groebner_steenrod_cpt.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_STEENROD_CPT_H
#define GROEBNER_STEENROD_CPT_H

#include "algebras/steenrod.h"

//#define MYDEPLOY

namespace steenrod {

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

class SteenrodMResConst
{
    using TypeIndices = std::vector<std::unordered_map<uint64_t, int1d>>;

private:
    int t_trunc_, s_trunc_; /* t < t_trunc_ or (t == t_trunc_ and s >= s_trunc_) */

    DataMResConst2d gb_;
    MMod2d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int2d basis_degrees_; /* `basis_degrees[s][i]` is the degree of v_{s,i} */

public:
    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    SteenrodMResConst(int t_trunc, int s_trunc, DataMResConst2d data, int2d basis_degrees);
    static SteenrodMResConst load(const std::string& filename);

public:
    int t_trunc() const
    {
        return t_trunc_;
    }

    int t_begin() const
    {
        int t = 0;
        for (size_t s = 0; s < gb_.size(); ++s) {
            if (!gb_[s].empty()) {
                auto lead = gb_[s].back().x1.GetLead();
                int t1 = lead.deg_m() + basis_degrees_[s][lead.v()];
                if (t1 > t)
                    t = t1;
            }
        }
        return t;
    }

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

void compute_products_ind(int t_trunc, std::string& db_in, std::string& db_out);
void compute_products(int t_trunc, std::string& db_in, std::string& db_out);

}  // namespace steenrod

#endif