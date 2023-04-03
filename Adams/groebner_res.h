/** \file groebner_steenrod.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_RES_H
#define GROEBNER_RES_H

#include "algebras/benchmark.h"
#include "algebras/groebner_steenrod.h"
#include <map>
#include <unordered_map>
#include <unordered_set>

using namespace steenrod;

struct Filtr
{
    uint64_t data;
    Filtr() : data(~0) {}
    Filtr(MMod m) : data(m.v_raw() | m.w_may()) {}
    Filtr(uint64_t data_) : data(data_) {}
    bool operator<(Filtr rhs) const
    {
        return data < rhs.data;
    }
    bool operator==(Filtr rhs) const
    {
        return data == rhs.data;
    }
    Filtr operator+(uint64_t w_may) const
    {
        return Filtr(data + w_may);
    }
};

/********************************************************
 *                    class GroebnerX2m
 ********************************************************/

class GroebnerX2m
{
    using TypeIndices = std::vector<std::unordered_map<uint64_t, int1d>>;

private:
    int t_trunc_, stem_trunc_;

    CriMilnors1d criticals_; /* Groebner basis of critical pairs */

    Mod2d gb_;
    MMod2d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int2d basis_degrees_; /* `basis_degrees_x2m_[s][i]` is the degree of w_{s,i} */

public:
    GroebnerX2m(int t_trunc, int stem_trunc, Mod2d data, int2d basis_degrees, std::map<int, int>& latest_st);

public:
    void resize(size_t s)
    {
        if (basis_degrees_.size() < s + 1)
            basis_degrees_.resize(s + 1);
        if (gb_.size() < s) {
            gb_.resize(s);
            leads_.resize(s);
            indices_.resize(s);
            while (criticals_.size() < s)
                criticals_.push_back(CriMilnors(t_trunc_));
        }
    }

    Mod new_gen(size_t s, int t)
    {
        basis_degrees_[s].push_back(t);
        return MMod(MMilnor(), basis_degrees_[s].size() - 1);
    }

    void push_back(Mod g, size_t s)
    {
        MMod m = g.GetLead();
        criticals_[s].AddToBuffers(leads_[s], m, basis_degrees_[s][m.v()]);

        leads_[s].push_back(m);
        indices_[s][m.v_raw()].push_back((int)gb_[s].size());
        gb_[s].push_back(std::move(g));
    }

    const int1d& basis_degrees(size_t s)
    {
        return basis_degrees_[s];
    }

    const int2d& basis_degrees()
    {
        return basis_degrees_;
    }

    const auto& data() const
    {
        return gb_;
    }

    Mod Reduce(Mod x2m, size_t s) const;
    Mod Reduce(const CriMilnor& p, size_t s) const;
    Mod1d AddRels(size_t s, int t);
};

/********************************************************
 *                    class AdamsRes
 ********************************************************/

struct DataMRes
{
    Mod x1, x2, x2m;
    Filtr fil;
    DataMRes() {}
    DataMRes(Mod x1_, Mod x2_, Mod x2m_) : x1(std::move(x1_)), x2(std::move(x2_)), x2m(std::move(x2m_))
    {
        if (x1)
            fil = Filtr(x1.GetLead());
    }
    void iaddP(const DataMRes& rhs, Mod& tmp_Mod)
    {
#ifndef NDEBUG
        if (!rhs.valid_x2m())
            throw MyException(0x2ae8baa3U, "Add only when rhs.x2m is valid");
#endif
        if (valid_x2m() && fil == rhs.fil)
            x2m.iaddP(rhs.x2m, tmp_Mod);
        x1.iaddP(rhs.x1, tmp_Mod);
        x2.iaddP(rhs.x2, tmp_Mod);
    }
    bool valid_x2m() const
    {
        if (x1)
            return fil == Filtr(x1.GetLead());
        if (x2m)
            return false;
        return true;
    }
};

using DataMRes1d = std::vector<DataMRes>;
using DataMRes2d = std::vector<DataMRes1d>;

inline void Reduce(Mod& x, const DataMRes1d& y, Mod& tmp)
{
    for (size_t i = 0; i < y.size(); ++i)
        if (std::binary_search(x.data.begin(), x.data.end(), y[i].x1.GetLead()))
            x.iaddP(y[i].x1, tmp);
}

class AdamsRes
{
    using TypeIndices = std::vector<std::unordered_map<uint64_t, int1d>>;

private:
    int t_trunc_;
    int stem_trunc_;

    CriMilnors1d criticals_; /* Groebner basis of critical pairs */

    DataMRes2d gb_;
    MMod2d leads_;        /* Leading monomials */
    TypeIndices indices_; /* Cache for fast divisibility test */

    int2d basis_degrees_; /* `basis_degrees[s][i]` is the degree of v_{s,i} */

    GroebnerX2m gb_x2m_;

public:
    /* Initialize from `polys` which already forms a Groebner basis. Must not add more relations. */
    AdamsRes(int t_trunc, int stem_trunc, DataMRes2d data, int2d basis_degrees, Mod2d data_x2m, int2d basis_degrees_x2m, std::map<int, int>& latest_st);
    static AdamsRes load(const std::string& db_filename, const std::string& tablename, int t_trunc, int stem_trunc);

public:
    int t_trunc() const
    {
        return t_trunc_;
    }

    const int1d& basis_degrees(size_t s) const
    {
        return basis_degrees_[s];
    }

    void set_v_degrees(int1d v_degs)
    {
        if (basis_degrees_.empty())
            basis_degrees_.push_back(std::move(v_degs));
        else
            basis_degrees_[0] = std::move(v_degs);
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

    const int1d& basis_degrees_x2m(size_t s)
    {
        return gb_x2m_.basis_degrees(s);
    }

    const int2d& basis_degrees_x2m()
    {
        return gb_x2m_.basis_degrees();
    }

    void resize_gb(size_t s)
    {
        if (basis_degrees_.size() < s + 1)
            basis_degrees_.resize(s + 1);
        if (gb_.size() < s) {
            gb_.resize(s);
            leads_.resize(s);
            indices_.resize(s);
            while (criticals_.size() < s)
                criticals_.push_back(CriMilnors(t_trunc_));
            if (s > 0)
                gb_x2m_.resize(s - 1);
        }
    }

    Mod new_gen(size_t s, int t)
    {
        basis_degrees_[s].push_back(t);
        return MMod(MMilnor(), basis_degrees_[s].size() - 1);
    }

    Mod new_gen_x2m(size_t s, int t)
    {
        return gb_x2m_.new_gen(s, t);
    }

    void push_back(DataMRes g, size_t s)
    {
        MMod m = g.x1.GetLead();
        criticals_[s].AddToBuffers(leads_[s], m, basis_degrees_[s][m.v()]);

        leads_[s].push_back(m);
        indices_[s][m.v_raw()].push_back((int)gb_[s].size());
        gb_[s].push_back(std::move(g));
    }

    void push_back_x2m(Mod g, size_t s)
    {
        gb_x2m_.push_back(std::move(g), s);
    }

    const auto& data() const
    {
        return gb_;
    }

    const auto& data_x2m() const
    {
        return gb_x2m_.data();
    }

    CriMilnor1d Criticals(size_t s, int t, Mod1d& rels_x2m);
    DataMRes Reduce(const CriMilnor& cp, size_t s) const;
    void ReduceBatch(const CriMilnor1d& cp, DataMRes1d& results, size_t s) const;
    Mod Reduce(Mod x, size_t s) const;
    Mod ReduceX2m(const CriMilnor& cp, size_t s) const;
    Mod ReduceX2m(Mod x2m, size_t s) const
    {
        return gb_x2m_.Reduce(std::move(x2m), s);
    }
};

/**
 * Comsume relations from 'rels` and `gb.criticals_` in degree `<= deg`.
 *
 * return the dimension of the calculated range for debugging.
 */
void Resolve(AdamsRes& gb, const Mod1d& rels, const int1d& v_degs, int t_max, int stem_max, const std::string& db_filename, const std::string& tablename);

void ResetDb(const std::string& filename, const std::string& tablename);

#endif