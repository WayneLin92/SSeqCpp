#include "groebner_steenrod.h"
#include "algebras/benchmark.h"  ////
#include "algebras/dbalg.h"
#include <cstring>
//#include <immintrin.h>

namespace steenrod {

/********************************************************
 *                    class CriMilnors
 ********************************************************/

CriMilnor1d CriMilnors::Criticals(int t)
{
    CriMilnor1d result;
    if (!buffer_singles_.empty() && buffer_singles_.begin()->first == t) {
        std::swap(result, buffer_singles_.begin()->second);
        buffer_singles_.erase(t);
    }
    if (!buffer_min_pairs_.empty() && buffer_min_pairs_.begin()->first == t) {
        auto& b_min_pairs_t = buffer_min_pairs_.begin()->second;
        for (size_t j = 0; j < b_min_pairs_t.size(); ++j)
            for (auto& pair : b_min_pairs_t[j])
                if (pair.i2 != -1)
                    result.push_back(pair);
        buffer_min_pairs_.erase(buffer_min_pairs_.begin());
    }

    return result;
}

void CriMilnors::Minimize(const MMod1d& leads, int t)
{
    if (buffer_redundent_pairs_.empty() || buffer_redundent_pairs_.begin()->first != t)
        return;
    if (buffer_min_pairs_.empty() || buffer_min_pairs_.begin()->first != t)
#ifdef NDEBUG
        return;
#else
        buffer_min_pairs_[t];
#endif

    constexpr uint64_t NULL_J = ~uint64_t(0);
    auto& b_min_pairs_t = buffer_min_pairs_.begin()->second;
    for (uint64_t ij : buffer_redundent_pairs_.begin()->second) {
        uint64_t i, j;
        ut::UnBind(ij, i, j);
        while (j != NULL_J) {
            MMilnor gcd = gcdLF(leads[i].m_no_weight(), leads[j].m_no_weight());
            MMilnor m2 = divLF(leads[i].m(), gcd);
            if (j < (int)b_min_pairs_t.size()) {
                auto p = std::find_if(b_min_pairs_t[j].begin(), b_min_pairs_t[j].end(), [&m2](const CriMilnor& c) { return c.m2 == m2; });
                /* Mark it to be removed from `buffer_min_pairs_` */
                if (p != b_min_pairs_t[j].end()) {
                    p->i2 = -1;
                    break;
                }
            }

            /* Reduce (i, j) */
            auto c = gb_[j].begin();
            auto end = gb_[j].end();
            for (; c < end; ++c) {
                if (divisibleLF(c->m2, m2)) {
                    MMilnor m1 = divLF(leads[j].m(), gcd);
                    if (gcdLF(c->m1, m1)) {
                        j = NULL_J;
                        break;
                    }
                    else {
                        j = c->i1;
                        if (i > j)
                            std::swap(i, j);
                        break;
                    }
                }
            }
#ifndef NDEBUG
            if (c == end) {
                std::cout << "i=" << i << '\n';
                std::cout << "j=" << j << '\n';
                std::cout << "leads[i]=" << leads[i].Str() << '\n';
                std::cout << "leads[j]=" << leads[j].Str() << '\n';
                std::cout << "m2=" << Milnor(m2) << '\n';
                for (size_t k = 0; k < gb_[j].size(); ++k)
                    std::cout << "gb_[j][k].m2" << Milnor(gb_[j][k].m2) << '\n';
                throw MyException(0xfa5db14U, "Should not happen because gb_ is groebner");
            }
#endif
        }
    }

    /* Delete `buffer_redundent_pairs_[t]` */
    buffer_redundent_pairs_.erase(buffer_redundent_pairs_.begin());
}

/**
 * Populate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
 * `buffer_min_pairs_` will be reduced later by `CPMilnors::Minimize()`.
 */
void CriMilnors::AddToBuffers(const MMod1d& leads, MMod mon, int t_v)
{
    size_t lead_size = leads.size();
    std::vector<std::pair<int, CriMilnor>> new_pairs(lead_size);

    /* Populate `new_pairs` */
    for (size_t i = 0; i < leads.size(); ++i) {
        new_pairs[i].first = -1;
        if (leads[i].v_raw() == mon.v_raw()) {
            int d_pair = lcmLF(leads[i].m_no_weight(), mon.m_no_weight()).deg() + t_v;
            if (d_pair <= t_trunc_) {
                new_pairs[i].first = d_pair;
                CriMilnor::SetFromLM(new_pairs[i].second, leads[i].m(), mon.m(), (int)i, (int)lead_size);
            }
        }
    }

    /* Remove some new pairs to form Groebner basis and discover redundent pairs */
    for (size_t j = 1; j < new_pairs.size(); ++j) {
        if (new_pairs[j].first != -1) {
            for (size_t i = 0; i < j; ++i) {
                if (new_pairs[i].first != -1) {
                    if (divisibleLF(new_pairs[i].second.m2, new_pairs[j].second.m2)) {
                        new_pairs[j].first = -1;
                        break;
                    }
                    else if (divisibleLF(new_pairs[j].second.m2, new_pairs[i].second.m2)) {
                        new_pairs[i].first = -1;
                    }
                    else if (!gcdLF(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                        int dij = lcmLF(leads[i].m_no_weight(), leads[j].m_no_weight()).deg() + t_v;
                        if (dij <= t_trunc_)
                            buffer_redundent_pairs_[dij].insert(ut::Bind(i, j));
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < new_pairs.size(); ++i) {
        if (new_pairs[i].first != -1) {
            gb_.resize(lead_size + 1);
            gb_[lead_size].push_back(new_pairs[i].second);
            buffer_min_pairs_[new_pairs[i].first].resize(lead_size + 1);
            buffer_min_pairs_[new_pairs[i].first][lead_size].push_back(new_pairs[i].second);
        }
    }

    /* Populate `buffer_singles_` */
    int t_mon = mon.deg_m() + t_v;
    for (int i : mon.m_no_weight()) {
        MMilnor m = MMilnor::FromIndex(i);
        int t = t_mon + m.deg();
        if (t <= t_trunc_)
            buffer_singles_[t].push_back(CriMilnor::Single(m, (int)lead_size));
    }
}

void CriMilnors::init(const MMod1d& leads, const array& basis_degrees)
{
    MMod1d tmp_leads;
    for (auto& mon : leads) {
        int t_v = basis_degrees[mon.v()];
        AddToBuffers(tmp_leads, mon, t_v);
        tmp_leads.push_back(mon);
    }

    int t_min_buffer = leads.back().deg_m() + basis_degrees[leads.back().v()] + 1;
    for (auto it = buffer_min_pairs_.begin(); it != buffer_min_pairs_.end(); ++it) {
        if (it->first < t_min_buffer)
            it = buffer_min_pairs_.erase(it);
        else
            ++it;
    }
    for (auto it = buffer_redundent_pairs_.begin(); it != buffer_redundent_pairs_.end(); ++it) {
        if (it->first < t_min_buffer)
            it = buffer_redundent_pairs_.erase(it);
        else
            ++it;
    }
    for (auto it = buffer_singles_.begin(); it != buffer_singles_.end(); ++it) {
        if (it->first < t_min_buffer)
            it = buffer_singles_.erase(it);
        else
            ++it;
    }
}

/********************************************************
 *                    class GroebnerX2m
 ********************************************************/

/** Return i such that mon divides leads[i].
 *
 * Return -1 if not found.
 * @param indices indices[v_raw] stores the indices of elements of leads with the given v_raw.
 */
int IndexOfDivisibleLeading(const MMod1d& leads, const std::unordered_map<uint64_t, array>& indices, MMod mon)
{
    auto key = mon.v_raw();
    auto p = indices.find(key);
    if (p != indices.end())
        for (int k : p->second)
            if (divisibleLF(leads[k], mon))
                return k;
    return -1;
}

Mod GroebnerX2m::Reduce(Mod x2m, size_t s) const
{
    Mod tmp;
    tmp.data.reserve(50);

    size_t index;
    index = 0;
    while (index < x2m.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], x2m.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(x2m.data[index], gb_[s][gb_index].data[0]);
            x2m.iaddmulMay(m, gb_[s][gb_index], tmp);
        }
        else
            ++index;
    }

    return x2m;
}

Mod GroebnerX2m::Reduce(const CriMilnor& p, size_t s) const
{
    Mod result;

    Mod tmp;
    tmp.data.reserve(50);

    if (p.i1 >= 0)
        result.iaddmulMay(p.m1, gb_[s][p.i1], tmp).iaddmulMay(p.m2, gb_[s][p.i2], tmp);
    else
        result.iaddmulMay(p.m2, gb_[s][p.i2], tmp);

    size_t index;
    index = 0;
    while (index < result.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], result.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.data[index], gb_[s][gb_index].data[0]);
            result.iaddmulMay(m, gb_[s][gb_index], tmp);
        }
        else
            ++index;
    }

    return result;
}

void GroebnerX2m::AddRels(size_t s, int t)
{
    /* Populate `rels_tmp` */
    resize(s + 1);
    Mod1d rels_tmp;

    criticals_[s].Minimize(leads_[s], t);
    CriMilnor1d pairs_st = criticals_[s].Criticals(t);
    if (!pairs_st.empty()) {
        rels_tmp.resize(pairs_st.size());
        ut::for_each_seq((int)rels_tmp.size(), [&](size_t i) { rels_tmp[i] = Reduce(pairs_st[i], s); });
    }

    if (rels_tmp.empty())
        return;

    /* Triangulate these relations */
    Mod1d rels_st;
    for (size_t i = 0; i < rels_tmp.size(); ++i) {
        for (size_t j = 0; j < rels_st.size(); ++j)
            if (std::binary_search(rels_tmp[i].data.begin(), rels_tmp[i].data.end(), rels_st[j].GetLead()))
                rels_tmp[i] += rels_st[j];
        if (rels_tmp[i])
            rels_st.push_back(std::move(rels_tmp[i]));
    }

    /* Add these relations */
    for (size_t i = 0; i < rels_st.size(); ++i)
        push_back(rels_st[i], s);

    /*db.begin_transaction();
    db.save_generators("SteenrodMRes", kernel_sp1, (int)s + 2, t);
    if (s >= 0)
        db.save_relations("SteenrodMRes", data_st, gb.basis_degs()[s], s);
    db.save_relations("SteenrodMRes", rels_splus, gb.basis_degs()[size_t(s + 1)], s + 1);
    db.end_transaction();*/
}

/********************************************************
 *                    class GroebnerMRes
 ********************************************************/

GroebnerMRes::GroebnerMRes(int t_trunc, DataMRes2d data, array2d basis_degrees, Mod2d data_x2m, array2d basis_degrees_x2m)
    : t_trunc_(t_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees)), gb_x2m_(t_trunc, std::move(data_x2m), std::move(basis_degrees_x2m))
{
    if (basis_degrees_.empty())
        basis_degrees_.push_back({0});
    if (basis_degrees_[0].empty())
        basis_degrees_[0].push_back(0);

    leads_.resize(gb_.size());
    indices_.resize(gb_.size());

    for (size_t s = 0; s < gb_.size(); ++s) {
        for (int j = 0; j < (int)gb_[s].size(); ++j) {
            leads_[s].push_back(gb_[s][j].x1.GetLead());
            indices_[s][gb_[s][j].x1.GetLead().v_raw()].push_back(j);
        }
    }

    for (size_t s = 0; s < gb_.size(); ++s) {
        criticals_.push_back(CriMilnors(t_trunc_));
        criticals_.back().init(leads_[s], basis_degrees_[s]);
    }

    ////x2m
}

CriMilnor1d GroebnerMRes::Criticals(size_t s, int t)
{
    gb_x2m_.AddRels(s, t);
    criticals_[s].Minimize(leads_[s], t);
    CriMilnor1d cris = criticals_[s].Criticals(t);
    std::vector<Filtr> fils(cris.size());
    for (size_t i = 0; i < cris.size(); ++i)
        fils[i] = gb_[s][cris[i].i2].fil + cris[i].m2.w_may();
    auto indices = ut::size_t_range(cris.size());
    std::sort(indices.begin(), indices.end(), [&fils](size_t i, size_t j) { return fils[j] < fils[i]; });
    CriMilnor1d result;
    result.reserve(cris.size());
    for (size_t i = 0; i < cris.size(); ++i)
        result.push_back(cris[indices[i]]);

    Mod1d x2ms;
    for (auto& cp : result) {
        Mod x2m = ReduceX2m(cp, s);
        for (size_t j = 0; j < x2ms.size(); ++j)
            if (std::binary_search(x2m.data.begin(), x2m.data.end(), x2ms[j].GetLead()))
                x2m += x2ms[j];
        if (x2m)
            x2ms.push_back(std::move(x2m));
        else {
            bench::Counter(4);
            cp.i2 = -1;
        }
    }
    ut::RemoveIf(result, [](const CriMilnor& cp) { return cp.i2 == -1; });
    return result;
}

DataMRes GroebnerMRes::Reduce(const CriMilnor& cp, size_t s) const
{
    DataMRes result;

    Milnor tmp_a;
    Mod tmp_m1;
    Mod tmp_m2;
    tmp_a.data.reserve(50);
    tmp_m1.data.reserve(100);
    tmp_m2.data.reserve(100);

    if (cp.i1 >= 0) {
        result.x1.iaddmul(cp.m1, gb_[s][cp.i1].x1, tmp_a, tmp_m1, tmp_m2).iaddmul(cp.m2, gb_[s][cp.i2].x1, tmp_a, tmp_m1, tmp_m2);
        result.x2.iaddmul(cp.m1, gb_[s][cp.i1].x2, tmp_a, tmp_m1, tmp_m2).iaddmul(cp.m2, gb_[s][cp.i2].x2, tmp_a, tmp_m1, tmp_m2);
        result.x2m.iaddmulMay(cp.m1, gb_[s][cp.i1].x2m, tmp_m1).iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_m1);
    }
    else {
        result.x1.iaddmul(cp.m2, gb_[s][cp.i2].x1, tmp_a, tmp_m1, tmp_m2);
        result.x2.iaddmul(cp.m2, gb_[s][cp.i2].x2, tmp_a, tmp_m1, tmp_m2);
        result.x2m.iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_m1);
    }
    result.fil = gb_[s][cp.i2].fil + cp.m2.w_may();

    size_t index;
    index = 0;
    while (index < result.x1.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], result.x1.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.x1.data[index], gb_[s][gb_index].x1.data[0]);
            if (result.valid_x2m() && result.fil == Filtr(result.x1.data[index]))
                result.x2m.iaddmulMay(m, gb_[s][gb_index].x2m, tmp_m1);
            result.x1.iaddmul(m, gb_[s][gb_index].x1, tmp_a, tmp_m1, tmp_m2);
            result.x2.iaddmul(m, gb_[s][gb_index].x2, tmp_a, tmp_m1, tmp_m2);
        }
        else
            ++index;
    }
    size_t sp1 = s + 1;
    index = 0;
    while (index < result.x2.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[sp1], indices_[sp1], result.x2.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.x2.data[index], gb_[sp1][gb_index].x1.data[0]);
            result.x2.iaddmul(m, gb_[sp1][gb_index].x1, tmp_a, tmp_m1, tmp_m2);
        }
        else
            ++index;
    }
    result.x2m = gb_x2m_.Reduce(std::move(result.x2m), s);

    return result;
}

Mod GroebnerMRes::ReduceX2m(const CriMilnor& cp, size_t s) const
{
    Mod result;

    Mod tmp_m2;
    tmp_m2.data.reserve(100);

    if (cp.i1 >= 0)
        result.iaddmulMay(cp.m1, gb_[s][cp.i1].x2m, tmp_m2).iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_m2);
    else
        result.iaddmulMay(cp.m2, gb_[s][cp.i2].x2m, tmp_m2);

    return gb_x2m_.Reduce(result, s);
}

class DbSteenrod : public myio::DbAlg
{
    using Statement = myio::Statement;

public:
    DbSteenrod() = default;
    explicit DbSteenrod(const std::string& filename) : myio::DbAlg(filename) {}

public:
    void create_generators(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, name TEXT UNIQUE, diff BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (id INTEGER PRIMARY KEY, gs BLOB, gsplus BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_generators_and_delete(const std::string& table_prefix) const
    {
        create_generators(table_prefix);
        delete_from(table_prefix + "_generators");
    }
    void create_relations_and_delete(const std::string& table_prefix) const
    {
        create_relations(table_prefix);
        delete_from(table_prefix + "_relations");
    }

    void save_generators(const std::string& table_prefix, const Mod1d& kernel, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (diff, s, t) VALUES (?1, ?2, ?3);");

        for (auto& k : kernel) {
            stmt.bind_blob(1, k.data.data(), (int)k.data.size() * sizeof(MMod));
            stmt.bind_int(2, s);
            stmt.bind_int(3, t);
            stmt.step_and_reset();
        }
    }
    void save_relations(const std::string& table_prefix, const DataMRes1d& rels, const array& basis_degrees, int s) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (gs, gsplus, s, t) VALUES (?1, ?2, ?3, ?4);");

        for (auto& rel : rels) {
            stmt.bind_blob(1, rel.x1.data.data(), (int)rel.x1.data.size() * sizeof(MMod));
            stmt.bind_blob(2, rel.x2.data.data(), (int)rel.x2.data.size() * sizeof(MMod));
            stmt.bind_int(3, s);
            int t = rel.x1.GetLead().deg_m() + basis_degrees[rel.x1.GetLead().v()];
            stmt.bind_int(4, t);
            stmt.step_and_reset();
        }
    }

    array2d load_basis_degrees(const std::string& table_prefix) const
    {
        array2d basis_degrees;
        Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            int s = stmt.column_int(0), t = stmt.column_int(1);
            if (basis_degrees.size() <= (size_t)s)
                basis_degrees.resize(size_t(s + 1));
            basis_degrees[s].push_back(t);
        }
        return basis_degrees;
    }
    DataMRes2d load_data(const std::string& table_prefix) const
    {
        DataMRes2d data;
        Statement stmt(*this, "SELECT gs, gsplus, s FROM " + table_prefix + "_relations ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            DataMRes x;
            const void* x1_data = stmt.column_blob(0);
            const void* x2_data = stmt.column_blob(1);
            int x1_bytes = stmt.column_blob_size(0);
            int x2_bytes = stmt.column_blob_size(1);
            size_t x1_size = (size_t)x1_bytes / sizeof(MMod);
            size_t x2_size = (size_t)x2_bytes / sizeof(MMod);
            x.x1.data.resize(x1_size);
            x.x2.data.resize(x2_size);
            memcpy(x.x1.data.data(), x1_data, x1_bytes);
            memcpy(x.x2.data.data(), x2_data, x2_bytes);

            int s = stmt.column_int(2);
            if (data.size() <= (size_t)s)
                data.resize(size_t(s + 1));
            data[s].push_back(std::move(x));
        }
        return data;
    }
};

size_t AddRelsMRes(GroebnerMRes& gb, const Mod1d& rels, int deg)
{
    size_t dim = 0;
    int deg_max = gb.t_trunc();
    if (deg > deg_max)
        throw MyException(0xb2474e19U, "deg is bigger than the truncation degree.");

    /*DbSteenrod db("AdamsE2.db");
    db.create_generators("SteenrodMRes");
    db.create_relations("SteenrodMRes");*/

    /* Group `rels` by degree */
    std::map<int, array> rels_graded;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            const auto& lead = rels[i].GetLead();
            int t = lead.deg_m();
            if (t <= deg)
                rels_graded[t].push_back((int)i);
        }
    }

    int t_begin = gb.t_begin();
    for (int t = t_begin; t <= deg; ++t) {
        std::cout << "t=" << t << "               \n";
        for (int s = t - 2; s >= -1; --s) {
            std::cout << "  s=" << s + 2 << "                \r";

            /* Populate `data_tmp` */
            gb.resize_gb(size_t(s + 2));
            DataMRes1d data_tmp;
            Mod1d x2m_st_tmp;
            if (s == -1) {
                auto p_rels_d = rels_graded.find(t);
                if (p_rels_d != rels_graded.end()) {
                    data_tmp.resize(p_rels_d->second.size());
                    auto& rels_d = p_rels_d->second;
                    ut::for_each_seq((int)data_tmp.size(), [&](size_t i) { data_tmp[i] = DataMRes(Mod(), rels[rels_d[i]], Mod()); });
                }
            }
            else {
                CriMilnor1d cris_st = gb.Criticals(s, t);
                if (!cris_st.empty()) {
                    data_tmp.resize(cris_st.size());
                    ut::for_each_seq((int)data_tmp.size(), [&](size_t i) { data_tmp[i] = gb.Reduce(cris_st[i], s); });
                }
            }
            if (data_tmp.empty())
                continue;

            /* Triangulate these relations */
            DataMRes1d data_st;
            Mod1d kernel_sp1_tmp;
            for (size_t i = 0; i < data_tmp.size(); ++i) {
                for (size_t j = 0; j < data_st.size(); ++j)
                    if (std::binary_search(data_tmp[i].x1.data.begin(), data_tmp[i].x1.data.end(), data_st[j].x1.GetLead()))
                        data_tmp[i] += data_st[j];
                if (data_tmp[i].x1) {
                    /* Determine if x2m aligns with x1 */
                    if (!data_tmp[i].valid_x2m()) {
                        x2m_st_tmp.push_back(gb.new_gen_x2m(s, t));
                        std::swap(x2m_st_tmp.back(), data_tmp[i].x2m);
                        data_tmp[i].fil = Filtr(data_tmp[i].x1.GetLead());
                    }
                    data_st.push_back(std::move(data_tmp[i]));
                    bench::Counter(0);
                }
                else {
                    if (data_tmp[i].x2)
                        kernel_sp1_tmp.push_back(std::move(data_tmp[i].x2));
                    else
                        bench::Counter(1);
                    if (data_tmp[i].x2m)
                        x2m_st_tmp.push_back(std::move(data_tmp[i].x2m));
                }
            }
            Mod1d kernel_sp1;
            for (size_t i = 0; i < kernel_sp1_tmp.size(); ++i) {
                for (size_t j = 0; j < kernel_sp1.size(); ++j)
                    if (std::binary_search(kernel_sp1_tmp[i].data.begin(), kernel_sp1_tmp[i].data.end(), kernel_sp1[j].GetLead()))
                        kernel_sp1_tmp[i] += kernel_sp1[j];
                if (kernel_sp1_tmp[i]) {
                    kernel_sp1.push_back(std::move(kernel_sp1_tmp[i]));
                    bench::Counter(2);
                }
                else
                    bench::Counter(3);
            }
            Mod1d x2m_st;
            for (size_t i = 0; i < x2m_st_tmp.size(); ++i) {
                for (size_t j = 0; j < x2m_st.size(); ++j)
                    if (std::binary_search(x2m_st_tmp[i].data.begin(), x2m_st_tmp[i].data.end(), x2m_st[j].GetLead()))
                        x2m_st_tmp[i] += x2m_st[j];
                if (x2m_st_tmp[i])
                    x2m_st.push_back(std::move(x2m_st_tmp[i]));
            }

            /* Add these relations */
            for (size_t i = 0; i < data_st.size(); ++i)
                gb.push_back(data_st[i], s);
            DataMRes1d rels_splus;
            for (size_t i = 0; i < kernel_sp1.size(); ++i)
                gb.push_back_kernel(kernel_sp1[i], rels_splus, size_t(s + 1), t);
            for (size_t i = 0; i < x2m_st.size(); ++i)
                gb.push_back_x2m(x2m_st[i], s);

            /*db.begin_transaction();
            db.save_generators("SteenrodMRes", kernel_sp1, (int)s + 2, t);
            if (s >= 0)
                db.save_relations("SteenrodMRes", data_st, gb.basis_degs()[s], s);
            db.save_relations("SteenrodMRes", rels_splus, gb.basis_degs()[size_t(s + 1)], s + 1);
            db.end_transaction();*/

            if (size_t size_k = kernel_sp1.size()) {
                std::cout << "  s=" << s + 2 << " dim=" << size_k << '\n';
                dim += size_k;
            }
        }
    }
    return dim;
}

GroebnerMRes GroebnerMRes::load(const std::string& filename, int t_trunc)
{
    DbSteenrod db("AdamsE2.db");
    db.create_generators_and_delete("SteenrodMRes");  ////
    db.create_relations_and_delete("SteenrodMRes");   ////
    array2d basis_degrees = db.load_basis_degrees("SteenrodMRes");
    DataMRes2d data = db.load_data("SteenrodMRes");
    return GroebnerMRes(t_trunc, std::move(data), std::move(basis_degrees), Mod2d(), array2d());
}

}  // namespace steenrod