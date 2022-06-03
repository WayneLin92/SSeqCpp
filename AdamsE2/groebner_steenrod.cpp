#include "groebner_steenrod.h"
#include "algebras/benchmark.h"
#include "algebras/database.h"
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

void CriMilnors::init(const MMod1d& leads, const array& basis_degrees, int t_min_buffer)
{
    MMod1d tmp_leads;
    for (auto& mon : leads) {
        int t_v = basis_degrees[mon.v()];
        AddToBuffers(tmp_leads, mon, t_v);
        tmp_leads.push_back(mon);
    }

    for (auto it = buffer_min_pairs_.begin(); it != buffer_min_pairs_.end();) {
        if (it->first < t_min_buffer)
            it = buffer_min_pairs_.erase(it);
        else
            ++it;
    }
    for (auto it = buffer_redundent_pairs_.begin(); it != buffer_redundent_pairs_.end();) {
        if (it->first < t_min_buffer)
            it = buffer_redundent_pairs_.erase(it);
        else
            ++it;
    }
    for (auto it = buffer_singles_.begin(); it != buffer_singles_.end();) {
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

GroebnerX2m::GroebnerX2m(int t_trunc, Mod2d data, array2d basis_degrees, int latest_s, int latest_t) : t_trunc_(t_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees))
{
    leads_.resize(gb_.size());
    indices_.resize(gb_.size());

    for (size_t s = 0; s < gb_.size(); ++s) {
        for (int j = 0; j < (int)gb_[s].size(); ++j) {
            leads_[s].push_back(gb_[s][j].GetLead());
            indices_[s][gb_[s][j].GetLead().v_raw()].push_back(j);
        }
    }

    for (size_t s = 0; s < gb_.size(); ++s) {
        criticals_.push_back(CriMilnors(t_trunc_));
        criticals_.back().init(leads_[s], basis_degrees_[s], s >= latest_s ? latest_t + 1 : latest_t);
    }
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

Mod1d GroebnerX2m::AddRels(size_t s, int t)
{
    Mod tmp_Mod;

    /* Populate `rels_tmp` */
    resize(s + 1);
    Mod1d rels_tmp;

    criticals_[s].Minimize(leads_[s], t);
    CriMilnor1d pairs_st = criticals_[s].Criticals(t);
    if (!pairs_st.empty()) {
        rels_tmp.resize(pairs_st.size());
        ut::for_each_seq((int)rels_tmp.size(), [&](size_t i) { rels_tmp[i] = Reduce(pairs_st[i], s); });
    }

    Mod1d rels_st;
    if (rels_tmp.empty())
        return rels_st;

    /* Triangulate these relations */
    for (size_t i = 0; i < rels_tmp.size(); ++i) {
        steenrod::Reduce(rels_tmp[i], rels_st, tmp_Mod);
        if (rels_tmp[i])
            rels_st.push_back(std::move(rels_tmp[i]));
    }

    /* Add these relations */
    for (size_t i = 0; i < rels_st.size(); ++i)
        push_back(rels_st[i], s);
    return rels_st;
}

/********************************************************
 *                    class GroebnerMRes
 ********************************************************/

GroebnerMRes::GroebnerMRes(int t_trunc, DataMRes2d data, array2d basis_degrees, Mod2d data_x2m, array2d basis_degrees_x2m, int latest_s, int latest_t)
    : t_trunc_(t_trunc), gb_(std::move(data)), basis_degrees_(std::move(basis_degrees)), gb_x2m_(t_trunc, std::move(data_x2m), std::move(basis_degrees_x2m), latest_s, latest_t)
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
        criticals_.back().init(leads_[s], basis_degrees_[s], s >= latest_s ? latest_t + 1 : latest_t);
    }
}

CriMilnor1d GroebnerMRes::Criticals(size_t s, int t, Mod1d& rels_x2m)
{
    rels_x2m = gb_x2m_.AddRels(s, t);
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
    Mod tmp_Mod;
    for (auto& cp : result) {
        Mod x2m = ReduceX2m(cp, s);
        steenrod::Reduce(x2m, x2ms, tmp_Mod);
        if (x2m)
            x2ms.push_back(std::move(x2m));
        else {
#ifndef MYDEPLOY
            bench::Counter(1);
#endif
            cp.i2 = -1;
        }
    }
    ut::RemoveIf(result, [](const CriMilnor& cp) { return cp.i2 == -1; });
    return result;
    return cris;
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

Mod GroebnerMRes::Reduce(Mod x, size_t s) const
{
    size_t index;
    index = 0;
    Milnor tmp_a;
    Mod tmp_x1, tmp_x2;
    while (index < x.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], x.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(x.data[index], gb_[s][gb_index].x1.data[0]);
            x.iaddmul(m, gb_[s][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
            x += m * gb_[s][gb_index].x1;
        }
        else
            ++index;
    }

    return x;
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

class DbSteenrod : public myio::Database
{
    using Statement = myio::Statement;

private:
    std::future<void> f_;

public:
    DbSteenrod() = default;
    explicit DbSteenrod(const std::string& filename) : Database(filename) {}

public:
    void create_generators(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_generators (id INTEGER PRIMARY KEY, diff BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_generators_x2m(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_X2m_generators (id INTEGER PRIMARY KEY, s SMALLINT, t SMALLINT);");
    }
    void create_relations(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_relations (id INTEGER PRIMARY KEY, x1 BLOB, x2 BLOB, x2m BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_relations_x2m(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_X2m_relations (id INTEGER PRIMARY KEY, x2m BLOB, s SMALLINT, t SMALLINT);");
    }
    void create_time(const std::string& table_prefix) const
    {
        execute_cmd("CREATE TABLE IF NOT EXISTS " + table_prefix + "_time (id INTEGER PRIMARY KEY, s SMALLINT, t SMALLINT, time REAL, UNIQUE(s, t));");
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
    void create_generators_x2m_and_delete(const std::string& table_prefix) const
    {
        create_generators_x2m(table_prefix);
        delete_from(table_prefix + "_X2m_generators");
    }
    void create_relations_x2m_and_delete(const std::string& table_prefix) const
    {
        create_relations_x2m(table_prefix);
        delete_from(table_prefix + "_X2m_relations");
    }
    void create_time_and_delete(const std::string& table_prefix) const
    {
        create_time(table_prefix);
        delete_from(table_prefix + "_time");
    }

    void save_generators(const std::string& table_prefix, const DataMRes1d& kernel, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_generators (diff, s, t) VALUES (?1, ?2, ?3);");

        for (auto& k : kernel) {
            stmt.bind_blob(1, k.x1.data.data(), int(k.x1.data.size() * sizeof(MMod)));
            stmt.bind_int(2, s);
            stmt.bind_int(3, t);
            stmt.step_and_reset();
        }
    }

    /* insert v_{0, 0} */
    void save_v0(const std::string& table_prefix) const
    {
        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_generators (id, diff, s, t) VALUES (?1, ?2, ?3, ?4);");
        Mod v0;
        v0.data.reserve(1);
        stmt.bind_int(1, 1); //TODO: should start with zero
        stmt.bind_blob(2, v0.data.data(), int(v0.data.size() * sizeof(MMod)));
        stmt.bind_int(3, 0);
        stmt.bind_int(4, 0);
        stmt.step_and_reset();
    }

    void save_relations(const std::string& table_prefix, const DataMRes1d& rels, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_relations (x1, x2, x2m, s, t) VALUES (?1, ?2, ?3, ?4, ?5);");

        for (auto& rel : rels) {
            stmt.bind_blob(1, rel.x1.data.data(), (int)rel.x1.data.size() * sizeof(MMod));
            stmt.bind_blob(2, rel.x2.data.data(), (int)rel.x2.data.size() * sizeof(MMod));
            stmt.bind_blob(3, rel.x2m.data.data(), (int)rel.x2m.data.size() * sizeof(MMod));
            stmt.bind_int(4, s);
            stmt.bind_int(5, t);
            stmt.step_and_reset();
        }
    }

    void save_generators_x2m(const std::string& table_prefix, int num, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_X2m_generators (s, t) VALUES (?1, ?2);");

        for (int k = 0; k < num; ++k) {
            stmt.bind_int(1, s);
            stmt.bind_int(2, t);
            stmt.step_and_reset();
        }
    }

    void save_relations_x2m(const std::string& table_prefix, const Mod1d& rels, int s, int t) const
    {
        Statement stmt(*this, "INSERT INTO " + table_prefix + "_X2m_relations (x2m, s, t) VALUES (?1, ?2, ?3);");

        for (auto& rel : rels) {
            stmt.bind_blob(1, rel.data.data(), (int)rel.data.size() * sizeof(MMod));
            stmt.bind_int(2, s);
            stmt.bind_int(3, t);
            stmt.step_and_reset();
        }
    }

    void save_time(const std::string& table_prefix, int s, int t, double time)
    {
        Statement stmt(*this, "INSERT OR IGNORE INTO " + table_prefix + "_time (s, t, time) VALUES (?1, ?2, ?3);");

        stmt.bind_int(1, s);
        stmt.bind_int(2, t);
        stmt.bind_double(3, time);
        stmt.step_and_reset();
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

    array2d load_basis_degrees_x2m(const std::string& table_prefix) const
    {
        return load_basis_degrees(table_prefix + "_X2m");
    }

    DataMRes2d load_data(const std::string& table_prefix) const
    {
        DataMRes2d data;
        Statement stmt(*this, "SELECT x1, x2, x2m, s FROM " + table_prefix + "_relations ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            DataMRes x;
            const void* x1_data = stmt.column_blob(0);
            const void* x2_data = stmt.column_blob(1);
            const void* x2m_data = stmt.column_blob(2);
            int x1_bytes = stmt.column_blob_size(0);
            int x2_bytes = stmt.column_blob_size(1);
            int x2m_bytes = stmt.column_blob_size(2);
            size_t x1_size = (size_t)x1_bytes / sizeof(MMod);
            size_t x2_size = (size_t)x2_bytes / sizeof(MMod);
            size_t x2m_size = (size_t)x2m_bytes / sizeof(MMod);
            x.x1.data.resize(x1_size);
            x.x2.data.resize(x2_size);
            x.x2m.data.resize(x2m_size);
            memcpy(x.x1.data.data(), x1_data, x1_bytes);
            memcpy(x.x2.data.data(), x2_data, x2_bytes);
            memcpy(x.x2m.data.data(), x2m_data, x2m_bytes);

            x.fil = Filtr(x.x1.GetLead());

            size_t s = (size_t)stmt.column_int(3);
            if (data.size() <= s)
                data.resize(s + 1);
            data[s].push_back(std::move(x));
        }
        return data;
    }

    Mod2d load_data_x2m(const std::string& table_prefix) const
    {
        Mod2d data;
        Statement stmt(*this, "SELECT x2m, s FROM " + table_prefix + "_X2m_relations ORDER BY id;");
        while (stmt.step() == MYSQLITE_ROW) {
            Mod x;
            const void* x1_data = stmt.column_blob(0);
            int x1_bytes = stmt.column_blob_size(0);
            size_t x1_size = (size_t)x1_bytes / sizeof(MMod);
            x.data.resize(x1_size);
            memcpy(x.data.data(), x1_data, x1_bytes);

            size_t s = (size_t)stmt.column_int(1);
            if (data.size() <= s)
                data.resize(s + 1);
            data[s].push_back(std::move(x));
        }
        return data;
    }

    void latest_st(const std::string& table_prefix, int& s, int& t) const
    {
        t = get_int("SELECT coalesce(max(t), 0) FROM " + table_prefix + "_relations;");
        s = get_int("SELECT coalesce(min(s), 1000) FROM " + table_prefix + "_relations WHERE t=" + std::to_string(t) + ";");
        int t1 = get_int("SELECT coalesce(max(t), 0) FROM " + table_prefix + "_X2m_relations;");
        int s1 = get_int("SELECT coalesce(min(s), 1000) FROM " + table_prefix + "_X2m_relations WHERE t=" + std::to_string(t1) + ";");
        if (t1 > t || (t1 == t && s1 < s)) {
            s = s1;
            t = t1;
        }
    }

    void save(const DataMRes1d& data_sp1t, const DataMRes1d& data_st, const Mod1d& x2m_st1, const Mod1d& x2m_st, int num_x2m_sp1, int num_x2m_s, double time, int s, int t)
    {
        if (f_.valid())
            f_.wait();
        f_ = std::async(std::launch::async, [this, data_sp1t, data_st, x2m_st1, x2m_st, num_x2m_sp1, num_x2m_s, time, s, t]() {
            begin_transaction();
            save_generators("SteenrodMRes", data_sp1t, s + 2, t);
            save_relations("SteenrodMRes", data_sp1t, s + 1, t);
            save_generators_x2m("SteenrodMRes", num_x2m_sp1, s + 1, t);
            if (s >= 0) {
                save_relations("SteenrodMRes", data_st, s, t);
                save_generators_x2m("SteenrodMRes", num_x2m_s, s, t);
                save_relations_x2m("SteenrodMRes", x2m_st1, s, t);
                save_relations_x2m("SteenrodMRes", x2m_st, s, t);
            }
            save_time("SteenrodMRes", s + 2, t, time);
            end_transaction();
        });
    }
};

void AddRelsMRes(GroebnerMRes& gb, const Mod1d& rels, int deg)
{
    int deg_max = gb.t_trunc();
    if (deg > deg_max)
        throw MyException(0xb2474e19U, "deg is bigger than the truncation degree.");
    Mod tmp_Mod;

#ifdef MYDEPLOY
    DbSteenrod db("AdamsE2.db");
    db.create_generators("SteenrodMRes");
    db.create_relations("SteenrodMRes");
    db.create_generators_x2m("SteenrodMRes");
    db.create_relations_x2m("SteenrodMRes");
    db.create_time("SteenrodMRes");
    db.save_v0("SteenrodMRes");

    bench::Timer timer;
    timer.SuppressPrint();
#endif

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
        std::cout << "t=" << t << "               " << std::endl;
        for (int s = t - 2; s >= -1; --s) {
            std::cout << "  s=" << s + 2 << "                \r";
            size_t sp1 = size_t(s + 1);
            size_t sp2 = size_t(s + 2);

            /* Populate `data_tmp` */
            gb.resize_gb(sp2);
            size_t old_size_x2m_s = 0;
            size_t old_size_x2m_sp1 = gb.basis_degrees_x2m(sp1).size();

            DataMRes1d data_tmp;
            Mod1d x2m_st_tmp;
            Mod1d x2m_st1;
            if (s == -1) {
                auto p_rels_d = rels_graded.find(t);
                if (p_rels_d != rels_graded.end()) {
                    data_tmp.resize(p_rels_d->second.size());
                    auto& rels_d = p_rels_d->second;
                    ut::for_each_seq((int)data_tmp.size(), [&](size_t i) { data_tmp[i] = DataMRes(Mod(), gb.Reduce(rels[rels_d[i]], s + 1), Mod()); });
                }
            }
            else {
                old_size_x2m_s = gb.basis_degrees_x2m(s).size();
                CriMilnor1d cris_st = gb.Criticals(s, t, x2m_st1);
                if (!cris_st.empty()) {
                    data_tmp.resize(cris_st.size());
#ifdef MYDEPLOY
                    ut::for_each_par((int)data_tmp.size(), [&](size_t i) { data_tmp[i] = gb.Reduce(cris_st[i], s); });
#else
                    ut::for_each_seq((int)data_tmp.size(), [&](size_t i) { data_tmp[i] = gb.Reduce(cris_st[i], s); });
#endif
                }
            }

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
                }
                else {
                    if (data_tmp[i].x2)
                        kernel_sp1_tmp.push_back(std::move(data_tmp[i].x2));
#ifndef MYDEPLOY
                    else
                        bench::Counter(0);
#endif
                    if (data_tmp[i].x2m)
                        x2m_st_tmp.push_back(std::move(data_tmp[i].x2m));
                }
            }
            DataMRes1d data_sp1t;
            for (size_t i = 0; i < kernel_sp1_tmp.size(); ++i) {
                Reduce(kernel_sp1_tmp[i], data_sp1t, tmp_Mod);
                if (kernel_sp1_tmp[i])
                    data_sp1t.push_back(DataMRes(std::move(kernel_sp1_tmp[i]), gb.new_gen(sp2, t), gb.new_gen_x2m(sp1, t)));
#ifndef MYDEPLOY
                else
                    bench::Counter(0);
#endif
            }
            Mod1d x2m_st;
            for (size_t i = 0; i < x2m_st_tmp.size(); ++i) {
                Reduce(x2m_st_tmp[i], x2m_st, tmp_Mod);
                if (x2m_st_tmp[i])
                    x2m_st.push_back(std::move(x2m_st_tmp[i]));
            }

#ifdef MYDEPLOY
            double time = timer.Elapsed();
            int num_x2m_sp1 = int(gb.basis_degrees_x2m(sp1).size() - old_size_x2m_sp1);
            int num_x2m_s = s >= 0 ? int(gb.basis_degrees_x2m(s).size() - old_size_x2m_s) : 0;
            db.save(data_sp1t, data_st, x2m_st1, x2m_st, num_x2m_sp1, num_x2m_s, time, s, t);
            timer.Reset();
#endif

            for (size_t i = 0; i < data_st.size(); ++i)
                gb.push_back(data_st[i], s);
            for (size_t i = 0; i < data_sp1t.size(); ++i)
                gb.push_back(data_sp1t[i], sp1);
            for (size_t i = 0; i < x2m_st.size(); ++i)
                gb.push_back_x2m(x2m_st[i], s);

            if (size_t size_k = data_sp1t.size())
#ifdef MYDEPLOY
                std::cout << "  s=" << s + 2 << " dim=" << size_k << ' ' << time << std::endl;
#else
                std::cout << "  s=" << s + 2 << " dim=" << size_k << std::endl;
#endif
        }
    }
}

GroebnerMRes GroebnerMRes::load(const std::string& filename, int t_trunc)
{
    DbSteenrod db(filename);
    db.create_generators("SteenrodMRes");
    db.create_relations("SteenrodMRes");
    db.create_generators_x2m("SteenrodMRes");
    db.create_relations_x2m("SteenrodMRes");
    DataMRes2d data = db.load_data("SteenrodMRes");
    array2d basis_degrees = db.load_basis_degrees("SteenrodMRes");
    Mod2d data_x2m = db.load_data_x2m("SteenrodMRes");
    array2d basis_degrees_x2m = db.load_basis_degrees_x2m("SteenrodMRes");
    int latest_s = 0, latest_t = 0;
    db.latest_st("SteenrodMRes", latest_s, latest_t);
    return GroebnerMRes(t_trunc, std::move(data), std::move(basis_degrees), std::move(data_x2m), std::move(basis_degrees_x2m), latest_s, latest_t);
}

void GroebnerMRes::reset(const std::string& filename)
{
    DbSteenrod db(filename);
    db.create_generators_and_delete("SteenrodMRes");
    db.create_relations_and_delete("SteenrodMRes");
    db.create_generators_x2m_and_delete("SteenrodMRes");
    db.create_relations_x2m_and_delete("SteenrodMRes");
    db.create_time_and_delete("SteenrodMRes");
}

}  // namespace steenrod