#include "groebner_steenrod.h"
#include "algebras/benchmark.h"  ////
#include "algebras/dbalg.h"
#include <cstring>
//#include <immintrin.h>

namespace steenrod {

void CPMilnors::Minimize(const MMod1d& leads, int t)
{
    if (buffer_redundent_pairs_.empty() || buffer_redundent_pairs_.begin()->first != t)
        return;
#ifndef NDEBUG
    if (buffer_min_pairs_.empty() || buffer_min_pairs_.begin()->first != t)
        throw MyException(0xb1b36dceU, "buffer_min_pairs_ is expected to start with t");
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
                auto p = std::find_if(b_min_pairs_t[j].begin(), b_min_pairs_t[j].end(), [&m2](const CPMilnor& c) { return c.m2 == m2; });
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
void CPMilnors::AddToBuffers(const MMod1d& leads, MMod mon, int t_v)
{
    size_t lead_size = leads.size();
    std::vector<std::pair<int, CPMilnor>> new_pairs(lead_size);

    /* Populate `new_pairs` */
    for (size_t i = 0; i < leads.size(); ++i) {
        new_pairs[i].first = -1;
        if (leads[i].v_raw() == mon.v_raw()) {
            int d_pair = lcmLF(leads[i].m_no_weight(), mon.m_no_weight()).deg() + t_v;
            if (d_pair <= t_trunc_) {
                new_pairs[i].first = d_pair;
                CPMilnor::SetFromLM(new_pairs[i].second, leads[i].m(), mon.m(), (int)i, (int)lead_size);
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
            buffer_singles_[t].push_back(CPMilnor::Single(m, (int)lead_size));
    }
}

void CPMilnors::init(const MMod1d& leads, const array& basis_degrees)
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

            /* Populate `rels_tmp` */
            gb.resize_gb(size_t(s + 2));
            DataMRes1d rels_tmp;
            if (s == -1) {
                auto p_rels_d = rels_graded.find(t);
                if (p_rels_d != rels_graded.end()) {
                    rels_tmp.resize(p_rels_d->second.size());
                    auto& rels_d = p_rels_d->second;
                    ut::for_each_seq((int)rels_tmp.size(), [&](size_t i) { rels_tmp[i] = DataMRes(Mod(), rels[rels_d[i]]); });
                }
            }
            else {
                gb.MinimizePairs(s, t);
                CPMilnor1d pairs_st = gb.pairs(s, t);
                if (!pairs_st.empty()) {
                    rels_tmp.resize(pairs_st.size());
                    ut::for_each_seq((int)rels_tmp.size(), [&](size_t i) { rels_tmp[i] = gb.Reduce(pairs_st[i], s); });
                }
            }
            if (rels_tmp.empty())
                continue;

            /* Triangulate these relations */
            DataMRes1d rels_st;
            Mod1d kernel_sp1_tmp;
            for (size_t i = 0; i < rels_tmp.size(); ++i) {
                for (size_t j = 0; j < rels_st.size(); ++j)
                    if (std::binary_search(rels_tmp[i].x1.data.begin(), rels_tmp[i].x1.data.end(), rels_st[j].x1.GetLead()))
                        rels_tmp[i] += rels_st[j];
                if (rels_tmp[i].x1) {
                    rels_st.push_back(std::move(rels_tmp[i]));
                    bench::Counter(0);
                }
                else if (rels_tmp[i].x2)
                    kernel_sp1_tmp.push_back(std::move(rels_tmp[i].x2));
                else
                    bench::Counter(1);
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

            /* Add these relations */
            for (size_t i = 0; i < rels_st.size(); ++i)
                gb.push_back(rels_st[i], s);
            DataMRes1d rels_splus;
            for (size_t i = 0; i < kernel_sp1.size(); ++i) {
                gb.push_back_kernel(kernel_sp1[i], rels_splus, (int)s + 1, t);
            }

            /*db.begin_transaction();
            db.save_generators("SteenrodMRes", kernel_sp1, (int)s + 2, t);
            if (s >= 0)
                db.save_relations("SteenrodMRes", rels_st, gb.basis_degs()[s], s);
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
    return GroebnerMRes(t_trunc, std::move(data), std::move(basis_degrees));
}

}  // namespace steenrod