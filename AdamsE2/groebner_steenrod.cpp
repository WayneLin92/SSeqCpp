#include "groebner_steenrod.h"
#include "algebras/benchmark.h"  ////
#include "algebras/dbalg.h"
#include <cstring>
//#include <immintrin.h>  

namespace steenrod {

void GbCriPairsMRes::Minimize(const MMod1d& leads, AdamsDeg st)
{
    /* Add to the Groebner basis of critical pairs */
    auto& b_min_pairs_st = buffer_min_pairs_.at(st);
    if (pairs_[st.s].size() < b_min_pairs_st.size())
        pairs_[st.s].resize(b_min_pairs_st.size());
    array old_sizes(b_min_pairs_st.size());
    for (int j = 0; j < (int)b_min_pairs_st.size(); ++j) {
        old_sizes[j] = (int)pairs_[st.s][j].size();
        for (const auto& pair : b_min_pairs_st[j])
            pairs_[st.s][j].push_back(pair);
    }

    /* Minimize `buffer_min_pairs_` */
    for (uint64_t ij : buffer_redundent_pairs_[st]) {
        int i, j;
        ut::GetPair(ij, i, j);
        while (j != -1) {
            MMilnor gcd = gcdLF(leads[i].m(), leads[j].m());
            MMilnor m2 = divLF(leads[i].m(), gcd);
            if (j < (int)b_min_pairs_st.size()) {
                auto p = std::find_if(b_min_pairs_st[j].begin(), b_min_pairs_st[j].end(), [&m2](const CriPairMRes& c) { return c.m2 == m2; });
                /* Remove it from `buffer_min_pairs_` */
                if (p != b_min_pairs_st[j].end()) {
                    p->i2 = -1;
                    break;
                }
            }

            /* Reduce (i, j) */
            auto c = pairs_[st.s][j].begin();
            auto end = pairs_[st.s][j].end();
            for (; c < end; ++c) {
                if (divisibleLF(c->m2, m2)) {
                    MMilnor m1 = divLF(leads[j].m(), gcd);
                    if (gcdLF(c->m1, m1)) {
                        j = -1;
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
                for (size_t k = 0; k < pairs_[st.s][j].size(); ++k)
                    std::cout << "pairs_[j][k].m2" << Milnor(pairs_[st.s][j][k].m2) << '\n';
                throw MyException(0xfa5db14U, "Should not happen because pairs_ is groebner");
            }
#endif
        }
    }

    /* Delete `bufbuffer_redundent_pairs_[t]` */
    buffer_redundent_pairs_.erase(st);
}

/* Propogate `buffer_redundent_pairs_` and `buffer_min_pairs_`.
** `buffer_min_pairs_` will become a Groebner basis at this stage.
*/
void GbCriPairsMRes::AddToBuffers(const MMod1d& leads, MMod mon, int d_mon_base, int s)
{
    std::vector<std::pair<int, CriPairMRes>> new_pairs(leads.size());
    size_t lead_size = leads.size();
    ut::Range range(0, (int)leads.size());

    for (size_t i = 0; i < leads.size(); ++i) {
        new_pairs[i].first = -1;
        if (leads[i].v() == mon.v()) {
            int d_pair = lcmLF(leads[i].m(), mon.m()).deg() + d_mon_base;
            if (d_pair <= deg_trunc_) {
                new_pairs[i].first = d_pair;
                CriPairMRes::SetFromLM(new_pairs[i].second, leads[i].m(), mon.m(), (int)i, (int)lead_size);
            }
        }
    }

    /* Compute sigma_j */
    int d_mon = mon.m().deg() + d_mon_base;
    for (int i : mon.m()) {
        MMilnor m = MMilnor::FromIndex(i);
        AdamsDeg st{s, d_mon + m.deg()};
        if (st.t <= deg_trunc_)
            buffer_singles_[st].push_back(CriPairMRes::Single(m, (int)lead_size));
    }

    /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
    for (size_t j = 1; j < new_pairs.size(); ++j) {
        if (new_pairs[j].first != -1) {
            for (size_t i = 0; i < j; ++i) {
                if (new_pairs[i].first != -1) {
                    if (divisibleLF(new_pairs[i].second.m2, new_pairs[j].second.m2)) {
                        new_pairs[j].first = -1;
                        break;
                    }
                    else if (divisibleLF(new_pairs[j].second.m2, new_pairs[i].second.m2))
                        new_pairs[i].first = -1;
                    else if (!gcdLF(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                        int dij = lcmLF(leads[i].m(), leads[j].m()).deg() + d_mon_base;
                        if (dij <= deg_trunc_)
                            buffer_redundent_pairs_[AdamsDeg{s, dij}].insert(ut::BindPair((uint32_t)i, (uint32_t)j));
                    }
                }
            }
        }
    }
    for (size_t i = 0; i < new_pairs.size(); ++i) {
        if (new_pairs[i].first != -1) {
            buffer_min_pairs_[AdamsDeg{s, new_pairs[i].first}].resize(lead_size + 1);
            buffer_min_pairs_[AdamsDeg{s, new_pairs[i].first}][lead_size].push_back(new_pairs[i].second);
        }
    }
}

void GbCriPairsMRes::init(const MMod2d& leads, const array2d& basis_degrees)
{
    int t_max = (int)basis_degrees.size() - 1;
    int s_sep = 0;
    while (s_sep < t_max) {
        if (basis_degrees[s_sep].back() == t_max)
            break;
        ++s_sep;
    }

    for (int s = 0; (size_t)s < leads.size(); ++s) {
        int t_min_buffer = s + 2 >= s_sep ? t_max + 1 : t_max;
        for (size_t k = 0; k < leads[s].size(); ++k) {
            auto& mon = leads[s][k];
            int d_mon_base = basis_degrees[s][leads[s][k].v()];
            int d_mon = mon.m().deg() + d_mon_base;

            std::vector<std::pair<int, CriPairMRes>> new_pairs(k);
            ut::Range range(0, (int)leads.size());

            for (size_t i = 0; i < k; ++i) {
                if (leads[s][i].v() == mon.v()) {
                    int d_pair = lcmLF(leads[s][i].m(), mon.m()).deg() + d_mon_base;
                    if (d_pair <= deg_trunc_) {
                        new_pairs[i].first = d_pair;
                        // CriPairMRes::SetFromLM(new_pairs[i].second, leads[s][i].m(), mon.m(), (int)i, (int)k);
                    }
                    else
                        new_pairs[i].first = -1;
                }
                else
                    new_pairs[i].first = -1;
            }

            /* Compute sigma_j */
            for (int i : mon.m()) {
                MMilnor m = MMilnor::FromIndex(i);
                AdamsDeg st{s, d_mon + m.deg()};
                if (st.t <= deg_trunc_ && st.t >= t_min_buffer)
                    ;
                // buffer_singles_[st].push_back(CriPairMRes::Single(m, (int)k));
            }

            /* Remove some critical pairs to form Groebner basis and discover redundent pairs */
            for (size_t j = 1; j < new_pairs.size(); ++j) {
                if (new_pairs[j].first != -1) {
                    for (size_t i = 0; i < j; ++i) {
                        if (new_pairs[i].first != -1) {
                            if (divisibleLF(new_pairs[i].second.m2, new_pairs[j].second.m2)) {
                                new_pairs[j].first = -1;
                                break;
                            }
                            else if (divisibleLF(new_pairs[j].second.m2, new_pairs[i].second.m2))
                                new_pairs[i].first = -1;
                            else if (!gcdLF(new_pairs[i].second.m1, new_pairs[j].second.m1)) {
                                int dij = lcmLF(leads[s][i].m(), leads[s][j].m()).deg() + d_mon_base;
                                if (dij <= deg_trunc_ && dij >= t_min_buffer)
                                    buffer_redundent_pairs_[AdamsDeg{s, dij}].insert(ut::BindPair((uint32_t)i, (uint32_t)j));
                            }
                        }
                    }
                }
            }

            /* Add critical pairs to buffer and `pairs_` */
            pairs_[s].resize(k + 1);
            for (size_t i = 0; i < new_pairs.size(); ++i) {
                if (new_pairs[i].first != -1) {
                    if (new_pairs[i].first >= t_min_buffer) {
                        buffer_min_pairs_[AdamsDeg{s, new_pairs[i].first}].resize(k + 1);
                        buffer_min_pairs_[AdamsDeg{s, new_pairs[i].first}][k].push_back(new_pairs[i].second);
                    }
                    else {
                        pairs_[s][k].push_back(new_pairs[i].second);
                    }
                }
            }
        }
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
            int t = rel.x1.GetLead().m().deg() + basis_degrees[rel.x1.GetLead().v()];
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

std::ostream& operator<<(std::ostream& sout, const MMilnorE& x)
{
    for (size_t i = 0; i < MMILNOR_INDEX_NUM; ++i)
        if (x.data[i + 1] > 0)
            std::cout << "P_{" << std::to_string(MMILNOR_GEN_I[i]) << std::to_string(MMILNOR_GEN_J[i]) << '}' << "^" << std::to_string(x.data[i + 1]);
    return sout;
}

void AddRelsMRes(GroebnerMRes& gb, const Mod1d& rels, int deg)
{
    int deg_max = gb.deg_trunc();
    if (deg > deg_max)
        throw MyException(0xb2474e19U, "deg is bigger than the truncation degree.");

    /*DbSteenrod db("AdamsE2.db");
    db.create_generators("SteenrodMRes");
    db.create_relations("SteenrodMRes");*/

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, array> rels_graded;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            const auto& lead = rels[i].GetLead();
            int t = lead.m().deg();
            if (t <= deg)
                rels_graded[t].push_back((int)i);
        }
    }
    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    auto next_st = gb.gb_pairs().next_st();
    /* Check if the program was terminated exact when adding elements from `rels` */
    if (gb.basis_degs().size() > 1 && next_st.t == gb.basis_degs()[1].back() * 2 + 1)
        --next_st.t;
    for (int t = next_st.t; t <= deg && (t <= deg_max_rels || !gb.gb_pairs().empty_pairs_for_gb()); ++t) {  ////
        std::cout << "t=" << t << "               \n";
        for (int s = t - 2; s >= -1; --s) {
            std::cout << "  s=" << s + 2 << "                \r";

            gb.resize_data(s + 1);
            auto st = AdamsDeg{s, t};
            gb.MinimizePairs(st);

            CriPairMRes1d pairs_st = gb.pairs(st);
            DataMRes1d rels_tmp;// TODO: change type to (Mod, Mod)
            if (s == -1) {
                auto p_rels_d = rels_graded.find(t);
                if (p_rels_d != rels_graded.end()) {
                    rels_tmp.resize(p_rels_d->second.size());
                    auto& rels_d = p_rels_d->second;
                    ut::for_each_seq((int)rels_tmp.size(), [&](int i) { rels_tmp[i] = DataMRes(Mod(), rels[rels_d[i]], Mod(), rels[rels_d[i]]); });
                }
            }
            else if (!pairs_st.empty()) {
                Mod1d x2_LFs(pairs_st.size());
                 ut::for_each_seq((int)x2_LFs.size(), [&](int i) { x2_LFs[i] = gb.ReduceX2LF(pairs_st[i], s); });
                array ind;
                for (size_t i = 0; i < x2_LFs.size(); ++i) {
                    for (size_t j = 0; j < ind.size(); ++j)
                        if (std::binary_search(x2_LFs[i].data.begin(), x2_LFs[i].data.end(), x2_LFs[ind[j]].GetLead(), cmpLF))
                            x2_LFs[i].iaddLF(x2_LFs[ind[j]]);
                    if (x2_LFs[i])
                        ind.push_back((int)i);
                    else
                        bench::Counter(4);
                }

                rels_tmp.resize(ind.size());
                ut::for_each_seq((int)rels_tmp.size(), [&](int i) { rels_tmp[i] = gb.Reduce(pairs_st[ind[i]], s); });
            }
            if (rels_tmp.empty())
                continue;

            /* Triangulate these relations */
            DataMRes1d rels_st;
            Mod1d kernel_splus_tmp;
            for (size_t i = 0; i < rels_tmp.size(); ++i) {
                for (size_t j = 0; j < rels_st.size(); ++j) {
                    if (std::binary_search(rels_tmp[i].x1.data.begin(), rels_tmp[i].x1.data.end(), rels_st[j].x1.GetLead(), gb.GetCmpMMod(s))) {
                        rels_tmp[i].x1.iadd(rels_st[j].x1, gb.GetCmpMMod(s));
                        rels_tmp[i].x2.iadd(rels_st[j].x2, gb.GetCmpMMod(s + 1));
                    }
                }
                if (rels_tmp[i].x1) {
                    rels_st.push_back(std::move(rels_tmp[i]));
                    bench::Counter(0);
                }
                else if (rels_tmp[i].x2)
                    kernel_splus_tmp.push_back(std::move(rels_tmp[i].x2));
                else
                    bench::Counter(1);
            }

            Mod1d kernel_sp1;
            for (size_t i = 0; i < kernel_splus_tmp.size(); ++i) {
                for (size_t j = 0; j < kernel_sp1.size(); ++j) {
                    if (std::binary_search(kernel_splus_tmp[i].data.begin(), kernel_splus_tmp[i].data.end(), kernel_sp1[j].GetLead(), gb.GetCmpMMod(s + 1)))
                        kernel_splus_tmp[i].iadd(kernel_sp1[j], gb.GetCmpMMod(s + 1));
                }
                if (kernel_splus_tmp[i]) {
                    kernel_sp1.push_back(std::move(kernel_splus_tmp[i]));
                    bench::Counter(2);
                }
                else
                    bench::Counter(3);
            }
            
            /* Add these relations */
            for (size_t i = 0; i < rels_st.size(); ++i)
                gb.push_back(rels_st[i], s, t);
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

            if (size_t size_k = kernel_sp1.size())
                std::cout << "  s=" << s + 2 << " dim=" << size_k << '\n';
        }
    }
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