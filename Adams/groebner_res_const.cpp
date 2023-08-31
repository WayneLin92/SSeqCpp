#include "groebner_res_const.h"
#include "algebras/database.h"

/********************************************************
 *                    class AdamsResConst
 ********************************************************/

/** Return i such that mon divides leads[i].
 *
 * Return -1 if not found.
 * @param indices indices[v_raw] stores the indices of elements of leads with the given v_raw.
 */
inline int IndexOfDivisibleLeading(const MMod1d& leads, const std::unordered_map<uint64_t, int1d>& indices, MMod mon)
{
    auto key = mon.v_raw();
    auto p = indices.find(key);
    if (p != indices.end())
        for (int k : p->second)
            if (divisibleLF(leads[k], mon))
                return k;
    return -1;
}

AdamsResConst::AdamsResConst(DataMResConst2d data, int2d basis_degrees) : gb_(std::move(data)), basis_degrees_(std::move(basis_degrees))
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
}

Mod AdamsResConst::DiffInv(Mod x, size_t s) const
{
    Milnor tmp_a;
    Mod result, tmp_x1, tmp_x2;
    tmp_a.data.reserve(50);
    tmp_x1.data.reserve(100);
    tmp_x2.data.reserve(100);

    size_t index;
    index = 0;
    while (index < x.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[s], indices_[s], x.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(x.data[index], gb_[s][gb_index].x1.data[0]);
            x.iaddmulP(m, gb_[s][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
            result.iaddmulP(m, gb_[s][gb_index].x2, tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }
    if (x) {
        fmt::print("Something is wrong: d_inv(x) not well defined. s={}\n", s);
        std::exit(-2);
    }
    size_t sp1 = s + 1;
    index = 0;
    while (index < result.data.size()) {
        int gb_index = IndexOfDivisibleLeading(leads_[sp1], indices_[sp1], result.data[index]);
        if (gb_index != -1) {
            MMilnor m = divLF(result.data[index], gb_[sp1][gb_index].x1.data[0]);
            result.iaddmulP(m, gb_[sp1][gb_index].x1, tmp_a, tmp_x1, tmp_x2);
        }
        else
            ++index;
    }

    return result;
}

struct IndexMMod
{
    MMod m;
    unsigned i, index;
    bool operator<(IndexMMod rhs) const
    {
        return rhs.m < m;
    }
};
using IndexMMod1d = std::vector<IndexMMod>;

void AdamsResConst::DiffInvBatch(Mod1d xs, Mod1d& result, size_t s) const
{
    Mod tmp_x, prod_x1, prod_x2;
    Milnor tmp_a;
    tmp_x.data.reserve(64);
    prod_x1.data.reserve(64);
    prod_x2.data.reserve(64);
    tmp_a.data.reserve(64);

    IndexMMod1d heap;
    heap.reserve(64);
    for (size_t i = 0; i < xs.size(); ++i)
        if (xs[i])
            heap.push_back(IndexMMod{xs[i].data[0], (unsigned)i, 0});
    std::make_heap(heap.begin(), heap.end());

    while (!heap.empty()) {
        MMod term = heap.front().m;
        int gb_index = s < leads_.size() ? IndexOfDivisibleLeading(leads_[s], indices_[s], term) : -1;
        if (gb_index != -1) {
            MMilnor m = divLF(term, gb_[s][gb_index].x1.data[0]);
            mulP(m, gb_[s][gb_index].x1, prod_x1, tmp_a);
            mulP(m, gb_[s][gb_index].x2, prod_x2, tmp_a);

            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());

                xs[i].iaddP(prod_x1, tmp_x);
                result[i].iaddP(prod_x2, tmp_x);

                if (index < xs[i].data.size()) {
                    heap.back() = IndexMMod{xs[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
        else {
            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());
                if (++index < xs[i].data.size()) {
                    heap.back() = IndexMMod{xs[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
    }

    for (size_t i = 0; i < xs.size(); ++i)
        if (xs[i]) {
            fmt::print("Something is wrong: d_inv(x) not well defined. s={}\n", s);
            std::exit(-2);
        }

    size_t sp1 = s + 1;
    heap.clear();
    for (size_t i = 0; i < result.size(); ++i)
        if (result[i])
            heap.push_back(IndexMMod{result[i].GetLead(), (unsigned)i, 0});
    std::make_heap(heap.begin(), heap.end());

    while (!heap.empty()) {
        MMod term = heap.front().m;
        int gb_index = sp1 < leads_.size() ? IndexOfDivisibleLeading(leads_[sp1], indices_[sp1], term) : -1;
        if (gb_index != -1) {
            MMilnor m = divLF(term, gb_[sp1][gb_index].x1.data[0]);
            mulP(m, gb_[sp1][gb_index].x1, prod_x1, tmp_a);

            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());

                result[i].iaddP(prod_x1, tmp_x);

                if (index < result[i].data.size()) {
                    heap.back() = IndexMMod{result[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
        else {
            while (!heap.empty() && heap.front().m == term) {
                unsigned i = heap.front().i, index = heap.front().index;
                std::pop_heap(heap.begin(), heap.end());
                if (++index < result[i].data.size()) {
                    heap.back() = IndexMMod{result[i].data[index], i, index};
                    std::push_heap(heap.begin(), heap.end());
                }
                else
                    heap.pop_back();
            }
        }
    }
}

int2d DbAdamsResLoader::load_basis_degrees(const std::string& table_prefix, int t_trunc) const
{
    int2d result;
    Statement stmt(*this, "SELECT s, t FROM " + table_prefix + "_generators WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
    while (stmt.step() == MYSQLITE_ROW) {
        int s = stmt.column_int(0), t = stmt.column_int(1);
        ut::get(result, s).push_back(t);
    }
    return result;
}

/* vid_num[s][stem] is the number of generators in (<=stem, s) */
void DbAdamsResLoader::load_generators(const std::string& table_prefix, std::vector<std::pair<int, AdamsDegV2>>& id_st, int2d& vid_num, std::map<AdamsDegV2, Mod1d>& diffs, int t_trunc, int stem_trunc) const
{
    Statement stmt(*this, fmt::format("SELECT id, s, t, diff FROM {}_generators WHERE t<={} AND t-s<={} ORDER BY id;", table_prefix, t_trunc, stem_trunc));
    AdamsDegV2 d_prev = AdamsDegV2(-1, -1);
    ut::map_seq2d<int, 0> gen_num; /* gen_num[(s,stem)] is the number of generators in (stem, s) */
    while (stmt.step() == MYSQLITE_ROW) {
        int id = stmt.column_int(0);
        AdamsDegV2 d = AdamsDegV2(stmt.column_int(1), stmt.column_int(2));
        Mod diff;
        diff.data = stmt.column_blob_tpl<MMod>(3);
        diffs[d].push_back(std::move(diff));
        ++gen_num(d.s, d.stem());
        if (d.s != d_prev.s || d.t != d_prev.t) {
            id_st.push_back(std::make_pair(id, d));
            d_prev = d;
        }
    }
    std::sort(id_st.begin(), id_st.end(), [](std::pair<int, AdamsDegV2> p1, std::pair<int, AdamsDegV2> p2) { return p1.second < p2.second; });

    vid_num.resize(size_t(t_trunc + 1));
    for (size_t s = 0; s <= t_trunc; ++s) {
        vid_num[s].resize(t_trunc - s + 1);
        for (size_t n = 0; n <= t_trunc - s; ++n) {
            if (n > 0)
                vid_num[s][n] += vid_num[s][n - 1];
            vid_num[s][n] += gen_num(s, n);
        }
    }
}

DataMResConst2d DbAdamsResLoader::load_data(const std::string& table_prefix, int t_trunc) const
{
    DataMResConst2d data;
    Statement stmt(*this, "SELECT x1, x2, s FROM " + table_prefix + "_relations WHERE t<=" + std::to_string(t_trunc) + " ORDER BY id;");
    while (stmt.step() == MYSQLITE_ROW) {
        DataMResConst x;
        x.x1.data = stmt.column_blob_tpl<MMod>(0);
        x.x2.data = stmt.column_blob_tpl<MMod>(1);

        size_t s = (size_t)stmt.column_int(2);
        ut::get(data, s).push_back(std::move(x));
    }
    return data;
}

DataMResConst1d DbAdamsResLoader::load_data_s(const std::string& table_prefix, int s, int t_trunc) const
{
    DataMResConst1d data;
    Statement stmt(*this, fmt::format("SELECT x1, x2 FROM {}_relations WHERE s={} AND t<={} ORDER BY id;", table_prefix, s, t_trunc));
    while (stmt.step() == MYSQLITE_ROW) {
        DataMResConst x;
        x.x1.data = stmt.column_blob_tpl<MMod>(0);
        x.x2.data = stmt.column_blob_tpl<MMod>(1);
        data.push_back(std::move(x));
    }
    return data;
}