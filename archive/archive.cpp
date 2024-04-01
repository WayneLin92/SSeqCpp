

/* Find the indices of leading lm from [begin, end) such that when ignoring the exponents
 * lm[-j-1] divides mon[i] and lm[:-j-1] divides mon[:i].
 */
int IndexOfDivisibleLeading_(const Mon& mon, size_t i, size_t j, IterIndex begin, IterIndex end) const
{
    IterIndex begin1 = std::lower_bound(begin, end, mon[i].gen, [j](const EleIndex& p, int gen) { return p.first[j] < gen; });
    IterIndex end1 = std::lower_bound(begin1, end, mon[i].gen + 1, [j](const EleIndex& p, int gen) { return p.first[j] < gen; });
    if (begin1 != end1) {
        if (begin1->first.size() == j + 1) {
            for (int i : begin1->second) {
                if (divisible(data[i].GetLead(), mon))
                    return i;
            }
            ++begin1;
        }
        for (size_t l = 0; l < i; ++l) {
            int ind = IndexOfDivisibleLeading_(mon, l, j + 1, begin1, end1);
            if (ind != -1)
                return ind;
        }
    }
    return -1;
}
/* Milnor's multiplication formula.
 * `result.data` is unordered and may contain duplicates.
 * R: row, S: column
 */
void MulMilnor(MMay lhs, MMay rhs, May& result)  // TODO: improve
{
    constexpr size_t N = MILNOR_BUFFER_SIZE;
    MMilnor R = lhs.ToMMilnor();
    MMilnor S = rhs.ToMMilnor();

    std::array<int, (N + 1) * (N + 1)> X = {};
    std::array<int, (N + 1) * (N + 1)> XR = {};
    std::array<int, (N + 1) * (N + 1)> XS = {};
    std::array<int, (N + 1) * (N + 1)> XT = {};
    for (size_t i = 1; i <= N; ++i)
        XR[i * (N + 1) + (N - i)] = R.data[i - 1];
    for (size_t i = 1; i <= N; ++i)
        XS[(N - i) * (N + 1) + i] = S.data[i - 1];

    size_t i = N, j = 0;
    bool decrease = false;
    while (true) {
        bool move_right = false;
        if (j) {
            size_t index = size_t(i * (N + 1) + j);
            size_t index_up = size_t((i - 1) * (N + 1) + j);
            size_t index_up_right = size_t((i - 1) * (N + 1) + j + 1);
            size_t index_left = size_t(i * (N + 1) + j - 1);
            if (i == 1) {
                if (decrease) {
                    X[index] = (X[index] - 1) & ~(XT[index] | X[index_up_right]);
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index] >> j, XS[index]), XT[index] | X[index_up_right]);
                X[index_up] = XS[index] - X[index];
                if (X[index_up] & XT[index_left]) {
                    if (X[index])
                        decrease = true;
                    else
                        move_right = true;
                }
                else {
                    XR[size_t(i * (N + 1) + j - 1)] = XR[index] - (X[index] << j);
                    XT[index_up_right] = XT[index] | X[index] | X[index_up_right];
                    --j;
                }
            }
            else {
                if (decrease) {
                    X[index] = (X[index] - 1) & ~XT[index];
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index] >> j, XS[index]), XT[index]);
                XR[size_t(i * (N + 1) + j - 1)] = XR[index] - (X[index] << j);
                XS[index_up] = XS[index] - X[index];
                XT[index_up_right] = XT[index] | X[index];
                --j;
            }
        }
        else {
            if (i == 1) {
                if (!(XR[N + 1] & X[1])) { /* Add to result. */
                    XT[1] = XR[N + 1] | X[1];
                    uint64_t data = 0;
                    uint64_t w = 0;
                    for (int d = 1; d <= N; ++d) {
                        for (int n = XT[d], i = 0; n; n >>= 1, ++i) {
                            if (n & 1) {
                                int j = i + d;
                                data |= MMay::rawP(i, j);
                                w += 2 * uint64_t(d) - 1;
                            }
                        }
                    }
                    result.data.push_back(MMay(data + (w << MMAY_INDEX_NUM)));
                }
                move_right = true;
            }
            else {
                size_t index = size_t(i * (N + 1));
                size_t index_up_right = size_t((i - 1) * (N + 1) + 1);
                XT[index_up_right] = XR[index];
                j = N - (--i);
            }
        }
        if (move_right) {
            size_t index = i * (N + 1) + j + 1;
            while (index <= N * (N + 1) && X[index] == 0)
                ++index;
            if (index > N * (N + 1))
                break;
            i = index / (N + 1);
            j = index % (N + 1);
            decrease = true;
        }
    }
}

/* R: column, S: row */
void MulMilnor(int R[MILNOR_BUFFER_SIZE], int S[MILNOR_BUFFER_SIZE], May& result)
{
    constexpr size_t N = MILNOR_BUFFER_SIZE;

    std::array<int, (N + 1) * (N + 1)> X = {};  // Improve
    std::array<int, (N + 1) * (N + 1)> XR;
    std::array<int, (N + 1) * (N + 1)> XS;
    std::array<int, (N + 1) * (N + 1)> XT;
    for (size_t i = 1; i <= N; ++i)
        XR[(N - i) * (N + 1) + i] = R[i - 1];
    for (size_t i = 1; i <= N; ++i)
        XS[i * (N + 1) + (N - i)] = S[i - 1];
    for (size_t i = 1; i <= N; ++i)
        XT[i * (N + 1)] = 0;

    size_t i = N, j = 0;
    bool decrease = false;
    while (true) {
        bool move_right = false;
        if (j) {
            size_t index = size_t(i * (N + 1) + j);
            size_t index_up = size_t((i - 1) * (N + 1) + j);
            size_t index_up_right = size_t((i - 1) * (N + 1) + j + 1);
            size_t index_left = size_t(i * (N + 1) + j - 1);
            if (i == 1) {
                if (decrease) {
                    X[index] = (X[index] - 1) & ~(XT[index] | X[index_up_right]);
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index] >> i, XS[index]), XT[index] | X[index_up_right]);
                X[index_up] = XR[index] - (X[index] << i);
                if (X[index_up] & XT[index_left]) {
                    if (X[index])
                        decrease = true;
                    else
                        move_right = true;
                }
                else {
                    XS[size_t(i * (N + 1) + j - 1)] = XS[index] - X[index];
                    XT[index_up_right] = XT[index] | X[index] | X[index_up_right];
                    --j;
                }
            }
            else {
                if (decrease) {
                    X[index] = (X[index] - 1) & ~XT[index];
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index] >> i, XS[index]), XT[index]);
                XS[size_t(i * (N + 1) + j - 1)] = XS[index] - X[index];
                XR[index_up] = XR[index] - (X[index] << i);
                XT[index_up_right] = XT[index] | X[index];
                --j;
            }
        }
        else {
            if (i == 1) {
                if (!(XS[N + 1] & X[1])) { /* Add to result. */
                    XT[1] = XS[N + 1] | X[1];
                    uint64_t data = 0;
                    uint64_t w = 0;
                    for (int d = 1; d <= N; ++d) {
                        for (int n = XT[d], i = 0; n; n >>= 1, ++i) {
                            if (n & 1) {
                                int j = i + d;
                                data |= MMay::rawP(i, j);
                                w += 2 * uint64_t(d) - 1;
                            }
                        }
                    }
                    result.data.push_back(MMay(data + (w << MMAY_INDEX_NUM)));
                }
                move_right = true;
            }
            else {
                size_t index = size_t(i * (N + 1));
                size_t index_up_right = size_t((i - 1) * (N + 1) + 1);
                XT[index_up_right] = XS[index];
                j = N - (--i);
            }
        }
        if (move_right) {
            size_t index = i * (N + 1) + j + 1;
            while (index <= N * (N + 1) && X[index] == 0)
                ++index;
            if (index > N * (N + 1))
                break;
            i = index / (N + 1);
            j = index % (N + 1);
            decrease = true;
        }
    }
}

/* Reverse the bits */
inline uint64_t Reverse(uint64_t b)
{
    b = (b & 0xFFFFFFFF00000000) >> 32 | (b & 0x00000000FFFFFFFF) << 32;
    b = (b & 0xFFFF0000FFFF0000) >> 16 | (b & 0x0000FFFF0000FFFF) << 16;
    b = (b & 0xFF00FF00FF00FF00) >> 8 | (b & 0x00FF00FF00FF00FF) << 8;
    b = (b & 0xF0F0F0F0F0F0F0F0) >> 4 | (b & 0x0F0F0F0F0F0F0F0F) << 4;
    b = (b & 0xCCCCCCCCCCCCCCCC) >> 2 | (b & 0x3333333333333333) << 2;
    b = (b & 0xAAAAAAAAAAAAAAAA) >> 1 | (b & 0x5555555555555555) << 1;
    return b;
}

void generate_basis(const std::string& filename, int t_trunc)
{
    using namespace steenrod;
    std::cout << "Constructing basis..." << std::endl;

    std::map<int, steenrod::MMilnor1d> basis;
    basis[0].push_back(MMilnor(0));

    /* Add new basis */
    for (int t = 0; t <= t_trunc; ++t) {
        std::cout << "t=" << t << "          \r";
        for (size_t gen_id = steenrod::MMILNOR_INDEX_NUM; gen_id-- > 0;) {
            MMilnor m = MMilnor::FromIndex(gen_id);
            int t1 = t - MMILNOR_GEN_DEG[gen_id];
            if (t1 >= 0) {
                for (MMilnor m1 : basis[t1]) {
                    if (!m1 || gen_id < *m1.begin()) {
                        MMilnor prod = mulLF(m, m1);
                        basis[t].push_back(prod);
                    }
                }
            }
        }
        std::sort(basis[t].begin(), basis[t].end());
    }

    DbSteenrod db(filename);
    db.execute_cmd("CREATE TABLE IF NOT EXISTS Steenrod_basis (id INTEGER PRIMARY KEY, base BIGINT, t SMALLINT, w SMALLINT);");

    db.begin_transaction();
    myio::Statement stmt(db, "INSERT INTO Steenrod_basis (base, t, w) VALUES (?1, ?2, ?3);");
    for (int t = 0; t <= t_trunc; ++t) {
        for (MMilnor m : basis[t]) {
            stmt.bind_int64(1, (int64_t)m.data());
            stmt.bind_int(2, t);
            stmt.bind_int(3, m.weight());
            stmt.step_and_reset();
        }
    }
    db.end_transaction();
}

Mod operator*(MMilnor m, const Mod& x)
{
    /*Mod result;
    for (MMod m2 : x.data)
        MulMilnor(m, m2, result);
    SortMod2(result.data);
    return result;*/

    Mod result;
    Milnor prod;
    auto p = x.data.begin();
    while (p != x.data.end()) {
        auto v_raw = p->v_raw();
        auto p1 = std::upper_bound(p, x.data.end(), v_raw, [](uint64_t v_raw, MMod x) { return v_raw < x.v_raw(); });
        for (; p != p1; ++p)
            MulMilnor(m, p->m(), prod);
        SortMod2(prod.data);
        for (MMilnor T : prod.data)
            result.data.push_back(MMod::FromRaw(T.data(), v_raw));
        prod.data.clear();
    }
    return result;
}

template <uint32_t d>
constexpr std::array<uint64_t, d >= 5 ? (1 << (9 - d)) : (1 << 4)> MMilnorXiLow()
{
    constexpr uint32_t N = d >= 5 ? (1 << (9 - d)) : (1 << 4);
    std::array<uint64_t, N> result = {};
    for (uint32_t m = 0; m < N; ++m) {
        for (uint32_t n = m, i = 0; n; n >>= 1, ++i) {
            if (n & 1) {
                uint32_t j = i + d;
                uint32_t index = j * (j + 1) / 2 - i - 1;
                result[m] |= MMILNOR_ONE >> index;
            }
        }
    }
    return result;
}
template <uint32_t d>
constexpr std::array<uint64_t, 1 << (5 - d)> MMilnorXiHigh()
{
    constexpr uint32_t N = 1 << (5 - d);
    std::array<uint64_t, N> result = {};
    for (uint32_t m = 0; m < N; ++m) {
        for (uint32_t n = (m << 4), i = 0; n; n >>= 1, ++i) {
            if (n & 1) {
                uint32_t j = i + d;
                uint32_t index = j * (j + 1) / 2 - i - 1;
                result[m] |= MMILNOR_ONE >> index;
            }
        }
    }
    return result;
}

inline constexpr auto MMILNOR_Xi0_h = detail::MMilnorXiHigh<1>();
inline constexpr auto MMILNOR_Xi0 = detail::MMilnorXiLow<1>();
inline constexpr auto MMILNOR_Xi1_h = detail::MMilnorXiHigh<2>();
inline constexpr auto MMILNOR_Xi1 = detail::MMilnorXiLow<2>();
inline constexpr auto MMILNOR_Xi2_h = detail::MMilnorXiHigh<3>();
inline constexpr auto MMILNOR_Xi2 = detail::MMilnorXiLow<3>();
inline constexpr auto MMILNOR_Xi3_h = detail::MMilnorXiHigh<4>();
inline constexpr auto MMILNOR_Xi3 = detail::MMilnorXiLow<4>();
inline constexpr auto MMILNOR_Xi4 = detail::MMilnorXiLow<5>();
inline constexpr auto MMILNOR_Xi5 = detail::MMilnorXiLow<6>();
inline constexpr auto MMILNOR_Xi6 = detail::MMilnorXiLow<7>();

static MMilnor Xi(const uint32_t* xi)
{
    const uint64_t w_raw = uint64_t(1 * ut::popcount(xi[0]) + 3 * ut::popcount(xi[1]) + 5 * ut::popcount(xi[2]) + 7 * ut::popcount(xi[3]) + 9 * ut::popcount(xi[4]) + 11 * ut::popcount(xi[5]) + 13 * ut::popcount(xi[6])) << MMILNOR_E_BITS;
    const uint64_t e = MMILNOR_Xi0_h[xi[0] >> 4] | MMILNOR_Xi0[xi[0] & 0xf] | MMILNOR_Xi1_h[xi[1] >> 4] | MMILNOR_Xi1[xi[1] & 0xf] | MMILNOR_Xi2_h[xi[2] >> 4] | MMILNOR_Xi2[xi[2] & 0xf] | MMILNOR_Xi3_h[xi[3] >> 4] | MMILNOR_Xi3[xi[3] & 0xf]
                       | MMILNOR_Xi4[xi[4]] | MMILNOR_Xi5[xi[5]] | MMILNOR_Xi6[xi[6]];
    return MMilnor(w_raw + e);
}

void subs_batch(const GenMRes1d& gens, size_t i_start, const Mod1d& map, Mod1d::iterator result)
{
    Milnor tmp_a;
    tmp_a.data.reserve(128);
    Mod prod, tmp_x;
    prod.data.reserve(128);
    tmp_x.data.reserve(128);

    IndexMMod1d heap;
    for (size_t i = i_start; i < gens.size(); ++i)
        if (gens[i].diff)
            heap.push_back(IndexMMod{gens[i].diff.data[0], (unsigned)i, 0});
    std::make_heap(heap.begin(), heap.end());

    while (!heap.empty()) {
        MMod term = heap.front().m;
        prod.data.clear();
        mulP(term.m(), map[term.v()], tmp_a, prod);
        int b = 0;
        while (!heap.empty() && heap.front().m == term) {
            if (b++)
                bench::Counter(2);
            else
                bench::Counter(3);

            unsigned i = heap.front().i, index = heap.front().index;
            std::pop_heap(heap.begin(), heap.end());

            result[i].iaddP(prod, tmp_x);

            if (++index < gens[i].diff.data.size()) {
                heap.back() = IndexMMod{gens[i].diff.data[index], i, index};
                std::push_heap(heap.begin(), heap.end());
            }
            else
                heap.pop_back();
        }
    }
}

/**********************************************************
 * Algorithms that use Groebner basis
 **********************************************************/

inline constexpr int GEN_IDEAL = 0x4000000;
inline constexpr int GEN_SUB = 0x2000000;

namespace detail {
/* For monomials in the homological groebner basis get the degree of the first part consisting of original generators and set `p` to be
 * the end location of the first part.
 */
template <typename TypeMonIter>
int GetPart1Deg(const TypeMonIter mon_begin, const TypeMonIter mon_end, const int1d& gen_degs, TypeMonIter& p)
{
    int result = 0;
    for (p = mon_begin; p != mon_end && !(p->g() & GEN_SUB); ++p)
        result += gen_degs[p->g()] * p->e();
    return result;
};
}  // namespace detail

/**
 * Revlex on extra generators and FnCmp on the rest
 */
template <typename FnCmp>
struct CmpIdeal
{
    using submo = FnCmp;
    static constexpr std::string_view name = "Ideal";
    template <typename Type1, typename Type2>
    static bool cmp(const Type1& m1, const Type2& m2)
    {
        auto mid1 = std::lower_bound(m1.begin(), m1.end(), GEN_IDEAL, [](const GenPow& p, int g) { return p.gen < g; });
        auto mid2 = std::lower_bound(m2.begin(), m2.end(), GEN_IDEAL, [](const GenPow& p, int g) { return p.gen < g; });
        if (CmpRevlex::cmp_ranges(mid1, m1.end(), mid2, m2.end()))
            return true;
        else if (std::equal(mid1, m1.end(), mid2, m2.end())) {
            if (FnCmp::template cmp_ranges(m1.begin(), mid1, m2.begin(), mid2))
                return true;
        }
        return false;
    }
};

/*
 * The monomial ordering that computes the homology.
 */
template <typename FnCmp>
struct CmpHomology
{
    using submo = FnCmp;
    static constexpr std::string_view name = "Homology";
    inline static int1d gen_degs;
    template <typename Type1, typename Type2>
    static bool cmp(const Type1& m1, const Type2& m2)
    {
        typename Type1::const_iterator m1_mid1;
        typename Type2::const_iterator m2_mid1;
        int d1 = detail::template GetPart1Deg(m1.begin(), m1.end(), gen_degs, m1_mid1);
        int d2 = detail::template GetPart1Deg(m2.begin(), m2.end(), gen_degs, m2_mid1);
        if (d1 > d2)
            return true;
        if (d1 < d2)
            return false;

        if (FnCmp::template cmp_ranges(m1.begin(), m1_mid1, m2.begin(), m2_mid1))
            return true;
        if (FnCmp::template cmp_ranges(m2.begin(), m2_mid1, m1.begin(), m1_mid1))
            return false;

        auto m1_mid2 = std::lower_bound(m1_mid1, m1.end(), GEN_IDEAL + GEN_SUB, [](const GenPow& p, int g) { return p.gen < g; });
        auto m2_mid2 = std::lower_bound(m2_mid1, m2.end(), GEN_IDEAL + GEN_SUB, [](const GenPow& p, int g) { return p.gen < g; });
        if (CmpRevlex::cmp_ranges(m1_mid2, m1.end(), m2_mid2, m2.end()))
            return true;
        if (CmpRevlex::cmp_ranges(m2_mid2, m2.end(), m1_mid2, m1.end()))
            return false;

        if (CmpRevlex::cmp_ranges(m1_mid1, m1_mid2, m2_mid1, m2_mid2))
            return true;
        return false;
    }
};

namespace detail {
template <typename FnCmp>
void AddRelsIdeal(Groebner& gb, const std::vector<Poly>& rels, int deg, const int1d& gen_degs, const int1d& gen_degs_y)
{
    TplAddRels(
        gb, rels, deg, [](const Mon& m1, const Mon& m2) { return !(m1.back().gen & GEN_IDEAL) && !(m2.back().gen & GEN_IDEAL); }, [&gen_degs, &gen_degs_y](int i) { return ((i & GEN_IDEAL) ? gen_degs_y[i - GEN_IDEAL] : gen_degs[i]); });
}
template <typename FnCmp>
void AddRelsModule(Groebner& gb, const std::vector<Poly>& rels, int deg, const int1d& gen_degs, const int1d& gen_degs_y)
{
    TplAddRels(
        gb, rels, deg, [](const Mon& m1, const Mon& m2) { return !((m1.back().gen & GEN_IDEAL) && (m2.back().gen & GEN_IDEAL) && (m1.back().gen != m2.back().gen)); },
        [&gen_degs, &gen_degs_y](int i) { return ((i & GEN_IDEAL) ? gen_degs_y[i - GEN_IDEAL] : gen_degs[i]); });
}
template <typename FnCmp>
void AddRelsHomology(Groebner& gb, const std::vector<Poly>& rels, int deg, const MayDeg1d& gen_degs, const MayDeg1d& gen_degs_x, const MayDeg1d& gen_degs_b)
{
    TplAddRels(
        gb, rels, deg, [](const Mon& m1, const Mon& m2) { return true; },
        [&gen_degs, &gen_degs_x, &gen_degs_b](int i) {
            if (i & GEN_IDEAL)
                return gen_degs_b[i - GEN_IDEAL - GEN_SUB].t;
            else if (i & GEN_SUB)
                return gen_degs_x[i - GEN_SUB].t;
            else
                return gen_degs[i].t;
        });
}
}  // namespace detail

/**
 * Compute the minimal generating set of `vectors` inplace.
 *
 * Each element of vectors is considered as an element of $R^n$ where $R$ is the algebra
 * determined by the Groebner basis `gb`.
 *
 * The vectors should all be nontrivial.
 *
 * The degree of the basis of $R^n$ is determined by `basis_degs`.
 */
template <typename FnCmp>
std::vector<std::vector<Poly>>& Indecomposables(const Groebner& gb, std::vector<std::vector<Poly>>& vectors, const int1d& gen_degs, const int1d& basis_degs)
{
    using Poly = Poly;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner;
    using PolyI = Poly<CmpIdeal<FnCmp>>;
    using PolyI1d = std::vector<PolyI>;
    using GbI = alg::Groebner<CmpIdeal<FnCmp>>;

    if (vectors.empty())
        return vectors;

    /* Convert each vector v into a relation \\sum v_iy_i */
    PolyI1d rels;
    int1d degs;
    for (auto& v : vectors) {
        PolyI rel;
        for (size_t i = 0; i < basis_degs.size(); ++i)
            if (v[i])
                rel += PolyI::Sort((v[i] * GenPow::Mon(GEN_IDEAL + (int)i)).data);
        for (size_t i = 0; i < basis_degs.size(); ++i) {
            if (v[i]) {
                degs.push_back(v[i].GetDeg(gen_degs) + basis_degs[i]);
                break;
            }
        }
        if (rel)
            rels.push_back(std::move(rel));
        else
            v.clear();
    }
    ut::RemoveEmptyElements(vectors);
    if (!vectors.empty()) {
        int1d indices = ut::int_range((int)vectors.size());
        std::sort(indices.begin(), indices.end(), [&degs](int i, int j) { return degs[i] < degs[j]; });

        /* Add relations ordered by degree to gb1 */
        int deg_max = degs[indices.back()];
        GbI gbI = gb.template ExtendMO<CmpIdeal<FnCmp>>();
        gbI.set_deg_trunc(deg_max);
        gbI.set_mode(CRI_ON);
        for (int i : indices) {
            detail::AddRelsModule(gbI, {}, degs[i], gen_degs, basis_degs);
            PolyI rel = gbI.Reduce(rels[i]);
            if (rel)
                detail::AddRelsModule(gbI, {std::move(rel)}, degs[i], gen_degs, basis_degs);
            else
                vectors[i].clear();
        }

        /* Keep only the indecomposables in `vectors` */
        ut::RemoveEmptyElements(vectors);
    }
    return vectors;
}

/**
 * Compute the generating set of linear relations among `polys`.
 *
 * The result is truncated by `deg<=deg_max`.
 */
template <typename FnCmp>
std::vector<std::vector<Poly>> AnnSeq(const Groebner& gb, const std::vector<Poly>& polys, const int1d& gen_degs, int deg_max)
{
    using Poly = Poly;
    using Poly1d = std::vector<Poly>;
    using Poly2d = std::vector<Poly1d>;
    using Gb = alg::Groebner;
    using PolyI = Poly<CmpIdeal<FnCmp>>;
    using PolyI1d = std::vector<PolyI>;
    using GbI = alg::Groebner<CmpIdeal<FnCmp>>;

    Poly2d result;
    if (polys.empty())
        return result;
    PolyI1d rels;
    int1d gen_degs_y;
    int n = (int)polys.size();

    /* Add relations Xi=polys[i] to gb to obtain gb1 */
    for (int i = 0; i < n; ++i) {
        PolyI p{polys[i].data};
        gen_degs_y.push_back(p.GetDeg(gen_degs));
        p += PolyI::Gen(GEN_IDEAL + i);
        rels.push_back(std::move(p));
    }
    GbI gbI = gb.template ExtendMO<CmpIdeal<FnCmp>>();
    gbI.set_deg_trunc(deg_max);
    gbI.set_mode(CRI_ON);
    detail::AddRelsIdeal(gbI, rels, deg_max, gen_degs, gen_degs_y);

    /* Extract linear relations from gb1 */
    for (const PolyI& g : gbI.data()) {
        if (g.GetLead().back().gen & GEN_IDEAL) {
            Poly1d ann;
            ann.resize(n);
            for (const Mon& m : g.data) {
                auto mid = std::lower_bound(m.begin(), m.end(), GEN_IDEAL, [](const GenPow& p, int g) { return p.gen < g; });
                Mon m1(m.begin(), mid), m2(mid, m.end());  // TODO: improve
                ann[m2.front().gen - GEN_IDEAL] += gb.Reduce(SubsMGb({div(m2, GenPow::Mon(m2.front().gen))}, gb, [&polys](int i) { return polys[i - GEN_IDEAL]; }) * m1);
            }
            result.push_back(std::move(ann));
        }
    }

    /* Add commutators to linear relations */
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (gen_degs_y[i] + gen_degs_y[j] <= deg_max) {
                Poly1d result_i;
                result_i.resize(n);
                result_i[i] = polys[j];
                result_i[j] = polys[i];
                result.push_back(std::move(result_i));
            }
        }
    }

    Indecomposables(gb, result, gen_degs, gen_degs_y);
    return result;
}

/* Propergate `basis` to degree t=`t` */
template <typename FnCmp>
std::map<MayDeg, Mon1d> ExtendBasis(const Mon2d& leadings, const MayDeg1d& gen_degs, std::map<MayDeg, Mon1d>& basis, int t)
{
    std::map<MayDeg, Mon1d> basis_t;
    if (t == 0) {
        basis[MayDeg{0, 0, 0}].push_back({});
        return basis_t;
    }
    for (int gen_id = 0; gen_id < (int)gen_degs.size(); ++gen_id) {
        int t1 = t - gen_degs[gen_id].t;
        if (t1 >= 0) {
            auto basis_t1_begin = basis.lower_bound(MayDeg{0, t1, 0});
            auto basis_t1_end = basis.lower_bound(MayDeg{0, t1 + 1, 0});
            for (auto p_basis = basis_t1_begin; p_basis != basis_t1_end; ++p_basis) {
                for (auto p_m = p_basis->second.begin(); p_m != p_basis->second.end(); ++p_m) {
                    if (p_m->empty() || gen_id >= p_m->back().gen) {
                        Mon mon = mul(*p_m, GenPow::Mon(gen_id));
                        if (std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(), [&mon](const Mon& _m) { return divisible(_m, mon); }))
                            basis_t[p_basis->first + gen_degs[gen_id]].push_back(std::move(mon));
                    }
                }
            }
        }
    }

    for (auto& p : basis_t)
        std::sort(p.second.begin(), p.second.end(), FnCmp::template cmp<Mon, Mon>);
    return basis_t;
}

template <typename FnCmp>
void Homology(const Groebner& gb, const MayDeg1d& gen_degs, /*std::vector<std::string>& gen_names, */ const std::vector<Poly>& gen_diffs, const MayDeg& deg_diff, GroebnerRevlex& gb_h, MayDeg1d& gen_degs_h, std::vector<std::string>& gen_names_h,
              std::vector<Poly>& gen_repr_h, int t_trunc)
{
    using Poly = Poly;
    using Poly1d = std::vector<Poly>;
    using PolyH = Poly<CmpHomology<FnCmp>>;
    using PolyH1d = std::vector<PolyH>;
    using Gb = Groebner;
    using GbH = Groebner<CmpHomology<FnCmp>>;

    CmpHomology<FnCmp>::gen_degs.clear();
    for (const auto& deg : gen_degs)
        CmpHomology<FnCmp>::gen_degs.push_back(deg.t);

    if (gb.deg_trunc() < t_trunc)
        throw MyException(0x12ddc7ccU, "The truncation degree of the homology should be smaller");

    Gb gb_AoAZ = gb;
    gb_AoAZ.set_deg_trunc(t_trunc);
    gb_AoAZ.set_mode(CRI_ON);
    if (!(gb.mode() & CRI_ON))
        gb_AoAZ.CacheDegs(CmpHomology<FnCmp>::gen_degs);

    Mon2d leadings_AoAZ = gb_AoAZ.GetLeadings(gen_degs.size());
    size_t gb_AoAZ_old_size = gb_AoAZ.size();
    std::map<MayDeg, Mon1d> basis_AoAZ;
    std::vector<const Mon*> basis_AoAZ_flat;
    std::vector<Poly> diff_basis_AoAZ_flat;
    basis_AoAZ[MayDeg{0, 0, 0}].push_back({});

    Groebner<CmpHomology<FnCmp>> gb_H = gb.template ExtendMO<CmpHomology<FnCmp>>();
    gb_H.set_deg_trunc(t_trunc);
    gb_H.set_mode(CRI_ON);
    if (!(gb.mode() & CRI_ON))
        gb_H.CacheDegs(CmpHomology<FnCmp>::gen_degs);

    MayDeg1d gen_degs_b;
    Poly1d gen_reprs_b; /* dm_i = b_i */

    int deg_t_allX = 0;
    for (auto& deg : gen_degs)
        deg_t_allX += deg.t;
    // std::cout << "deg_t_allX=" << deg_t_allX << '\n';

    for (int t = 1; t <= t_trunc; ++t) {
        std::cout << "t=" << t << '\n';  //
        PolyH1d rels_h_t;
        Poly1d rels_AoAZ_t;
        if (t <= deg_t_allX) {
            /* Find m_i, b_i and x_i */
            for (size_t i = gb_AoAZ_old_size; i < gb_AoAZ.data().size(); ++i)
                leadings_AoAZ[gb_AoAZ[i].GetLead().back().gen].push_back(gb_AoAZ[i].GetLead());  // TODO: leads_
            gb_AoAZ_old_size = gb_AoAZ.size();
            std::map<MayDeg, Mon1d> basis_AoAZ_t = ExtendBasis<FnCmp>(leadings_AoAZ, gen_degs, basis_AoAZ, t);
            for (auto& [deg, b_deg] : basis_AoAZ_t) {
                int2d map_diff;
                for (const Mon& m : b_deg) {
                    Poly diff_m = gb.Reduce(GetDiff(m, gen_diffs));
                    map_diff.push_back(alg::hash1d(diff_m.data));
                }
                int2d image_diff, kernel_diff, g_diff;
                lina::SetLinearMap(map_diff, image_diff, kernel_diff, g_diff);

                /* Add x_i and the relations for x_i */
                for (const int1d& k : kernel_diff) {
                    gen_names_h.push_back("x_{" + std::to_string(gen_degs_h.size()) + "}");
                    gen_degs_h.push_back(deg);
                    gen_repr_h.push_back({Indices2Poly(k, b_deg)});
                    rels_AoAZ_t.push_back(gen_repr_h.back());
                    rels_h_t.push_back(PolyH({gen_repr_h.back().data}) + PolyH::Gen(GEN_SUB + (int)gen_degs_h.size() - 1));
                    b_deg[k.front()].clear();
                }
                ut::RemoveEmptyElements(b_deg);

                /* Add b_i, m_i and the relation dm_i = b_i */
                std::vector<Mon1d::iterator> toBeRemoved;
                MayDeg d_dm = deg + deg_diff;
                for (const Mon& m : b_deg) {
                    basis_AoAZ_flat.push_back(&m);
                    diff_basis_AoAZ_flat.push_back(gb.Reduce(GetDiff(m, gen_diffs)));
                    auto& dm = diff_basis_AoAZ_flat.back();
                    if (dm) {
                        rels_AoAZ_t.push_back(dm);
                        auto lead_dm = dm.GetLead();
                        auto ptr = std::lower_bound(basis_AoAZ_t[d_dm].begin(), basis_AoAZ_t[d_dm].end(), lead_dm, FnCmp::template cmp<Mon, Mon>);
                        if (ptr != basis_AoAZ_t[d_dm].end() && (*ptr == lead_dm))
                            toBeRemoved.push_back(ptr);
                    }
                    gen_degs_b.push_back(deg);
                    rels_h_t.push_back(PolyH({dm.data}) + PolyH::Gen(GEN_SUB + GEN_IDEAL + (int)gen_degs_b.size() - 1));
                }
                if (!toBeRemoved.empty()) {
                    for (auto p : toBeRemoved)
                        p->clear();
                    ut::RemoveEmptyElements(basis_AoAZ_t[d_dm]);
                }
            }
            basis_AoAZ.merge(basis_AoAZ_t);
            AddRels(gb_AoAZ, rels_AoAZ_t, t, CmpHomology<FnCmp>::gen_degs);
            rels_AoAZ_t.clear();
        }

        /* Find  y_i */
        size_t i_start_t = gb_H.size();
        detail::AddRelsHomology(gb_H, rels_h_t, t, gen_degs, gen_degs_h, gen_degs_b);
        rels_h_t.clear();

        for (size_t i = i_start_t; i < gb_H.size(); ++i) {
            const auto& lead = gb_H[i].GetLead();
            if ((lead.front().gen & GEN_SUB) && (lead.back().gen & GEN_IDEAL) && (lead.back().exp == 1) && (lead.size() == 1 || !(lead[lead.size() - 2].gen & GEN_IDEAL))) {
                PolyH y;
                for (const Mon& m : gb_H[i].data) {
                    auto mid = std::lower_bound(m.begin(), m.end(), GEN_IDEAL + GEN_SUB, [](const GenPow& p, int g) { return p.gen < g; });
                    Mon m1(m.begin(), mid), m2(mid, m.end());  // TODO: improve if bottlenet
                    y += gb_H.Reduce(PolyH::Mon_(mul(mul(m1, div(m2, GenPow::Mon(m2.front().gen))), *basis_AoAZ_flat[m2.front().gen - GEN_IDEAL - GEN_SUB])));
                }
                if (y && !(y.GetLead().front().gen & GEN_SUB)) {
                    Poly repr_y = SubsMGb(y.data, gb, [&gen_repr_h, &diff_basis_AoAZ_flat](int i) { return ((i & GEN_IDEAL) ? diff_basis_AoAZ_flat[i - GEN_IDEAL - GEN_SUB] : ((i & GEN_SUB) ? gen_repr_h[i - GEN_SUB] : Poly::Gen(i))); });
                    detail::AddRelsHomology(gb_H, {y + PolyH::Gen(GEN_SUB + (int)gen_degs_h.size())}, t, gen_degs, gen_degs_h, gen_degs_b);
                    gen_names_h.push_back("y_{" + std::to_string(gen_degs_h.size()) + "}");
                    gen_degs_h.push_back(repr_y.GetMayDeg(gen_degs));
                    gen_repr_h.push_back(std::move(repr_y));
                    rels_AoAZ_t.push_back(gen_repr_h.back());
                }
            }
        }
        AddRels(gb_AoAZ, rels_AoAZ_t, t, CmpHomology<FnCmp>::gen_degs);
        rels_AoAZ_t.clear();
    }
    /* Prepare relations for homology */
    PolyRevlex1d polys_h;
    for (auto& p : gb_H.data()) {
        if (p && (p.GetLead().front().gen & GEN_SUB)) {
            auto p_m = std::lower_bound(p.data.begin(), p.data.end(), 0, [](const Mon& m, int) { return !(m.back().gen & GEN_IDEAL); });
            if (p_m != p.data.begin()) {
                Mon1d p1(p.data.begin(), p_m);
                polys_h.push_back(alg::subs<alg::CmpRevlex>(p1, [](int i) { return PolyRevlex::Gen(i - GEN_SUB); }));
            }
        }
    }
    gb_h.InitFrom(polys_h);
}

#define TEMPLATE_EXAMPLES
#ifdef TEMPLATE_EXAMPLES
namespace template_examples {
using FnCmp = alg::CmpLex;
using FnPred = bool (*)(alg::Mon, alg::Mon);
using FnDeg = int (*)(int);

inline void TplAddRels_(Groebner& gb, const Poly1d<FnCmp>& rels, int deg, FnPred pred, FnDeg _gen_deg)
{
    return alg::TplAddRels(gb, rels, deg, pred, _gen_deg);
}

inline void Homology_(const Groebner& gb, const MayDeg1d& gen_degs, /*std::vector<std::string>& gen_names, */ const std::vector<Poly>& gen_diffs, const MayDeg& deg_diff, GroebnerRevlex& gb_h, MayDeg1d& gen_degs_h,
                      std::vector<std::string>& gen_names_h, std::vector<Poly>& gen_repr_h, int t_max)
{
    Homology(gb, gen_degs, gen_diffs, deg_diff, gb_h, gen_degs_h, gen_names_h, gen_repr_h, t_max);
}
}  // namespace template_examples
#endif

void benchmark_B9_Lex()
{
    std::cout << "benchmark_B9_Lex: \n";
    bench::Timer timer;

    using FnCmp = alg::CmpLex;
    using Poly = alg::Polynomial<FnCmp>;
    using Poly1d = std::vector<Poly>;
    using Gb = alg::Groebner<FnCmp>;
    constexpr auto GenExp = Poly::GenExp;

    int n_max = 9;
    alg::array gen_degs;
    for (int i = 0; i < n_max; ++i) {
        for (int j = i + 1; j <= n_max; ++j) {
            gen_degs.push_back((1 << j) - (1 << i));
        }
    }
    Poly1d rels;
    for (int d = 2; d <= n_max; d++) {
        for (int i = 0; i <= n_max - d; i++) {
            int j = i + d;
            Poly rel;
            for (int k = i + 1; k < j; k++) {
                int a = (1 << k) - (1 << i);
                int b = (1 << j) - (1 << k);
                auto p1 = std::find(gen_degs.begin(), gen_degs.end(), a);
                auto p2 = std::find(gen_degs.begin(), gen_degs.end(), b);
                int index1 = int(p1 - gen_degs.begin());
                int index2 = int(p2 - gen_degs.begin());
                rel += GenExp(index1, 1) * GenExp(index2, 1);
            }
            rels.push_back(std::move(rel));
        }
    }

    Gb gb(alg::DEG_MAX, alg::CRI_ON);
    std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) { return p1.GetDeg(gen_degs) < p2.GetDeg(gen_degs); });
    alg::AddRels(gb, std::move(rels), alg::DEG_MAX, gen_degs);
    size_t answer = 402;
    std::cout << "new: " << gb.size() << "==" << answer << '\n';
}

namespace mydetail {
template <typename T, std::size_t N>
constexpr void QuickSort(std::array<T, N>& array, std::size_t low, std::size_t high)
{
    if (high <= low)
        return;
    auto i = low, j = high + 1;
    auto key = array[low];
    for (;;) {
        while (array[++i] < key && i < high)
            ;
        while (key < array[--j] && j > low)
            ;
        if (i >= j)
            break;
        auto tmp = array[i];
        array[i] = array[j];
        array[j] = tmp;
    }

    auto tmp = array[low];
    array[low] = array[j];
    array[j] = tmp;

    if (j > 0)
        QuickSort(array, low, j - 1);
    QuickSort(array, j + 1, high);
}

template <typename T, std::size_t N>
constexpr std::array<T, N> QuickSort(std::array<T, N> array)
{
    QuickSort(array, 0, N - 1);
    return array;
}
}  // namespace mydetail