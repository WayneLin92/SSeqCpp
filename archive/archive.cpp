

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