

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

/****
*/

constexpr size_t N = MILNOR_BUFFER_SIZE;
using Tuple = MMilnor::TypeData;
using TE = std::array<int, N + 1>;
using TE1d = std::vector<TE>;
using TE2d = std::vector<TE1d>;

/* Return (x_0,x_1,...) such that $\sum 2^jx_j=r$ and $x_i<=s_{i-1}$ and not x_j & b_j */
void GetRows(TE1d& result, int r, const Tuple& S, size_t j_max, const TE& B)
{
    if (j_max == 0) {
        result.push_back(TE{r}); /* return (r,0,...) */
        return;
    }
    int s = S[j_max - 1];
    int x_j_max = std::min(s, r >> j_max);
    for (int x_j = x_j_max; x_j >= 0; --x_j) {
        if (!(x_j & B[j_max])) {
            size_t k = result.size();
            GetRows(result, r - (x_j << j_max), S, j_max - 1, B);
            for (; k < result.size(); ++k)
                result[k][j_max] = x_j;
        }
    }
}

TE2d GetXs(const Tuple& R, const Tuple& S, size_t i_max, const TE& B)
{
    TE2d result;
    if (i_max == 0) {
        bool all_true = true;
        for (size_t i = 0; i < N - 1; ++i) {
            if (S[i] & B[i + 1]) {
                all_true = false;
                break;
            }
        }
        if (all_true) {
            TE S1;
            S1[0] = 0;
            for (size_t i = 1; i < N + 1; ++i)
                S1[i] = S[i - 1];
            result.push_back({S1});
        }
        return result;
    }
    TE1d last_rows;
    GetRows(last_rows, R[i_max - 1], S, N, B);
    for (auto& row : last_rows) {
        TE B1;
        B1[0] = 0;
        for (size_t i = 1; i < N + 1; ++i)
            B1[i] = B[i - 1] + row[i - 1];
        Tuple S1;
        for (size_t i = 0; i < N; ++i)
            S1[i] = S[i] - row[i + 1];
        TE2d result1 = GetXs(R, S1, i_max - 1, B1);
        for (auto& X : result1) {
            X.push_back(row);
            result.push_back(std::move(X));
        }
    }
    return result;
}