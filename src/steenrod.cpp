#include "steenrod.h"

namespace steenrod {
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

/* Return the maximum number result such that 
** `result <= upper_bound` and `result & mask == 0`
*/
int max_mask(int upper_bound, int mask)
{
    int m = upper_bound & mask;
    int n = 0;
    while (m >>= 1)
        n |= m;
    return (upper_bound | n) & ~mask;
}

/* Milnor's multiplication formula.
 * `result.data` is unordered and may contain duplicates.
*/
void MulMilnor(MMay lhs, MMay rhs, May& result)
{
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

void SortAndReduce(MMay1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i] = MMay{0xffffffffffffffff};
            data[i + 1] = MMay{0xffffffffffffffff};
            ++i;
        }
    ut::RemoveIf(data, [](const MMay& m) { return m == MMay{0xffffffffffffffff}; });
}

MMay MMilnor::ToMMay() const
{
    uint64_t result = 0;
    uint64_t weight = 0;
    for (int d = 1; d <= data.size(); ++d) {
        for (int n = data[d - 1], i = 0; n; n >>= 1, ++i) {
            if (n & 1) {
                int j = i + d;
                result |= MMay::rawP(i, j);
                weight += 2 * uint64_t(d) - 1;
            }
        }
    }
    return MMay(result + (weight << MMAY_INDEX_NUM));
}

Milnor MMilnor::operator*(const MMilnor& rhs) const
{
    Milnor result;
    TE2d Xs = GetXs(data, rhs.data, N, TE{0});
    for (auto& X : Xs) {
        MMilnor T;
        for (size_t n = 1; n < N + 1; ++n)
            for (size_t i = 0; i <= n; ++i)
                T.data[n - 1] += X[i][n - i];
        result.data.push_back(T);
    }
    std::sort(result.data.begin(), result.data.end());
    for (size_t i = 0; i + 1 < result.data.size(); ++i)
        if (result.data[i].ToMMay() == result.data[i + 1].ToMMay()) {
            result.data[i].data[0] = -1;
            result.data[i + 1].data[0] = -1;
            ++i;
        }
    ut::RemoveIf(result.data, [](const MMilnor& m) { return m.data[0] == -1; });
    return result;
}

May Milnor::ToMay() const
{
    May result;
    for (auto& m : data)
        result.data.push_back(m.ToMMay());
    std::sort(result.data.begin(), result.data.end());
    return result;
}

Milnor May::ToMilnor() const
{
    Milnor result;
    result.data.reserve(data.size());
    for (MMay m : data)
        result.data.push_back(m.ToMMilnor());
    return result;
}

May May::operator*(const May& rhs) const
{
    May result;
    for (MMay R : this->data)
        for (MMay S : rhs.data)
            MulMilnor(R, S, result);
    SortAndReduce(result.data);
    return result;
}

May May::mul(MMay rhs) const
{
    May result;
    for (auto& r : this->data)
        MulMilnor(r, rhs, result);
    SortAndReduce(result.data);
    return result;
}

std::ostream& operator<<(std::ostream& sout, const Milnor& x)
{
    if (x.data.empty()) {
        sout << '0';
        return sout;
    }
    for (auto pr = x.data.begin(); pr != x.data.end(); ++pr) {
        if (pr != x.data.begin())
            sout << '+';
        sout << "Sq(";
        auto pr_end = pr->data.end();
        while (pr_end != pr->data.begin() && *(pr_end - 1) == 0)
            --pr_end;
        for (auto pi = pr->data.begin(); pi != pr_end; ++pi) {
            if (pi != pr->data.begin())
                sout << ',';
            sout << *pi;
        }
        sout << ')';
    }
    return sout;
}

std::ostream& operator<<(std::ostream& sout, const May& x)
{
    if (x.data.empty()) {
        sout << '0';
        return sout;
    }
    for (auto pm = x.data.begin(); pm != x.data.end(); ++pm) {
        if (pm != x.data.begin())
            sout << '+';
        if (!*pm)
            sout << '1';
        else {
            for (int i : *pm)
                sout << "P_{" << MMAY_GEN_I[i] << MMAY_GEN_J[i] << '}';
        }
    }
    return sout;
}

std::ostream& operator<<(std::ostream& sout, const Mod& x)
{
    if (x.data.empty()) {
        sout << '0';
        return sout;
    }
    for (auto pm = x.data.begin(); pm != x.data.end(); ++pm) {
        if (pm != x.data.begin())
            sout << '+';
        for (int i : pm->m)
            sout << "P_{" << MMAY_GEN_I[i] << MMAY_GEN_J[i] << '}';
        sout << "v_{" << pm->v << '}';
    }
    return sout;
}

Mod operator*(const May& a, const Mod& x)
{
    Mod result;
    May prod;
    auto p = x.data.begin();
    while (p != x.data.end()) {
        int v = p->v;
        auto p1 = std::lower_bound(p, x.data.end(), v - 1, [](MMod x, int v) { return x.v > v; });
        for (MMay R : a.data)
            for (; p != p1; ++p)
                MulMilnor(R, p->m, prod);
        SortAndReduce(prod.data);
        for (MMay T : prod.data)
            result.data.push_back(MMod{T, v});
        prod.data.clear();
    }

    return result;
}

}  // namespace steenrod