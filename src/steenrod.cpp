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

MMay MMilnor::ToMMay() const
{
    uint64_t result = 0;
    uint64_t weight = 0;
    for (int d = 0; d < data.size(); ++d) {
        for (int n = data[d], i = 0; n; n >>= 1, ++i) {
            if (n & 1) {
                int j = i + d + 1;
                result |= MMay::rawP(i, j);
                weight += 2 * uint64_t(d) + 1;
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
    for (MMay r : this->data)
        for (MMay r1 : rhs.data) {
            TE2d Xs = GetXs(r.ToMMilnor().data, r1.ToMMilnor().data, N, TE{0});
            for (auto& X : Xs) {
                MMilnor T;
                for (size_t n = 1; n < N + 1; ++n)
                    for (size_t i = 0; i <= n; ++i)
                        T.data[n - 1] += X[i][n - i];
                result.data.push_back(T.ToMMay());
            }
        }
    std::sort(result.data.begin(), result.data.end());
    for (size_t i = 0; i + 1 < result.data.size(); ++i) {
        if (result.data[i] == result.data[i + 1]) {
            result.data[i] = MMay(MMAY_NULL);
            result.data[i + 1] = MMay(MMAY_NULL);
            ++i;
        }
    }
    ut::RemoveIf(result.data, [](MMay m) { return m == MMay(MMAY_NULL); });
    return result;
}

May May::mul(MMay rhs) const
{
    May result;
    for (auto& r : this->data) {
        TE2d Xs = GetXs(r.ToMMilnor().data, rhs.ToMMilnor().data, N, TE{0});
        for (auto& X : Xs) {
            MMilnor T;
            for (size_t n = 1; n < N + 1; ++n)
                for (size_t i = 0; i <= n; ++i)
                    T.data[n - 1] += X[i][n - i];
            result.data.push_back(T.ToMMay());
        }
    }
    std::sort(result.data.begin(), result.data.end());
    for (size_t i = 0; i + 1 < result.data.size(); ++i)
        if (result.data[i] == result.data[i + 1]) {
            result.data[i] = MMay{0xffffffffffffffff};
            result.data[i + 1] = MMay{0xffffffffffffffff};
            ++i;
        }
    ut::RemoveIf(result.data, [](const MMay& m) { return m == MMay{0xffffffffffffffff}; });
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
    for (MMod mv : x.data)
        result += Mod(a.mul(mv.m), mv.v);
    return result;
}

}  // namespace steenrod