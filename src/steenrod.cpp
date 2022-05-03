#include "steenrod.h"
#include "benchmark.h"  ////
#include "myio.h"

namespace steenrod {

/* Return the maximum number `result` such that
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

void MulMilnor(const std::array<int, XI_MAX>& R, const std::array<int, XI_MAX>& S, MMilnor1d& result)
{
    constexpr size_t N = 7;  // Support up to t=254

    if (R[N - 1] & S[N - 1])
        return;

    std::array<int, (N + 1) * (N + 1)> X, XR, XS, XT;
    for (size_t i = 1; i <= N; ++i)
        XR[(N - i) * (N + 1) + i] = R[i - 1];
    for (size_t i = 1; i <= N - 1; ++i)
        XS[i * (N + 1) + (N - i)] = S[i - 1];
    for (size_t i = 1; i <= N - 1; ++i)
        XT[i * (N + 1)] = 0;
    X[N] = R[N - 1];
    XT[(N - 1) * (N + 1) + 1] = R[N - 1] | S[N - 1];

    size_t i = N - 1, j = 1;
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
                    result.push_back(MMilnor::Xi(XT.data() + 1));
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
            size_t index = i * (N + 1) + j;
            do {
                if (i + j < N) {
                    ++j;
                    ++index;
                }
                else {
                    ++i;
                    index += (N + 2) - j;
                    j = 1;
                }
            } while (i < N && X[index] == 0);
            if (i >= N)
                break;
            decrease = true;
        }
    }
}

void MulMilnorV2(const std::array<int, XI_MAX>& R, const std::array<int, XI_MAX>& S, MMilnor1d& result)
{
    constexpr size_t N = 7;  // Support up to t=254

    std::array<int, (N + 1) * (N + 1)> X, XR, XS, XT;

    size_t R_floor[N + 1];
    R_floor[0] = 0;
    for (size_t row = 1; row <= N; ++row) {
        if (R[row - 1]) {
            R_floor[row] = row;
            XR[row * (N + 1) + (N - row)] = R[row - 1];
        }
        else
            R_floor[row] = R_floor[row - 1];
    }
    if (R_floor[N] == 0) { /* R = 1 */
        result.push_back(MMilnor::Xi(S.data()));
        return;
    }

    /* The minimal row number r corresponding to R[r - 1] != 0 */
    size_t r_min = 0, r_max = R_floor[N];
    while (r_min < N && R[r_min++] == 0)
        ;
    for (size_t col = 1; col <= N - r_min; ++col)
        XS[R_floor[N - col] * (N + 1) + col] = S[col - 1];
    for (size_t row = 1; row <= N; ++row)
        XT[R_floor[row] * (N + 1) + row - R_floor[row]] = 0;
    for (size_t col = (N + 1) - r_min; col <= N; ++col) {
        X[col] = S[col - 1];
        for (size_t row = r_max; row > 0; row = R_floor[row - 1]) {
            if (col >= row) {
                XT[row * (N + 1) + col - row] = S[col - 1];
                break;
            }
        }
    }

    size_t i = r_max, j = N - i;
    bool decrease = false;
    while (true) {
        bool move_right = false;
        if (j) {
            size_t index = size_t(i * (N + 1) + j);
            size_t index_up = size_t(R_floor[i - 1] * (N + 1) + j);
            size_t index_up_right = size_t(R_floor[i - 1] * N + j + i);
            size_t index_left = index - 1;
            if (i == r_min) {
                if (decrease) {
                    X[index] = (X[index] - 1) & ~(XT[index] | X[index_up_right]);
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index] >> j, XS[index]), XT[index] | X[index_up_right]);
                X[index_up] = XS[index] - X[index];
                if (j >= r_min && (X[index_up] & XT[index - r_min])) {
                    if (X[index])
                        decrease = true;
                    else
                        move_right = true;
                }
                else {
                    XR[index_left] = XR[index] - (X[index] << j);
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
                XR[index_left] = XR[index] - (X[index] << j);
                XS[index_up] = XS[index] - X[index];
                XT[index_up_right] = XT[index] | X[index];
                --j;
            }
        }
        else {
            if (i == r_min) {
                if (!(XR[i * (N + 1)] & X[i])) { /* Add to result. */
                    XT[i] = XR[i * (N + 1)] | X[i];
                    for (size_t i1 = 1; i1 < i; ++i1)
                        XT[i1] = X[i1];
                    result.push_back(MMilnor::Xi(XT.data() + 1));
                }
                move_right = true;
            }
            else {
                size_t index = size_t(i * (N + 1));
                size_t index_up_right = size_t(R_floor[i - 1] * N + i);
                if (XT[index] & XR[index])  //
                    move_right = true;
                else {
                    XT[index_up_right] = XR[index];
                    i = R_floor[i - 1];
                    j = N - i;
                }
            }
        }
        if (move_right) {
            size_t index = i * (N + 1) + j;
            do {
                if (i + j < N) {
                    ++j;
                    ++index;
                }
                else {
                    do {
                        ++i;
                    } while (i <= r_max && R_floor[i] != i);
                    if (i > r_max)
                        break;
                    j = 1;
                    index = i * (N + 1) + 1;
                }
            } while (X[index] == 0);
            if (i > r_max || i + j > N)
                break;
            decrease = true;
        }
    }
}

/* Milnor's multiplication formula.
 * `result.data` is unordered and may contain duplicates.
 */
void MulMilnor(MMilnor lhs, MMilnor rhs, Milnor& result)
{
    auto R = lhs.ToXi();
    auto S = rhs.ToXi();
    int nonzeroes = 0;
    for (int i : R)
        if (i)
            ++nonzeroes;
    if (nonzeroes <= 3)
        MulMilnorV2(R, S, result.data);
    else
        MulMilnor(R, S, result.data);
}

void MulMay(MMilnor lhs, MMilnor rhs, Milnor& result)  ////
{
    MulMilnor(lhs, rhs, result);
    auto w_may = lhs.w_may() + rhs.w_may();
    ut::RemoveIf(result.data, [w_may](MMilnor m) { return m.w_may() != w_may; });
}

/**
 * Sort the sequence and each time remove a pair of identical elements
 */
void SortMod2(MMilnor1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i] = MMilnor{MMILNOR_NULL};
            data[++i] = MMilnor{MMILNOR_NULL};
        }
    ut::RemoveIf(data, [](const MMilnor& m) { return m == MMilnor{MMILNOR_NULL}; });
}

/**
 * Sort the sequence and each time remove a pair of identical elements
 */
void SortMod2(MMod1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i] = MMod(MMILNOR_NULL);
            data[++i] = MMod(MMILNOR_NULL);
        }
    ut::RemoveIf(data, [](const MMod& m) { return m == MMod(MMILNOR_NULL); });
}

Milnor Milnor::operator*(const Milnor& rhs) const
{
    Milnor result;
    for (MMilnor R : this->data)
        for (MMilnor S : rhs.data)
            MulMilnor(R, S, result);
    SortMod2(result.data);
    return result;
}

std::string MMilnor::StrXi() const
{
    auto xi = ToXi();
    auto xi_end = xi.end();
    while (xi_end != xi.begin() && *(xi_end - 1) == 0)
        --xi_end;
    return myio::TplStrCont("Sq(", ",", ")", "1", xi.begin(), xi_end, [](int r) { return std::to_string(r); });
}

std::string MMilnor::Str() const
{
    std::string result;
    for (int i : *this)
        result += "P_{" + std::to_string(MMILNOR_GEN_I[i]) + std::to_string(MMILNOR_GEN_J[i]) + '}';
    return result;
}

std::string Milnor::StrXi() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMilnor m) { return m.StrXi(); });
}

std::string Milnor::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMilnor m) { return m.Str(); });
}

Mod mulLF(MMilnor m, const Mod& x)
{
    Mod result;
    for (MMod mx : x.data)
        if (!(m.data() & mx.e()))
            result.data.push_back(mulLF(m, mx));
    return result;
}

std::string MMod::Str() const
{
    std::string result;
    for (int i : m_no_weight())
        result += "P_{" + std::to_string(MMILNOR_GEN_I[i]) + std::to_string(MMILNOR_GEN_J[i]) + '}';
    result += "v_{" + std::to_string(v()) + '}';
    return result;
}

std::string MMod::StrXi() const
{
    auto s = m_no_weight().StrXi();
    return (s == "1" ? "" : s) + "v_{" + std::to_string(v()) + '}';
}

Mod operator*(MMilnor m, const Mod& x)
{
    Mod result;
    result.data.reserve(150);
    Milnor tmp;
    for (MMod m2 : x.data) {
        tmp.data.clear();
        MulMilnor(m, m2.m_no_weight(), tmp);
        auto v_raw = m2.v_raw();
        for (MMilnor m : tmp.data)
            result.data.push_back(MMod(m.data() + v_raw));
    }
    SortMod2(result.data);
    return result;
}

Mod MulMay(MMilnor m, const Mod& x)  ////
{
    Mod result;
    result.data.reserve(150);
    Milnor tmp;
    for (MMod m2 : x.data) {
        tmp.data.clear();
        MulMay(m, m2.m(), tmp);
        auto v_raw = m2.v_raw();
        for (MMilnor m : tmp.data)
            result.data.push_back(MMod(m.data() + v_raw));
    }
    SortMod2(result.data);
    return result;
}

void mul(MMilnor m, const Mod& x, Milnor& tmp, Mod& result)
{
    result.data.clear();
    for (MMod m2 : x.data) {
        tmp.data.clear();
        MulMilnor(m, m2.m_no_weight(), tmp);
        auto v_raw = m2.v_raw();
        for (MMilnor m : tmp.data)
            result.data.push_back(MMod(m.data() + v_raw));
    }
    SortMod2(result.data);
}

Mod& Mod::iaddmul(MMilnor m, const Mod& x, Milnor& tmp_a, Mod& tmp_m1, Mod& tmp_m2)
{
    mul(m, x, tmp_a, tmp_m1); /* `tmp_m1 = m * x` */
    tmp_m2.data.clear();
    std::swap(data, tmp_m2.data);
    std::set_symmetric_difference(tmp_m1.data.cbegin(), tmp_m1.data.cend(), tmp_m2.data.cbegin(), tmp_m2.data.cend(), std::back_inserter(data));
    return *this;
}

Mod& Mod::iaddmulMay(MMilnor m, const Mod& x, Mod& tmp)
{
    Mod mx = MulMay(m, x);  ////
    tmp.data.clear();
    std::swap(data, tmp.data);
    std::set_symmetric_difference(mx.data.cbegin(), mx.data.cend(), tmp.data.cbegin(), tmp.data.cend(), std::back_inserter(data));
    return *this;
}

std::string Mod::StrXi() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMod m) { return m.StrXi(); });
}

std::string Mod::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMod m) { return m.Str(); });
}

}  // namespace steenrod
