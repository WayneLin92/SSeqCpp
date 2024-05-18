#include "steenrod.h"
#include "benchmark.h"
#include "myio.h"

namespace steenrod {

/* Return the maximum number `result` such that
** `result <= upper_bound` and `result & mask == 0`
*/
uint32_t max_mask(uint32_t upper_bound, uint32_t mask)
{
    uint32_t m = upper_bound & mask;
    uint32_t n = 0;
    while (m >>= 1)
        n |= m;
    return (upper_bound | n) & ~mask;
}

constexpr size_t DIMX = XI_MAX + 1;
std::string ToStr(uint32_t x)
{
    return x < 10 ? std::to_string(x) : "?";
}

template<typename T>
std::string ToStr(const T& X, const T& XR, const T& XS, const T& XT)
{
    std::string result;
    for (size_t i = 0; i < DIMX; ++i) {
        for (size_t j = 0; j <= XI_MAX - i; ++j) {
            size_t index = i * DIMX + j;
            fmt::format_to(std::back_inserter(result), "({},{},{},{}) ", ToStr(X[index]), ToStr(XR[index]), ToStr(XS[index]), ToStr(XT[index]));
        }
        result += '\n';
    }
    return result;
}

void MulMilnor(const std::array<uint32_t, XI_MAX>& R, const std::array<uint32_t, XI_MAX>& S, MMilnor1d& result_app)
{
    constexpr size_t N = XI_MAX_MULT;
    constexpr size_t N1 = N + 1;

    if (R[N - 1] & S[N - 1])
        return;

    /* XR[i,j], XS[i,j] controls the upper bound of X[i,j] when we consider R(X)=R and S(X)=S.
     * XT[i,j] controls the upper bound of X[i,j] binary digits.
     */
    std::array<uint32_t, N1 * N1 - N> X, XT;  // TODO: use less memory
    std::array<uint32_t, N1 * N1 - N - N1> XR, XS;
    for (size_t i = 1; i <= N - 1; ++i)
        XR[(N - i) * N1 + i - N1] = R[i - 1];
    for (size_t i = 1; i <= N - 1; ++i)
        XS[i * N1 + (N - i) - N1] = S[i - 1];
    for (size_t i = 1; i <= N - 1; ++i)
        XT[i * N1] = 0;
    X[N] = R[N - 1];
    XT[(N - 1) * N1 + 1] = R[N - 1] | S[N - 1];

    size_t i = N - 1, j = 1;
    bool decrease = false;
    while (true) {
        bool move_right = false;
        if (j) {
            const size_t index = i * N1 + j;
            const size_t index1 = index - N1;
            const size_t index_up = (i - 1) * N1 + j;
            const size_t index_up_right = (i - 1) * N1 + j + 1;
            const size_t index_left = i * N1 + j - 1;
            if (i == 1) {
                /* Row 1 is special because X[index_up_right] is determined before X[index] is determined */
                if (decrease) {
                    X[index] = (X[index] - 1) & ~(XT[index] | X[index_up_right]);
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index1] >> i, XS[index1]), XT[index] | X[index_up_right]);
                X[index_up] = XR[index1] - (X[index] << 1);
                if (X[index_up] & XT[index_left]) {
                    if (X[index])
                        decrease = true;
                    else
                        move_right = true;
                }
                else {
                    XS[index_left - N1] = XS[index1] - X[index];
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
                    X[index] = max_mask(std::min(XR[index1] >> i, XS[index1]), XT[index]);
                XS[index_left - N1] = XS[index1] - X[index];
                XR[index_up - N1] = XR[index1] - (X[index] << i);
                XT[index_up_right] = XT[index] | X[index];
                --j;
            }
        }
        else {
            if (i == 1) {
                if (!(XS[N1 - N1] & X[1])) { /* Add to result. */
                    XT[1] = XS[N1 - N1] | X[1];
                    result_app.push_back(MMilnor::Xi(XT.data() + 1));
                }
                move_right = true;
            }
            else {
                const size_t index = i * N1;
                const size_t index1 = index - N1;
                const size_t index_up_right = (i - 1) * N1 + 1;
                XT[index_up_right] = XS[index1];
                j = N - (--i);
            }
        }
        if (move_right) {
            /* Find the next nonzero element. */
            size_t index = i * N1 + j;
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

void MulMilnorV2(const std::array<uint32_t, XI_MAX>& R, const std::array<uint32_t, XI_MAX>& S, MMilnor1d& result_app)
{
    constexpr size_t N = XI_MAX;
    constexpr size_t nCol = N + 1;
    /* XR[i,j], XS[i,j] controls the upper bound of X[i,j] when we consider R(X)=R and S(X)=S.
     * XT[i,j] controls the upper bound of X[i,j] binary digits.
     */
    std::array<uint32_t, nCol * nCol> X, XR, XS, XT;

    if (R[N - 1] & S[N - 1])
        return;

    /* The minimal i such that R[i - 1] != 0 */
    const size_t i_min = [&R]() {
        size_t i = 0;
        while (i < N && R[i++] == 0)
            ;
        return i;
    }();
    if (i_min == N + 1) { /* Sq(R) == 1 */
        result_app.push_back(MMilnor::Xi(S.data()));
        return;
    }
    if (i_min == N) {
        for (size_t i = 0; i < N - 1; ++i)
            X[i] = S[i];
        X[N - 1] = R[N - 1] | S[N - 1];
        result_app.push_back(MMilnor::Xi(X.data()));
        return;
    }
    /* The maximal i such that R[i - 1] != 0 */
    size_t i_max = [&R]() {
        size_t i = N;
        while (i > 1 && R[i - 1] == 0)
            --i;
        return i;
    }();

    /* Initialization */
    {
        for (size_t i = i_min; i <= i_max; ++i)
            XR[i * nCol + (N - i)] = R[i - 1];
    }
    {
        for (size_t j = 1; j < N - i_max; ++j)
            XS[i_max * nCol + j] = S[j - 1];
        for (size_t j = std::max(1ull, N - i_max); j <= N - i_min; ++j)
            XS[(N - j) * nCol + j] = S[j - 1];
    }
    {
        for (size_t j = N - i_min + 1; j <= N; ++j)
            X[j] = S[j - 1];
    }
    {
        for (size_t i = i_min; i <= i_max; ++i)
            XT[i * nCol] = i > N - i_min ? S[i - 1] : 0;
        for (size_t j = 1; i_max + j <= N; ++j)
            XT[i_max * nCol + j] = (i_max + j) > N - i_min ? S[i_max + j - 1] : 0;
    }

    if (i_max == N) {
        i_max = N - 1;
        XT[(N - 1) * nCol + 1] = R[N - 1] | S[N - 1];
    }
    size_t i = i_max, j = N - i;
    size_t index = i * nCol + j;
    bool decrease = false;
    while (true) {
        if (decrease && X[index] == 0) {
            /* Find the next nonzero element. */
            do {
                if (i + j < N) {
                    ++j;
                    ++index;
                }
                else {
                    ++i;
                    index += (nCol + 1) - j;
                    j = 1;
                }
            } while (i < N && X[index] == 0);
            if (i > i_max) {
                // fmt::print("{}\n", ToStr(X, XR, XS, XT, i_min, i_max));
                break;
            }
        }
        const size_t index_left = index - 1;
        if (i > i_min) {
            const size_t index_up = index - nCol;
            const size_t index_up_right = index_up + 1;
            if (decrease) {
                X[index] = (X[index] - 1) & ~XT[index];
                decrease = false;
            }
            else
                X[index] = max_mask(std::min(XR[index] >> j, XS[index]), XT[index]);
            XS[index_up] = XS[index] - X[index];
            XT[index_up_right] = XT[index] | X[index];

            if (j > 1) {
                XR[index_left] = XR[index] - (X[index] << j);
                --j;
                --index;
            }
            else {
                const auto X_index_left = XR[index] - (X[index] << j);
                if (XT[index_left] & X_index_left)
                    decrease = true;
                else {
                    XT[index_up] = XT[index_left] | (XR[index] - (X[index] << j));
                    index -= i + 1;
                    j = N - (--i);
                }
            }
        }
        else {
            /* Row 1 is special because X[index_up_right] is determined before X[index] is determined */
            const size_t index_up = j;
            const size_t index_up_right = j + i_min;
            if (decrease) {
                X[index] = (X[index] - 1) & ~(XT[index] | X[index_up_right]);
                decrease = false;
            }
            else
                X[index] = max_mask(std::min(XR[index] >> j, XS[index]), XT[index] | X[index_up_right]);
            X[index_up] = XS[index] - X[index];
            if (j >= i_min && (X[index_up] & XT[index - i_min])) {
                decrease = true;
            }
            else {
                XT[index_up_right] = XT[index] | X[index] | X[index_up_right];
                if (j > 1) {
                    XR[index_left] = XR[index] - (X[index] << j);
                    --j;
                    --index;
                }
                else {
                    auto XR_index_left = XR[index] - (X[index] << 1);
                    if (!(XR_index_left & X[i_min])) { /* Add to result. */
                        XT[i_min] = XR_index_left | X[i_min];
                        for (size_t i1 = 1; i1 < i_min; ++i1)
                            XT[i1] = X[i1];
                        // fmt::print("{}\n", ToStr(X, XR, XS, XT, i_min, i_max));
                        result_app.push_back(MMilnor::Xi(XT.data() + 1));
                    }
                    decrease = true;
                }
            }
        }
    }
}

void MulMilnorV3(const std::array<uint32_t, XI_MAX>& R, const std::array<uint32_t, XI_MAX>& S, MMilnor1d& result)
{
    constexpr size_t N = XI_MAX_MULT;
    std::array<uint32_t, (N + 1) * (N + 1)> X, XR, XS, XT;

    /* R_floor[i] */
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
    if (R_floor[N] == 0) { /* R = Sq(0) */
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
            const size_t index = i * (N + 1) + j;
            const size_t index_up = R_floor[i - 1] * (N + 1) + j;
            const size_t index_up_right = R_floor[i - 1] * N + j + i;
            const size_t index_left = index - 1;
            if (i == r_min) {
                if (decrease) {
                    X[index] = (X[index] - 1) & ~(XT[index] | X[index_up_right]);
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index] >> j, XS[index]), XT[index] | X[index_up_right]);
                uint32_t X_index = X[index];
                X[index_up] = XS[index] - X_index;
                if (j >= r_min && (X[index_up] & XT[index - r_min])) {
                    if (X_index)
                        decrease = true;
                    else
                        move_right = true;
                }
                else {
                    XR[index_left] = XR[index] - (X_index << j);
                    XT[index_up_right] = XT[index] | X_index | X[index_up_right];
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
                uint32_t X_index = X[index];
                XR[index_left] = XR[index] - (X_index << j);
                XS[index_up] = XS[index] - X_index;
                XT[index_up_right] = XT[index] | X_index;
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
                uint32_t XT_index = XT[i * (N + 1)];
                uint32_t XR_index = XR[i * (N + 1)];
                if (XT_index & XR_index)
                    move_right = true;
                else {
                    /* Assign XT[index_up_right]. Fixed! */
                    XT[R_floor[i - 1] * N + i] = XT_index | XR_index;
                    i = R_floor[i - 1];
                    j = N - i;
                }
            }
        }
        if (move_right) {
            do {
                if (i + j < N)
                    ++j;
                else {
                    do {
                        ++i;
                    } while (i <= r_max && R_floor[i] != i);
                    if (i > r_max)
                        break;
                    j = 1;
                }
            } while (X[i * (N + 1) + j] == 0);
            if (i > r_max || i + j > N)
                break;
            decrease = true;
        }
    }
}

void MulMilnorV3(const std::array<uint32_t, XI_MAX>& R, const MMod s, MMod1d& result)
{
    constexpr size_t N = XI_MAX_MULT;
    auto S = s.m_no_weight().ToXi();
    std::array<uint32_t, (N + 1) * (N + 1)> X, XR, XS, XT;

    /* R_floor[i] */
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
    if (R_floor[N] == 0) { /* R = Sq(0) */
        result.push_back(MMod(MMilnor::Xi(S.data()).data() | s.v_raw()));
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
            const size_t index = i * (N + 1) + j;
            const size_t index_up = R_floor[i - 1] * (N + 1) + j;
            const size_t index_up_right = R_floor[i - 1] * N + j + i;
            const size_t index_left = index - 1;
            if (i == r_min) {
                if (decrease) {
                    X[index] = (X[index] - 1) & ~(XT[index] | X[index_up_right]);
                    decrease = false;
                }
                else
                    X[index] = max_mask(std::min(XR[index] >> j, XS[index]), XT[index] | X[index_up_right]);
                uint32_t X_index = X[index];
                X[index_up] = XS[index] - X_index;
                if (j >= r_min && (X[index_up] & XT[index - r_min])) {
                    if (X_index)
                        decrease = true;
                    else
                        move_right = true;
                }
                else {
                    XR[index_left] = XR[index] - (X_index << j);
                    XT[index_up_right] = XT[index] | X_index | X[index_up_right];
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
                uint32_t X_index = X[index];
                XR[index_left] = XR[index] - (X_index << j);
                XS[index_up] = XS[index] - X_index;
                XT[index_up_right] = XT[index] | X_index;
                --j;
            }
        }
        else {
            if (i == r_min) {
                if (!(XR[i * (N + 1)] & X[i])) { /* Add to result. */
                    XT[i] = XR[i * (N + 1)] | X[i];
                    for (size_t i1 = 1; i1 < i; ++i1)
                        XT[i1] = X[i1];
                    result.push_back(MMod(MMilnor::Xi(XT.data() + 1).data() | s.v_raw()));
                }
                move_right = true;
            }
            else {
                uint32_t XT_index = XT[i * (N + 1)];
                uint32_t XR_index = XR[i * (N + 1)];
                if (XT_index & XR_index)
                    move_right = true;
                else {
                    /* Assign XT[index_up_right]. Fixed! */
                    XT[R_floor[i - 1] * N + i] = XT_index | XR_index;
                    i = R_floor[i - 1];
                    j = N - i;
                }
            }
        }
        if (move_right) {
            do {
                if (i + j < N)
                    ++j;
                else {
                    do {
                        ++i;
                    } while (i <= r_max && R_floor[i] != i);
                    if (i > r_max)
                        break;
                    j = 1;
                }
            } while (X[i * (N + 1) + j] == 0);
            if (i > r_max || i + j > N)
                break;
            decrease = true;
        }
    }
}

/* Milnor's multiplication formula.
 * `result.data` is unordered and may contain duplicates.
 */
void MulMilnor(MMilnor lhs, MMilnor rhs, Milnor& result_app)
{
    auto R = lhs.ToXi();
    auto S = rhs.ToXi();
    MulMilnorV3(R, S, result_app.data);
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

void mulP(const Milnor& lhs, const Milnor& rhs, Milnor& result)
{
    for (MMilnor R : lhs.data)
        for (MMilnor S : rhs.data)
            MulMilnor(R, S, result);
    SortMod2(result.data);
}

std::string MMilnor::Str() const
{
    auto xi = ToXi();
    auto xi_end = xi.end();
    while (xi_end != xi.begin() && *(xi_end - 1) == 0)
        --xi_end;
    return myio::TplStrCont("Sq(", ",", ")", "1", xi.begin(), xi_end, [](int r) { return std::to_string(r); });
}

std::string MMilnor::StrP() const
{
    std::string result;
    for (int i : *this)
        result += "P_{" + std::to_string(MMILNOR_GEN_I[i]) + std::to_string(MMILNOR_GEN_J[i]) + '}';
    return result;
}

std::string Milnor::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMilnor m) { return m.Str(); });
}

std::string Milnor::StrP() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMilnor m) { return m.StrP(); });
}

Mod mulLF(MMilnor m, const Mod& x)
{
    Mod result;
    for (MMod mx : x.data)
        if (!(m.data() & mx.e()))
            result.data.push_back(mulLF(m, mx));
    return result;
}

std::string MMod::StrP() const
{
    std::string result;
    for (int i : m_no_weight())
        result += "P_{" + std::to_string(MMILNOR_GEN_I[i]) + std::to_string(MMILNOR_GEN_J[i]) + '}';
    result += "v_{" + std::to_string(v()) + '}';
    return result;
}

std::string MMod::Str() const
{
    auto s = m_no_weight().Str();
    return (s == "1" ? "" : s) + "v_{" + std::to_string(v()) + '}';
}

void MulMayP(MMilnor mon, const Mod& x, Mod& result)
{
    result.data.clear();
    auto R = mon.ToXi();
    for (MMod mx : x.data) {
        auto old_size = result.data.size();
        MulMilnorV3(R, mx, result.data);
        auto w_may = mon.w_may() + mx.w_may();
        result.data.erase(std::remove_if(result.data.begin() + old_size, result.data.end(), [w_may](MMod m) { return m.w_may() != w_may; }), result.data.end());
    }
    SortMod2(result.data);
}

Mod MulMay(MMilnor m, const Mod& x)
{
    Mod result;
    MulMayP(m, x, result);
    return result;
}

void mulP(MMilnor m, const Mod& x, Mod& result)
{
    result.data.clear();
    auto R = m.ToXi();
    for (MMod m1 : x.data)
        MulMilnorV3(R, m1, result.data);
    SortMod2(result.data);
}

Mod& Mod::iaddmulP(MMilnor m, const Mod& x, Mod& tmp_x1, Mod& tmp_x2)
{
    mulP(m, x, tmp_x1); /* `tmp_x1 = m * x` */
    tmp_x2.data.clear();
    std::set_symmetric_difference(tmp_x1.data.cbegin(), tmp_x1.data.cend(), data.cbegin(), data.cend(), std::back_inserter(tmp_x2.data));
    ut::copy(tmp_x2.data, data);
    return *this;
}

Mod& Mod::iaddmulMay(MMilnor m, const Mod& x, Mod& tmp)
{
    Mod mx = MulMay(m, x);  ////
    tmp.data.clear();
    std::set_symmetric_difference(data.cbegin(), data.cend(), mx.data.cbegin(), mx.data.cend(), std::back_inserter(tmp.data));
    ut::copy(tmp.data, data);
    return *this;
}

std::string Mod::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMod m) { return m.Str(); });
}

std::string Mod::StrP() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](MMod m) { return m.StrP(); });
}

}  // namespace steenrod
