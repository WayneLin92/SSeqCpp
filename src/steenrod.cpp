#include "steenrod.h"
#include "benchmark.h"  ////
#include "myio.h"

namespace steenrod {

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

void debug_print_X(int X[(XI_MAX + 1) * (XI_MAX + 1)])
{
    constexpr size_t N = XI_MAX;
    for (size_t i = 0; i <= N; ++i) {
        for (size_t j = 0; j <= N; ++j) {
            if (i == 0 && j == 0)
                std::cout << 0;
            else if (i + j >= N)
                std::cout << 0;
            else
                std::cout << X[i * (N + 1) + j];
            std::cout << ' ';
        }
        std::cout << '\n';
    }
    std::cout << '\n';
}

void MulMilnor(std::array<int, XI_MAX> R, std::array<int, XI_MAX> S, Milnor& result)
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
                    uint64_t data = 0;
                    uint64_t w = 0;
                    for (int d = 1; d <= N; ++d) {
                        for (int n = XT[d], i = 0; n; n >>= 1, ++i) {
                            if (n & 1) {
                                int j = i + d;
                                data |= MMilnor::rawP(i, j);
                                w += 2 * uint64_t(d) - 1;
                            }
                        }
                    }
                    result.data.push_back(MMilnor(data + (w << MMILNOR_INDEX_NUM)));
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

void MulMilnorV2(std::array<int, XI_MAX> R, std::array<int, XI_MAX> S, Milnor& result)
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
        uint64_t data = 0;
        uint64_t w = 0;
        for (int d = 1; d <= N; ++d) {
            for (int n = S[d - 1], i = 0; n; n >>= 1, ++i) {
                if (n & 1) {
                    int j = i + d;
                    data |= MMilnor::rawP(i, j);
                    w += 2 * uint64_t(d) - 1;
                }
            }
        }
        result.data.push_back(MMilnor(data + (w << MMILNOR_INDEX_NUM)));
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
                    uint64_t data = 0;
                    uint64_t w = 0;
                    for (int d = 1; d <= N; ++d) {
                        for (int n = XT[d], i_ = 0; n; n >>= 1, ++i_) {
                            if (n & 1) {
                                int j = i_ + d;
                                data |= MMilnor::rawP(i_, j);
                                w += 2 * uint64_t(d) - 1;
                            }
                        }
                    }
                    result.data.push_back(MMilnor(data + (w << MMILNOR_INDEX_NUM)));
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

void SortMod2(MMilnor1d& data);

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
        MulMilnorV2(R, S, result);
    else
        MulMilnor(R, S, result);
}

void SortMod2(MMilnor1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i] = MMilnor{0xffffffffffffffff};
            data[i + 1] = MMilnor{0xffffffffffffffff};
            ++i;
        }
    ut::RemoveIf(data, [](const MMilnor& m) { return m == MMilnor{0xffffffffffffffff}; });
}

void SortMod2(MModCpt1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i] = MModCpt(0);
            data[i + 1] = MModCpt(0);
            ++i;
        }
    ut::RemoveIf(data, [](const MModCpt& m) { return m == MModCpt(0); });
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

Milnor Milnor::mul(MMilnor rhs) const
{
    Milnor result;
    for (auto& r : this->data)
        MulMilnor(r, rhs, result);
    SortMod2(result.data);
    return result;
}

std::ostream& operator<<(std::ostream& sout, const Milnor& x)
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
                sout << "P_{" << MMILNOR_GEN_I[i] << MMILNOR_GEN_J[i] << '}';
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
            sout << "P_{" << MMILNOR_GEN_I[i] << MMILNOR_GEN_J[i] << '}';
        sout << "v_{" << pm->v << '}';
    }
    return sout;
}

Mod operator*(const Milnor& a, const Mod& x)
{
    Mod result;
    Milnor prod;
    auto p = x.data.begin();
    while (p != x.data.end()) {
        int v = p->v;
        auto p1 = std::lower_bound(p, x.data.end(), v - 1, [](MMod x, int v) { return x.v > v; });
        for (MMilnor m : a.data)
            for (; p != p1; ++p)
                MulMilnor(m, p->m, prod);
        SortMod2(prod.data);
        for (MMilnor T : prod.data)
            result.data.push_back(MMod{T, v});
        prod.data.clear();
    }

    return result;
}

std::ostream& operator<<(std::ostream& sout, const ModCpt& x)
{
    if (x.data.empty()) {
        sout << '0';
        return sout;
    }
    for (auto pm = x.data.begin(); pm != x.data.end(); ++pm) {
        if (pm != x.data.begin())
            sout << '+';
        for (int i : pm->m())
            sout << "P_{" << MMILNOR_GEN_I[i] << MMILNOR_GEN_J[i] << '}';
        sout << "v_{" << pm->v() << '}';
    }
    return sout;
}

ModCpt operator*(const Milnor& a, const ModCpt& x)
{
    ModCpt result;
    Milnor prod;
    auto p = x.data.begin();
    while (p != x.data.end()) {
        int v = p->v();
        auto p1 = std::lower_bound(p, x.data.end(), v - 1, [](MModCpt x, int v) { return x.v() > v; });
        for (MMilnor R : a.data)
            for (; p != p1; ++p)
                MulMilnor(R, p->m(), prod);
        SortMod2(prod.data);
        for (MMilnor T : prod.data)
            result.data.push_back(MModCpt(T, v));
        prod.data.clear();
    }

    return result;
}

ModCpt mulLF(MMilnor m, const ModCpt& x)
{
    ModCpt result;
    Milnor prod;
    auto p = x.data.begin();
    while (p != x.data.end()) {
        int v = p->v();
        auto p1 = std::lower_bound(p, x.data.end(), v - 1, [](MModCpt x, int v) { return x.v() > v; });
        for (; p != p1; ++p)
            if (!gcdLF(m, p->m()))
                prod.data.push_back(mulLF(m, p->m()));
        SortMod2(prod.data);
        for (MMilnor T : prod.data)
            result.data.push_back(MModCpt{T, v});
        prod.data.clear();
    }
    return result;
}

std::string MMilnor::StrXi()
{
    auto xi = ToXi();
    auto xi_end = xi.end();
    while (xi_end != xi.begin() && *(xi_end - 1) == 0)
        --xi_end;
    return myio::TplStrCont("Sq(", ",", ")", "1", xi.begin(), xi_end, [](int r) { return std::to_string(r); });
}

}  // namespace steenrod