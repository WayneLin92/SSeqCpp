#include "steenrod_sec.h"
#include "algebras/myio.h"

template <>
struct fmt::formatter<std::array<uint32_t, (steenrod::XI_MAX_MULT + 1) * (steenrod::XI_MAX_MULT + 1)>>
{
    using T = std::array<uint32_t, (steenrod::XI_MAX_MULT + 1) * (steenrod::XI_MAX_MULT + 1)>;
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template <typename FormatContext>
    auto format(const T& X, FormatContext& ctx)
    {
        for (size_t i = 0; i <= steenrod::XI_MAX_MULT; ++i) {
            for (size_t j = 0; j <= steenrod::XI_MAX_MULT - i; ++j) {
                auto n = X[i * (steenrod::XI_MAX_MULT + 1) + j];
                if (n < 100)
                    fmt::format_to(ctx.out(), "{} ", n);
                else
                    fmt::format_to(ctx.out(), "? ");
            }
            fmt::format_to(ctx.out(), "\n");
        }
        return ctx.out();
    }
};

namespace steenrod {
/* Compute the contraction C(\xi_i^{2^k}, Sq(R)) inplace */
int Contr(uint32_t i, uint32_t k, std::array<uint32_t, XI_MAX>& R)
{
    if (i == 0)
        return 1;
    --i;
    if (R[i] >> k) {
        R[i] -= 1 << k;
        return 1;
    }
    return 0;
}

/* Compute the contraction C(\xi_i^{2^k}, Sq(R)) */
int Contr(uint32_t i, uint32_t k, const std::array<uint32_t, XI_MAX>& R, std::array<uint32_t, XI_MAX>& result)
{
    result = R;
    return Contr(i, k, result);
}

/* Compute the contraction C(\xi_i^{2^k}\xi_j^{2^l}, m)
 */
int Contr(uint32_t i, uint32_t k, uint32_t j, uint32_t l, const std::array<uint32_t, XI_MAX>& R, std::array<uint32_t, XI_MAX>& result)
{
    if (Contr(i, k, R, result))
        return Contr(j, l, result);
    return 0;
}

/* Compute the contraction C(\xi_i^{2^k}, m)*/
int Contr(uint32_t i, uint32_t k, MMilnor m, MMilnor& result)
{
    auto xi = m.ToXi();
    if (Contr(i, k, xi)) {
        result = MMilnor::Xi(xi.data());
        return 1;
    }
    return 0;
}

/* Compute the contraction C(\xi_i^{2^k}\xi_j^{2^l}, m)
 */
int Contr(uint32_t i, uint32_t k, uint32_t j, uint32_t l, MMilnor m, MMilnor& result)
{
    auto xi = m.ToXi();
    if (Contr(i, k, xi)) {
        if (Contr(j, l, xi)) {
            result = MMilnor::Xi(xi.data());
            return 1;
        }
    }
    return 0;
}

MMilnorSec MMilnorSec::Y(int k, int l)
{
    if (k < l)
        return MMilnorSec(k, l);
    else if (k > l)
        return MMilnorSec(l, k);
    else
        return MMilnorSec(MMilnor::P(0, k + 1));
}

bool MMilnorSec::operator<(MMilnorSec other) const
{
    if (k >= 0 && other.k >= 0) {
        auto d1 = (1 << k) + (1 << l);
        auto d2 = (1 << other.k) + (1 << other.l);
        if (d1 < d2)
            return true;
        if (d1 > d2)
            return false;
        if (k < other.k)
            return true;
        if (k > other.l)
            return false;
        return m < other.m;
    }
    if (k < 0 && other.k < 0)
        return m < other.m;
    return k < 0;
};

std::string MMilnorSec::Str() const
{
    std::string result;
    if (k < 0) {
        if (l == 1)
            result = m.Str();
        else if (m)
            result = fmt::format("{}{}", l, m.Str());
        else
            result = fmt::format("{}", l);
    }
    else {
        if (m)
            result = fmt::format("Y_{{{},{}}}{}", k, l, m.Str());
        else
            result = fmt::format("Y_{{{},{}}}", k, l);
    }
    return result;
}

std::string MilnorSec::Str() const
{
    auto data1 = data;
    std::sort(data1.begin(), data1.end(), [](MMilnorSec m1, MMilnorSec m2) { return m1.Is2Tor() < m2.Is2Tor() || m1.Is2Tor() == m2.Is2Tor() && m1 < m2; });
    return myio::TplStrCont("", "+", "", "0", data1.begin(), data1.end(), [](MMilnorSec m) { return m.Str(); });
}

MilnorSec& MilnorSec::iaddP(const MilnorSec& other, MilnorSec& tmp)
{
    tmp.data.clear();
    auto first1 = data.begin(), last1 = data.end();
    auto first2 = other.data.begin(), last2 = other.data.end();
    auto p_result = std::back_inserter(tmp.data);
    while (first1 != last1) {
        if (first2 == last2) {
            std::copy(first1, last1, p_result);
            break;
        }
        if (*first1 < *first2) {
            *p_result = *first1;
            ++p_result;
            ++first1;
        }
        else {
            if (*first2 < *first1) {
                *p_result = *first2;
                ++p_result;
            }
            else {
                if (first1->k < 0) {
                    auto l = (first1->l + first2->l) % 4;
                    if (l) {
                        *p_result = MMilnorSec(-1, l, first1->m);
                        ++p_result;
                    }
                }
                ++first1;
            }
            ++first2;
        }
    }
    ut::copy(tmp.data, data);
    std::copy(first2, last2, std::back_inserter(data));
    return *this;
}

MilnorSec MilnorSec::operator+(const MilnorSec& other) const
{
    MilnorSec result;
    auto first1 = data.begin(), last1 = data.end();
    auto first2 = other.data.begin(), last2 = other.data.end();
    auto p_result = std::back_inserter(result.data);
    while (first1 != last1) {
        if (first2 == last2) {
            std::copy(first1, last1, p_result);
            break;
        }
        if (*first1 < *first2) {
            *p_result = *first1;
            ++p_result;
            ++first1;
        }
        else {
            if (*first2 < *first1) {
                *p_result = *first2;
                ++p_result;
            }
            else {
                if (first1->k < 0) {
                    auto l = (first1->l + first2->l) % 4;
                    if (l) {
                        *p_result = MMilnorSec(-1, l, first1->m);
                        ++p_result;
                    }
                }
                ++first1;
            }
            ++first2;
        }
    }
    std::copy(first2, last2, p_result);
    return result;
}

/* Binomial coefficient modulo 4 */
uint32_t BinomMod4(uint32_t n, uint32_t m);

void MulMilnor(const std::array<uint32_t, XI_MAX>& R, const std::array<uint32_t, XI_MAX>& S, MMilnor1d& result);

void MulMilnorMod4(uint32_t l, const std::array<uint32_t, XI_MAX>& R, const std::array<uint32_t, XI_MAX>& S, MMilnorSec1d& result, MMilnor1d& tmp)
{
    constexpr size_t N = XI_MAX_MULT;

    /* XR[i,j], XS[i,j] controls the upper bound of X[i,j] when we consider R(X)=R and S(X)=S.
     * XT[i,j] is the partial sum of diagonals of X[i,j].
     */
    std::array<uint32_t, (N + 1) * (N + 1)> X, XR, XS, XT, Xb;
    XT[N * (N + 1)] = S[N - 1] + R[N - 1];
    Xb[N * (N + 1)] = BinomMod4(XT[N * (N + 1)], R[N - 1]);
    if (Xb[N * (N + 1)] == 0)
        return;
    X[N] = R[N - 1];

    for (size_t i = 1; i <= N; ++i)
        XR[(N - i) * (N + 1) + i] = R[i - 1];
    for (size_t i = 1; i <= N - 1; ++i)
        XS[i * (N + 1) + (N - i)] = S[i - 1];

    size_t i = N - 1, j = 1;
    bool decrease = false;
    while (true) {
        bool move_right = false;
        if (j) {
            const size_t index = i * (N + 1) + j;
            const size_t index_up = (i - 1) * (N + 1) + j;
            const size_t index_up_right = (i - 1) * (N + 1) + j + 1;
            const size_t index_left = i * (N + 1) + j - 1;
            const size_t index_prev = i + j == N ? (i + 1) * (N + 1) : i * (N + 1) + j + 1;
            const size_t index_bottom_left = (i + 1) * (N + 1) + j - 1;

            if (decrease)
                --X[index];
            else
                X[index] = std::min(XR[index] >> i, XS[index]);
            /* Row 1 is special because X[index_up_right] is determined before X[index] is determined */
            XT[index] = XT[(i > 1 || j == N - 1) ? index_bottom_left : index_up_right] + X[index];

            while (X[index]) {
                uint32_t b = (Xb[index_prev] * BinomMod4(XT[index], X[index])) % 4;
                if (b) {
                    Xb[index] = b;
                    break;
                }
                --X[index];
                --XT[index];
            }
            if (X[index] == 0)
                Xb[index] = Xb[index_prev];
            decrease = false;

            if (i == 1) {
                X[index_up] = XR[index] - (X[index] << 1);
                if (j > 1) {
                    XT[index_up] = XT[(i + 1) * (N + 1) + j - 2] + X[index_up];
                    Xb[index] = (Xb[index] * BinomMod4(XT[index_up], X[index_up])) % 4;
                    if (Xb[index] == 0) {
                        if (X[index])
                            decrease = true;
                        else {
                            // fmt::print("X=\n{}\nXb=\n{}\n", X, Xb);
                            move_right = true;
                        }
                    }
                    else {
                        XS[index_left] = XS[index] - X[index];
                        --j;
                    }
                }
                else {
                    XS[index_left] = XS[index] - X[index];
                    --j;
                }
            }
            else {
                XS[index_left] = XS[index] - X[index];
                XR[index_up] = XR[index] - (X[index] << i);
                --j;
            }
        }
        else {
            if (i == 1) {
                XT[N + 1] = XS[N + 1] + X[1];
                uint32_t b = (Xb[N + 2] * BinomMod4(XT[N + 1], X[1])) % 4;
                // fmt::print("push: X=\n{}\nXS=\n{}\nXT=\n{}\nXb=\n{}\n", X, XS, XT, Xb);
                if (b) {
                    result.push_back(MMilnorSec(-1, (b * l) % 4, MMilnor::Xi(XT.data() + N + 1)));
                }
                move_right = true;
            }
            else {
                const size_t index = i * (N + 1);
                XT[index] = XS[index];
                Xb[index] = Xb[index + 1];
                j = N - (--i);
            }
        }
        if (move_right) {
            /* Find the next nonzero element. */
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

    /* Y part */
    // for (uint32_t k = 0; k <= XI_MAX_MULT - 1; ++k) {
    //     std::array<uint32_t, XI_MAX> S1;
    //     if (Contr(k + 1, 0, S, S1)) {
    //         for (uint32_t m = 0; m < XI_MAX_MULT; ++m) {
    //             for (uint32_t n = m + 1; n <= XI_MAX_MULT; ++n) {
    //                 std::array<uint32_t, XI_MAX> R1;
    //                 if (Contr(m, k, n, k, R, R1)) {
    //                     tmp.clear();
    //                     MulMilnor(R1, S1, tmp);
    //                     for (auto& mm : tmp)
    //                         result.push_back(MMilnorSec(m + k, n + k, mm));
    //                 }
    //             }
    //         }
    //     }
    // }
}

void MulMilnor(MMilnor lhs, MMilnor rhs, Milnor& result_app);

void MulMilnorMod4(MMilnorSec lhs, MMilnorSec rhs, MilnorSec& result_app, Milnor& tmp1, Milnor& tmp2)
{
    if (lhs.Is2Tor()) {
        if (rhs.Is2Tor())
            return;
        else { /* ex: YSq * 3Sq or 2Sq * 3Sq */
            tmp1.data.clear();
            MulMilnor(lhs.m, rhs.m, tmp1);
            for (auto& m : tmp1.data)
                result_app.data.push_back(MMilnorSec(lhs.k, lhs.l, m));
        }
    }
    else {
        if (rhs.Is2Tor()) {
            if (rhs.k >= 0) { /* ex: 3Sq * YSq */
                for (uint32_t i = 0; i <= XI_MAX_MULT; ++i) {
                    for (uint32_t j = 0; j <= XI_MAX_MULT; ++j) {
                        MMilnor m1;
                        if (Contr(i, rhs.k, j, rhs.l, lhs.m, m1)) {
                            uint32_t k = rhs.k + i, l = rhs.l + j;
                            if (k != l) {
                                if (k > l)
                                    std::swap(k, l);
                                tmp1.data.clear();
                                MulMilnor(m1, rhs.m, tmp1);
                                for (auto& m : tmp1.data)
                                    result_app.data.push_back(MMilnorSec(k, l, m));
                            }
                            else {
                                tmp1.data.clear();
                                MulMilnor(MMilnor::P(0, k + 1), m1, tmp1);
                                tmp2.data.clear();
                                for (auto& m : tmp1.data)
                                    MulMilnor(m, rhs.m, tmp2); /* tmp2 = Milnor::P(0, k + 1) * m1 * rhs.m; */
                                for (auto& m : tmp2.data)
                                    result_app.data.push_back(MMilnorSec(-1, 2, m));
                            }
                        }
                    }
                }
            }
            else { /* ex: 3Sq * 2Sq */
                tmp1.data.clear();
                MulMilnor(lhs.m, rhs.m, tmp1);  ////
                for (auto& m : tmp1.data)
                    result_app.data.push_back(MMilnorSec(-1, 2, m));
            }
        }
        else { /* ex: 3Sq * 3Sq */
            auto R = lhs.m.ToXi();
            auto S = rhs.m.ToXi();
            MulMilnorMod4((lhs.l * rhs.l) % 4, R, S, result_app.data, tmp1.data);
        }
    }
}

/**
 * Sort the sequence and each time remove a pair of identical elements
 */
void SortMod4(MMilnorSec1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i) {
        if (data[i].k == data[i + 1].k && data[i].m == data[i + 1].m) {
            if (data[i].k >= 0) {
                if (data[i].l == data[i + 1].l) {
                    data[i].k = -2;
                    data[++i].k = -2;
                }
            }
            else {
                size_t j = i + 1;
                data[i].l += data[j].l;
                data[j].k = -2;
                for (++j; j < data.size() && data[i].k == data[j].k && data[i].m == data[j].m; ++j) {
                    data[i].l += data[j].l;
                    data[j].k = -2;
                }
                data[i].l %= 4;
                if (data[i].l == 0)
                    data[i].k = -2;
                i = j - 1;
            }
        }
    }
    ut::RemoveIf(data, [](const MMilnorSec& m) { return m.k == -2; });
}

void mulP(const MilnorSec& lhs, const MilnorSec& rhs, MilnorSec& result, Milnor& tmp1, Milnor& tmp2)
{
    for (MMilnorSec R : lhs.data)
        for (MMilnorSec S : rhs.data)
            MulMilnorMod4(R, S, result, tmp1, tmp2);
    SortMod4(result.data);
}

MilnorSec MilnorSec::operator*(const MilnorSec& other) const
{
    MilnorSec result;
    Milnor tmp1, tmp2;
    mulP(*this, other, result, tmp1, tmp2);
    return result;
}

void SortMod2(MMilnor1d& data);

/* Return A(a,b) */
void A_sec(MMilnor a, MMilnorSec b, Milnor& result, Milnor& tmp1, Milnor& tmp2)
{
    // fmt::print("a={}, b={}, A(a, b)=", a.StrXi(), b);
    result.data.clear();
    if (!b.Is2Tor()) {
        fmt::print("b should be a two torsion\n");
        std::exit(-2);
    }
    if (b.k < 0) { /* A(a, 2Sq(R)) */
        MMilnor a1;
        if (Contr(1, 0, a, a1)) {
            MulMilnor(a1, b.m, result);
            SortMod2(result.data);
        }
    }
    else { /* A(a, Y_{k,l}Sq(R)) */
        std::array<uint32_t, XI_MAX> A = a.ToXi(), M1;
        for (uint32_t i = 0; i <= XI_MAX_MULT - b.k; ++i) {
            for (uint32_t j = 0; j <= XI_MAX_MULT && b.k + i >= b.l + j; ++j) {
                if (Contr(i, b.k, j, b.l, A, M1)) {
                    uint32_t k = b.k + i, l = b.l + j;
                    std::array<uint32_t, XI_MAX> M = {};
                    ++M[size_t(k - 1)];
                    ++M[size_t(l - 1)];
                    tmp1.data.clear();
                    MulMilnor(M, M1, tmp1.data);  // TODO: optimize by avoiding toXi() and Xi()
                    tmp2.data.clear();
                    for (auto& m_ : tmp1.data)
                        MulMilnor(m_, b.m, tmp2); /* tmp2 = m * m1 * b.m; */
                    SortMod2(tmp2.data);
                    result.iaddP(tmp2, tmp1);
                }
            }
        }
    }
    // fmt::print("{}\n", result.StrXi());
}

}  // namespace steenrod