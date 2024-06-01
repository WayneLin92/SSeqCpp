#include "steenrodR.h"
#include "main.h"
#include "algebras/myio.h"
#include "algebras/database.h"
#include <fmt/format.h>
#include <array>
#include <vector>

namespace steenrodR {

using ArrMB = std::array<uint32_t, XI_MAX>;

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

void MulMilnor(const ArrMB& R, const ArrMB& S, std::vector<ArrMB>& result_app)
{
    constexpr size_t N = XI_MAX;
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
                    result_app.push_back({});
                    for (size_t i = 0; i < N; ++i)
                        result_app.back()[i] = XT[i + 1];
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
}



int main_test(int argc, char** argv, int& index, const char* desc) {
    std::vector<int> R, S;
    int a = -1, b = -2;

    myio::CmdArg1d args = {{"R", &R}, {"S", &S}, {"a", &a}};
    myio::CmdArg1d op_args = {{"b", &b}};
    if (int error = myio::ParseArguments(argc, argv, index, PROGRAM, desc, VERSION, args, op_args))
        return error;

    fmt::print("R={}\n", myio::Serialize(R));
    fmt::print("S={}\n", myio::Serialize(S));
    fmt::print("a={}, b={}\n", a, b);

    return 0;
}
