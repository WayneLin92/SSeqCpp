/** \file steenrod.h
 * A Component for monomials and polynomials over $F_2$.
 * All are incapsulated in the `namespace alg`.
 */

#ifndef STEENROD_H
#define STEENROD_H

#include "myexception.h"
#include <iterator>
#include <array>
#include <vector>
#include <algorithm>
#include <iostream>

namespace alg {
using array = std::vector<int>;

struct SteenrodMilnor;

struct MonSteenrodMilnor
{
    static constexpr size_t BUFFER_SIZE = 9; /* Support degree < 2^(s+1) - 1 */
    using TypeData = std::array<int, BUFFER_SIZE>;
    TypeData data;

    static MonSteenrodMilnor P(int i, int j)  /* P_{ij} */
    {
        MonSteenrodMilnor result{{}};
        result.data[size_t(j - i - 1)] = (1 << i);
        return result;
    }

    bool operator<(const MonSteenrodMilnor& rhs) const
    {
        int w1 = weight(), w2 = rhs.weight();
        if (w1 < w2)
            return true;
        if (w2 < w1)
            return false;
        if (data < rhs.data)
            return true;
        return false;
    };

    int weight() const
    {
        int result = 0;
        for (size_t i = 0; i < BUFFER_SIZE; ++i)
            result += (2 * (int)i + 1) * std::_Popcount(unsigned(data[i]));
        return result;
    }
    int deg() const
    {
        int result = 0;
        for (size_t i = 0; i < BUFFER_SIZE; ++i)
            result += ((1 << (i + 1)) - 1) * data[i];
        return result;
    }

    SteenrodMilnor operator*(const MonSteenrodMilnor& rhs) const;
};
using MonSteenrodMilnor1d = std::vector<MonSteenrodMilnor>;

struct SteenrodMilnor
{
    MonSteenrodMilnor1d data;

    SteenrodMilnor() = default;
    SteenrodMilnor(const MonSteenrodMilnor& r) : data({r}) {}

    static SteenrodMilnor Unit()
    {
        return SteenrodMilnor{{{{
            0,
        }}}};
    }

    SteenrodMilnor operator+(const SteenrodMilnor& rhs) const
    {
        SteenrodMilnor result;
        std::set_symmetric_difference(data.begin(), data.end(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(result.data));
        return result;
    }
    SteenrodMilnor& operator+=(const SteenrodMilnor& rhs)
    {
        SteenrodMilnor tmp;
        std::swap(data, tmp.data);
        std::set_symmetric_difference(tmp.data.cbegin(), tmp.data.cend(), rhs.data.cbegin(), rhs.data.cend(), std::back_inserter(data));
        return *this;
    }
    SteenrodMilnor operator*(const SteenrodMilnor& rhs) const
    {
        SteenrodMilnor result;
        for (auto& r : data)
            for (auto& r1 : rhs.data)
                result += r * r1;
        return result;
    }
};
std::ostream& operator<<(std::ostream& sout, const SteenrodMilnor& x);

namespace detail {
    inline constexpr size_t INDEX_NUM = (MonSteenrodMilnor::BUFFER_SIZE + 1) * (MonSteenrodMilnor::BUFFER_SIZE + 2) / 2 - 1;
    inline constexpr std::array<int, INDEX_NUM> MonSteenrodMayI()
    {
        std::array<int, INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)MonSteenrodMilnor::BUFFER_SIZE + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < INDEX_NUM)
                    result[n++] = i;
        return result;
    }
    inline constexpr std::array<int, INDEX_NUM> MonSteenrodMayJ()
    {
        std::array<int, INDEX_NUM> result = {0};
        size_t n = 0;
        for (int j = 1; j <= (int)MonSteenrodMilnor::BUFFER_SIZE + 1; ++j)
            for (int i = j - 1; i >= 0; --i)
                if (n < INDEX_NUM)
                    result[n++] = j;
        return result;
    }
}
struct MonSteenrodMay
{
    static constexpr size_t INDEX_NUM = detail::INDEX_NUM;
    static constexpr std::array<int, INDEX_NUM> GEN_I = detail::MonSteenrodMayI();
    static constexpr std::array<int, INDEX_NUM> GEN_J = detail::MonSteenrodMayJ();

    array data; /*(i,j,...) means x_ix_j... */

    explicit operator SteenrodMilnor() const;
    MonSteenrodMay() = default;

    static int GenP(int i, int j) /* P_{ij} */
    {
#ifndef NDEBUG
        if (!(0 <= i && i < j))
            throw MyException(0x4bfc0d6fU, "BUG");
#endif  // !NDEBUG
        return j * (j + 1) / 2 - i - 1;
    }
    static MonSteenrodMay P(int i, int j) /* P_{ij} */
    {
        return MonSteenrodMay{{GenP(i, j)}};
    }

    bool operator<(const MonSteenrodMay& rhs) const
    {
        int w1 = weight(), w2 = rhs.weight();
        if (w1 < w2)
            return true;
        if (w2 < w1)
            return false;
        if (data > rhs.data) /* We want Revlex here */
            return true;
        return false;
    };

    int weight() const
    {
        int result = 0;
        for (int index : data)
            result += 2 * (GEN_J[index] - GEN_I[index]) - 1;
        return result;
    }
    int deg() const
    {
        int result = 0;
        for (int index : data)
            result += (1 << GEN_J[index]) - (1 << GEN_I[index]);
        return result;
    }

    static MonSteenrodMay LF(const MonSteenrodMilnor& r) /* Leading form according to weight */
    {
        MonSteenrodMay result;
        for (int d = 0; d < r.data.size(); ++d) {
            int n = r.data[d];
            int i = 0;
            while (n) {
                if (n & 1) {
                    int j = i + d + 1;
                    result.data.push_back(GenP(i, j));
                }
                n >>= 1;
                ++i;
            }
        }
        std::sort(result.data.begin(), result.data.end());
        return result;
    }
};
using MonSteenrodMay1d = std::vector<MonSteenrodMay>;

struct SteenrodMay
{
    MonSteenrodMay1d data;
    SteenrodMay() = default;
    SteenrodMay(const MonSteenrodMay& r) : data({r}) {}
    explicit SteenrodMay(SteenrodMilnor x) {
        while (!x.data.empty()) {
            int w = x.data.front().weight();
            MonSteenrodMay m = MonSteenrodMay::LF(x.data.front());
            data.push_back(m);
            x += SteenrodMilnor(m);
        }
    }
};

std::ostream& operator<<(std::ostream& sout, const SteenrodMay& x);

}  // namespace alg


#endif