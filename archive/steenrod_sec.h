/**
 * @file steenrod_sec.h
 * This file defines the secondary steenrod algebra.
 * We use the algorithm of Dexter Chua: https://arxiv.org/abs/2105.07628
 */

#ifndef STEENROD_SEC_H
#define STEENROD_SEC_H

#include "algebras/steenrod.h"

namespace steenrod {

/* Secondary steenrod algebra in the Milnor basis.
 * When 0<=k<l, this is $Y_{k,l}$.
 * When k=-1, this is l*Sq(m).
 */
class MMilnorSec
{
public:
    MMilnor m;
    int k, l;

public:
    constexpr MMilnorSec() : m(), k(-1), l(1) {}
    constexpr MMilnorSec(MMilnor m) : m(m), k(-1), l(1) {}
    constexpr MMilnorSec(int k, int l) : m(), k(k), l(l) {}
    constexpr MMilnorSec(int k, int l, MMilnor m) : m(m), k(k), l(l) {}

    static MMilnorSec Y(int k, int l);

public:
    int deg() const
    {
        return k >= 0 ? m.deg() + (1 << k) + (1 << l) - 1 : m.deg();
    }
    bool Is2Tor() const
    {
        return k >= 0 || (l % 2) == 0;
    }

    /* Warning: This does not compare the coefficients */
    bool operator<(MMilnorSec other) const;
    bool operator==(MMilnorSec other) const
    {
        return k == other.k && l == other.l && m == other.m;
    };

    std::string Str() const;
};
using MMilnorSec1d = std::vector<MMilnorSec>;

struct MilnorSec
{
    MMilnorSec1d data;
    MilnorSec() {}
    MilnorSec(MMilnorSec m) : data({m}) {}

    MilnorSec operator+(const MilnorSec& other) const;
    MilnorSec& iaddP(const MilnorSec& other, MilnorSec& tmp);
    MilnorSec& operator+=(const MilnorSec& other) {
        MilnorSec tmp;
        return iaddP(other, tmp);
    }
    MilnorSec operator*(const MilnorSec& other) const;
    std::string Str() const;
};

void mulP(const MilnorSec& lhs, const MilnorSec& rhs, MilnorSec& result, Milnor& tmp1, Milnor& tmp2);

void A_sec(MMilnor a, MMilnorSec b, Milnor& result, Milnor& tmp1, Milnor& tmp2);

}  // namespace steenrod

#endif
