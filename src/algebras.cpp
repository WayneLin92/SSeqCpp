#include "algebras.h"
#include "myexception.h"
#include <iterator>

#ifndef NDEBUG
#include <iostream>
#endif

namespace alg {

Mon mul(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            result.push_back(*k++);
        else if (k->gen > l->gen)
            result.push_back(*l++);
        else {
            result.emplace_back(k->gen, k->exp + l->exp);
            ++k;
            ++l;
        }
    }
    if (k != mon1.end())
        result.insert(result.end(), k, mon1.end());
    else
        result.insert(result.end(), l, mon2.end());
    return result;
}

void detail::mul(const Mon& mon1, const MonOnStack& mon2, Mon& result)
{
    result.clear();
    MonInd k = mon1.begin();
    auto l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            result.push_back(*k++);
        else if (k->gen > l->gen)
            result.push_back(*l++);
        else {
            result.emplace_back(k->gen, k->exp + l->exp);
            ++k;
            ++l;
        }
    }
    if (k != mon1.end())
        result.insert(result.end(), k, mon1.end());
    else
        result.insert(result.end(), l, mon2.end());
}

Mon div(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            result.push_back(*k++);
#ifndef NDEBUG /* DEBUG */
        else if (k->gen > l->gen)
            throw MyException(0x1227de8e, "mon1/mon2 not divisible!\n");
#endif
        else if (k->exp > l->exp) {
            result.emplace_back(k->gen, k->exp - l->exp);
            ++k;
            ++l;
        }
#ifdef NDEBUG
        else {
            ++k;
            ++l;
        }
#else /* DEBUG */
        else if (k->exp == l->exp) {
            ++k;
            ++l;
        }
        else
            throw MyException(0xa9c74ef9, "mon1/mon2 not divisible!\n");
#endif
    }
#ifndef NDEBUG /* DEBUG */
    if (l != mon2.end())
        throw MyException(0x6cdd66bd, "mon1/mon2 not divisible!\n");
    else
#endif
        result.insert(result.end(), k, mon1.end());
    return result;
}

bool divisible(const Mon& mon1, const Mon& mon2)
{
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            return false;
        else if (k->gen > l->gen)
            ++l;
        else if (k->exp > l->exp)
            return false;
        else {
            ++k;
            ++l;
        }
    }
    if (k != mon1.end())
        return false;
    return true;
}

Mon pow(const Mon& m, int e)
{
    if (e == 0)
        return {};
    else if (e == 1)
        return m;
    Mon result;
    for (MonInd p = m.begin(); p != m.end(); ++p)
        result.emplace_back(p->gen, p->exp * e);
    return result;
}

int log(const Mon& mon1, const Mon& mon2)
{
    if (mon2.empty()) {
        /* log with 0 base */
        throw "f50d7f56";
    }
    int q = -1;

    /* Compute q */
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            ++k;
        else if (k->gen > l->gen) {
            q = 0;
            break;
        }
        else if (k->exp < l->exp) {
            q = 0;
            break;
        }
        else {
            int q1 = k->exp / l->exp;
            if (q == -1 || q > q1)
                q = q1;
            ++k;
            ++l;
        }
    }
    if (l != mon2.end())
        q = 0;
    return q;
}

Mon GCD(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            ++k;
        else if (k->gen > l->gen)
            ++l;
        else {
            result.emplace_back(k->gen, std::min(k->exp, l->exp));
            ++k;
            ++l;
        }
    }
    return result;
}

Mon LCM(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            result.push_back(*k++);
        else if (k->gen > l->gen)
            result.push_back(*l++);
        else {
            result.emplace_back(k->gen, std::max(k->exp, l->exp));
            ++k;
            ++l;
        }
    }
    if (k != mon1.end())
        result.insert(result.end(), k, mon1.end());
    else
        result.insert(result.end(), l, mon2.end());
    return result;
}

array Poly2Indices(const Mon1d& poly, const Mon1d& basis)
{
    array result;
    for (const Mon& mon : poly) {
        auto p = std::lower_bound(basis.begin(), basis.end(), mon);
#ifndef NDEBUG
        if (p == basis.end() || mon < (*p)) {
            std::cout << "index not found\n";
            throw "178905cf";
        }
#endif
        result.push_back(int(p - basis.begin()));
    }
    return result;
}

Mon1d Indices2Poly(const array& indices, const Mon1d& basis)
{
    Mon1d result;
    for (int i : indices)
        result.push_back(basis[i]);
    return result;
}

}  // namespace alg
