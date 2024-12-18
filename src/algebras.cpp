#include "algebras.h"
#include "myexception.h"
#include "myio.h"
#include "utility.h"
#include <iterator>

#ifndef NDEBUG
#include <iostream>
#endif

namespace alg {

std::string GE::Str() const
{
    std::string result = "x_";
    std::string under = std::to_string(g());
    std::string upper = std::to_string(e_masked());
    if (under.size() > 1)
        result += '{' + under + '}';
    else
        result += under;

    if (upper.size() > 1)
        result += "^{" + upper + '}';
    else if (e_masked() != 1)
        result += '^' + upper;
    return result;
}

}  // namespace alg

namespace alg2 {

MonTrace Mon::Trace() const
{
    MonTrace result = 0;
    for (size_t i = 0; i < (size_t)size_; ++i) {
        const int bits_exp1 = 56;
        const int bits_exp2 = 64 - bits_exp1;
        result |= (MonTrace(1) << (data_[i].g() % bits_exp1));
        if (data_[i].e_masked() >= 2)
            result |= (MonTrace(1) << ((data_[i].g() % bits_exp2) + bits_exp1));
    }
    return result;
}

std::string Mon::Str() const
{
    return myio::TplStrCont("", "", "", "1", begin(), end(), [](GE p) { return p.Str(); });
}

Mon operator*(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            result.push_back(*k++);
        else if (k->g_raw() < l->g_raw())
            result.push_back(*l++);
        else {
            result.push_back(GE(k->data + l->e_masked()));
            ++k;
            ++l;
        }
    }
    if (k != mon1.end())
        result.insert(k, mon1.end());
    else
        result.insert(l, mon2.end());
    return result;
}

Mon operator/(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            result.push_back(*k++);
#ifndef NDEBUG /* DEBUG */
        else if (k->g_raw() < l->g_raw())
            throw ErrorIdMsg(0x1227de8e, "mon1/mon2 not divisible!\n");
#endif
        else if (k->e() > l->e()) {
            result.push_back(GE(k->data - l->e_masked()));
            ++k;
            ++l;
        }
#ifdef NDEBUG
        else {
            ++k;
            ++l;
        }
#else /* DEBUG */
        else if (k->e() == l->e()) {
            ++k;
            ++l;
        }
        else
            throw ErrorIdMsg(0xa9c74ef9, "mon1/mon2 not divisible!\n");
#endif
    }
#ifndef NDEBUG /* DEBUG */
    if (l != mon2.end())
        throw ErrorIdMsg(0x6cdd66bd, "mon1/mon2 not divisible!\n");
    else
#endif
        result.insert(k, mon1.end());
    return result;
}

bool divisible(const Mon& mon1, const Mon& mon2)
{
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            return false;
        else if (k->g_raw() < l->g_raw())
            ++l;
        else if (k->e() > l->e())
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

Mon pow(const Mon& mon, int e)
{
    Mon result;
    if (e == 0)
        return result;
    else if (e == 1) {
        result = mon;
        return result;
    }
    for (auto p = mon.begin(); p != mon.end(); ++p)
        result.push_back(GE(p->g_raw() | (p->e() * e)));
    return result;
}

int log(const Mon& mon1, const Mon& mon2)
{
    if (!mon2) {
        /* log with 0 base */
        throw ErrorIdMsg(0x88dc0b9bU, "Log with base 1");
    }
    int q = -1;

    /* Compute q */
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            ++k;
        else if (k->g_raw() < l->g_raw()) {
            q = 0;
            break;
        }
        else if (k->e() < l->e()) {
            q = 0;
            break;
        }
        else {
            int q1 = k->e_masked() / l->e_masked();
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
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            ++k;
        else if (k->g_raw() < l->g_raw())
            ++l;
        else {
            result.push_back(GE(k->g_raw() | std::min(k->e(), l->e())));
            ++k;
            ++l;
        }
    }
    return result;
}

Mon LCM(const Mon& mon1, const Mon& mon2)
{
    Mon result;
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            result.push_back(*k++);
        else if (k->g_raw() < l->g_raw())
            result.push_back(*l++);
        else {
            result.push_back(GE(k->g_raw() | std::max(k->e(), l->e())));
            ++k;
            ++l;
        }
    }
    if (k != mon1.end())
        result.insert(k, mon1.end());
    else
        result.insert(l, mon2.end());
    return result;
}

std::string Poly::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](const Mon& m) { return m.Str(); });
}

/**
 * Sort the sequence and each time remove a pair of identical elements
 */
void SortMod2(Mon1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i] = Mon::Null();
            data[++i] = Mon::Null();
        }
    ut::RemoveIf(data, [](const Mon& m) { return m.IsNull(); });
}

/**
 * Sort the sequence and each time remove a pair of identical elements
 */
void SortMod2(MMod1d& data)
{
    std::sort(data.begin(), data.end());
    for (size_t i = 0; i + 1 < data.size(); ++i)
        if (data[i] == data[i + 1]) {
            data[i] = MMod::Null();
            data[++i] = MMod::Null();
        }
    ut::RemoveIf(data, [](const MMod& m) { return m.IsNull(); });
}

void mulP(const Poly& poly, const Mon& mon, Poly& result)
{
    result.data.clear();
    result.data.reserve(poly.data.size());
    for (const Mon& m : poly.data)
        result.data.push_back(m * mon);
}

void mulP(const Poly& p1, const Poly& p2, Poly& result)
{
    result.data.clear();
    result.data.reserve(p1.data.size());
    for (size_t i = 0; i < p1.data.size(); ++i)
        for (size_t j = 0; j < p2.data.size(); ++j)
            result.data.push_back(p1.data[i] * p2.data[j]);
    SortMod2(result.data);
}

void powP(const Poly& poly, uint32_t n, Poly& result, Poly& tmp)
{
    result.data.clear();
    result.data.push_back(Mon());
    if (n == 0)
        return;
    Poly power = poly;
    while (n) {
        if (n & 1)
            result.imulP(power, tmp);
        n >>= 1;
        if (n) {
            power.frobP(tmp);
            std::swap(power, tmp);
        }
    }
}

Poly GetDiff(const Mon& mon, const Poly1d& diffs)
{
    Poly result, tmp_prod, tmp;
    for (auto k = mon.begin(); k != mon.end(); ++k) {
        if (k->e() % 2) {
            mulP(diffs[k->g()], mon / GE(k->g_raw() | 1), tmp_prod);
            result.iaddP(tmp_prod, tmp);
        }
    }
    return result;
}

Poly GetDiff(const Poly& poly, const Poly1d& diffs)
{
    Poly result, tmp_prod, tmp;
    for (const Mon& mon : poly.data) {
        for (auto k = mon.begin(); k != mon.end(); ++k) {
            if (k->e() % 2) {
                mulP(diffs[k->g()], mon / GE(k->g_raw() | 1), tmp_prod);
                result.iaddP(tmp_prod, tmp);
            }
        }
    }
    return result;
}

int1d Poly2Indices(const Poly& poly, const Mon1d& basis)
{
    int1d result;
    for (const Mon& mon : poly.data) {
        auto p = std::lower_bound(basis.begin(), basis.end(), mon);
#ifdef MYDEBUG
        if (p == basis.end() || mon < (*p)) {
            throw ErrorIdMsg(0x57f14e21U, "Index not found");
        }
#endif
        result.push_back(int(p - basis.begin()));
    }
    return result;
}

Poly Indices2Poly(const int1d& indices, const Mon1d& basis)
{
    Poly result;
    for (int i : indices)
        result.data.push_back(basis[i]);
    return result;
}

/********************************************************
 *                      Modules
 ********************************************************/

std::string MMod::Str() const
{
    std::string result = m.Str();
    if (result == "1")
        result = "v_";
    else
        result += "v_";
    std::string under = std::to_string(v);
    if (under.size() > 1)
        result += '{' + under + '}';
    else
        result += under;
    return result;
}

std::string Mod::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](const MMod& m) { return m.Str(); });
}

void mulP(const Mon& mon, const Mod& poly, Mod& result)
{
    result.data.clear();
    result.data.reserve(poly.data.size());
    for (const MMod& m : poly.data)
        result.data.emplace_back(m.m * mon, m.v);
}

void mulP(const Poly& poly, const Mod& mod, Mod& result)
{
    result.data.clear();
    result.data.reserve(poly.data.size());
    for (size_t i = 0; i < poly.data.size(); ++i)
        for (size_t j = 0; j < mod.data.size(); ++j)
            result.data.push_back(poly.data[i] * mod.data[j]);
    SortMod2(result.data);
}

int1d Mod2Indices(const Mod& x, const MMod1d& basis)
{
    int1d result;
    for (const MMod& mon : x.data) {
        auto p = std::lower_bound(basis.begin(), basis.end(), mon);
#ifndef NDEBUG
        if (p == basis.end() || mon < (*p)) {
            std::cerr << "MyException(0xc3389529U): Index not found\n";
            throw ErrorIdMsg(0xc3389529U, "Index not found");
        }
#endif
        result.push_back(int(p - basis.begin()));
    }
    return result;
}

Mod Indices2Mod(const int1d& indices, const MMod1d& basis)
{
    Mod result;
    for (int i : indices)
        result.data.push_back(basis[i]);
    return result;
}

}  // namespace alg2
