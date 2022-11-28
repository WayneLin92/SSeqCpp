#include "algebrasZ.h"
#include "myexception.h"
#include "myio.h"
#include "utility.h"
#include <iterator>

#ifndef NDEBUG
#include <iostream>
#endif

namespace algZ {

/**
 * Sort the sequence and combine the coefficients of monomials
 */
template <typename T>
void Merge(std::vector<T>& data, std::vector<T>& tmp)
{
    auto cmp = [](const T& m1, const T& m2) { return m2 < m1; };
    tmp.clear();
    std::make_heap(data.begin(), data.end(), cmp);
    while (!data.empty()) {
        if (data.front().IsUnKnown()) {
            tmp.push_back(data.front());
            break;
        }
        std::pop_heap(data.begin(), data.end(), cmp);
        if (data.size() > 1 && data.back() == data.front()) {
            data.pop_back();
            std::pop_heap(data.begin(), data.end(), cmp);
            if (data.back().Is2Torsion() || data.back().fil() >= FIL_MAX)
                data.pop_back();
            else {
                data.back().imul2();
                std::push_heap(data.begin(), data.end(), cmp);
            }
        }
        else {
            tmp.push_back(data.back());
            data.pop_back();
        }
    }
    ut::copy(tmp, data);
}

void Mon::SetFil(const AdamsDeg1d& gen_degs)
{
    fil_ = GetDegS(*this, gen_degs);
}

void Mon::SetFil(int fil)
{
    fil_ = fil;
}

std::string Mon::Str() const
{
    std::string result;
    if (c_ < 0) {
        result = "O(" + std::to_string(fil_) + ')';
        return result;
    }
    std::string str_c;
    if (c_ == 0)
        str_c = "";
    else if (c_ < 10)
        str_c = std::to_string(1 << c_);
    else
        str_c = "2^{" + std::to_string(c_) + '}';
    std::string str_m0 = myio::TplStrCont("", "", "", "", m0_.begin(), m0_.end(), [](GE p) { return p.Str(); });
    std::string str_m1 = myio::TplStrCont("", "", "", "", m1_.begin(), m1_.end(), [](GE p) { return p.Str(); });
    result = str_c + str_m1 + str_m0;
    if (result == "" || result == "-")
        result += "1";
    return result;
}

std::string Poly::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](const Mon& m) { return m.Str(); });
}

/* Determine the signed of the product in a graded-commutative algebra */
int get_sign(const Mon& mon1, const Mon& mon2, const Mon& prod)
{
    if (prod.Is2Torsion()) {
        if (prod.c() > 0)
            return 0;
        else
            return 1;
    }
    else {
        unsigned pairity = 0;
        auto& m1 = mon1.m1();
        auto& m2 = mon2.m1();
        unsigned s1 = m1.size(), s2 = m2.size();
        unsigned i1 = 0, i2 = 0;
        while (i1 < s1 && i2 < s2) {
            if (m1[i1].g() < m2[i2].g()) {
                ++i1;
            }
            else {
                pairity ^= (s1 - i1) & 1;
                ++i2;
            }
        }
        return pairity ? -1 : 1;
    }
}

Poly operator-(const Mon& mon)
{
    Poly result(mon);
    if (mon.IsUnKnown() || mon.Is2Torsion())
        return result;
    int c_max = mon.c() + (FIL_MAX - mon.fil());
    for (int c = mon.c() + 1; c <= c_max; ++c)
        result.data.emplace_back(c, mon.m0(), mon.m1(), mon.fil() + c - mon.c());
    return result;
}

Poly operator*(const Mon& mon1, const Mon& mon2)
{
    if (mon1.IsUnKnown() || mon2.IsUnKnown())
        return Mon::O(mon1.fil() + mon2.fil());
    Mon prod = mul_unsigned(mon1, mon2);
    int sign = get_sign(mon1, mon2, prod);
    if (prod.fil() > FIL_MAX || sign == 0)
        return Poly();
    else if (sign == 1)
        return prod;
    else
        return -prod;
}

void mulP(const Mon& mon, const Poly& poly, Poly& result)
{
    result.data.clear();
    Poly tmp;
    for (const Mon& m : poly.data)
        result.iaddP(mon * m, tmp);
}

void mulP(const Poly& p1, const Poly& p2, Poly& result)
{
    result.data.clear();
    result.data.reserve(p1.data.size());
    Poly tmp;
    for (size_t i = 0; i < p1.data.size(); ++i)
        for (size_t j = 0; j < p2.data.size(); ++j)
            result.iaddP(p1.data[i] * p2.data[j], tmp);
}

void append_neg(Mon1d& data, const Mon& mon, int c_null)
{
    data.push_back(mon);
    if (mon.IsUnKnown() || mon.Is2Torsion())
        return;
    int c_max = std::min(c_null - 1, mon.c() + (FIL_MAX - mon.fil()));
    for (int c = mon.c() + 1; c <= c_max; ++c)
        data.emplace_back(c, mon.m0(), mon.m1(), mon.fil() + c - mon.c());
}

void AppendMul(Poly& result, const Mon& mon1, const Mon& mon2, const int1d& gen_2tor_degs)
{
    if (mon1.IsUnKnown() || mon2.IsUnKnown()) {
        result.data.push_back(Mon::O(mon1.fil() + mon2.fil()));
        return;
    }
    Mon prod = mul_unsigned(mon1, mon2);
    int sign = get_sign(mon1, mon2, prod);
    if (prod.fil() <= FIL_MAX) {
        if (sign == 1)
            result.data.push_back(prod);
        else if (sign == -1)
            append_neg(result.data, prod, gen_2tor_degs[prod.frontg()]);
    }
}

void AppendNegMul(Poly& result, const Mon& mon1, const Mon& mon2, const int1d& gen_2tor_degs)
{
    if (mon1.IsUnKnown() || mon2.IsUnKnown()) {
        result.data.push_back(Mon::O(mon1.fil() + mon2.fil()));
        return;
    }
    Mon prod = mul_unsigned(mon1, mon2);
    int sign = get_sign(mon1, mon2, prod);
    if (prod.fil() <= FIL_MAX) {
        if (sign == -1)
            result.data.push_back(prod);
        else if (sign == 1)
            append_neg(result.data, prod, gen_2tor_degs[prod.frontg()]);
    }
}

Poly& Poly::iaddP(const Poly& rhs, Poly& tmp)
{
    data.insert(data.end(), rhs.data.begin(), rhs.data.end());
    Merge(data, tmp.data);
    return *this;
}

void Poly::iaddmulP(const Mon& mon, const Poly& poly, Poly& tmp, const int1d& gen_2tor_degs)
{
    for (const Mon& m : poly.data)
        AppendMul(*this, mon, m, gen_2tor_degs);
    Merge(data, tmp.data);
}

void Poly::isubmulP(const Mon& mon, const Poly& poly, Poly& tmp, const int1d& gen_2tor_degs)
{
    for (const Mon& m : poly.data)
        AppendNegMul(*this, mon, m, gen_2tor_degs);
    Merge(data, tmp.data);
}

Poly Poly::negP(Poly& tmp, const int1d& gen_2tor_degs) const
{
    Poly result;
    for (const Mon& m : data)
        append_neg(result.data, m, gen_2tor_degs[m.frontg()]);
    Merge(result.data, tmp.data);
    return result;
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
        if (n)
            power = power * power;
    }
}

int1d Poly2Indices(const Poly& poly, const Mon1d& basis)
{
    int1d result;
    for (const Mon& mon : poly.data) {
        auto p = std::lower_bound(basis.begin(), basis.end(), mon);
#ifndef NDEBUG
        if (p == basis.end() || mon < (*p)) {
            throw MyException(0x57f14e21U, "Index not found");
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

Mod operator-(const MMod& mon)
{
    return mon.IsUnKnown() ? Mod(mon) : Mod(-mon.m, mon.v, 0);
}

Mod operator*(const Mon& m, const MMod& x)
{
    Mod result;
    const Poly p = m * x.m;
    for (auto& mm : p.data)
        result.data.push_back(MMod(mm, x.v, 0));
    return result;
}

void append_neg(MMod1d& data, const MMod& mon, int c_null)
{
    data.push_back(mon);
    if (mon.IsUnKnown() || mon.Is2Torsion())
        return;
    int c_max = std::min(c_null - 1, mon.m.c() + (FIL_MAX - mon.fil()));
    for (int c = mon.m.c() + 1; c <= c_max; ++c)
        data.emplace_back(c, mon.m.m0(), mon.m.m1(), mon.fil() + c - mon.m.c(), mon.v);
}

void AppendMul(Mod& result, const Mon& m, const MMod& x, const int1d& gen_2tor_degs)
{
    if (m.IsUnKnown() || x.IsUnKnown()) {
        result.data.push_back(MMod::O(m.fil() + x.fil()));
        return;
    }
    Mon prod = mul_unsigned(m, x.m);
    int sign = get_sign(m, x.m, prod);
    if (prod.fil() <= FIL_MAX) {
        if (sign == 1)
            result.data.push_back(MMod(prod, x.v, 0));
        else if (sign == -1)
            append_neg(result.data, MMod(prod, x.v, 0), gen_2tor_degs[prod.frontg()]);
    }
}

void AppendNegMul(Mod& result, const Mon& m, const MMod& x, const int1d& gen_2tor_degs)
{
    if (m.IsUnKnown() || x.IsUnKnown()) {
        result.data.push_back(MMod::O(m.fil() + x.fil()));
        return;
    }
    Mon prod = mul_unsigned(m, x.m);
    int sign = get_sign(m, x.m, prod);
    if (prod.fil() <= FIL_MAX) {
        if (sign == -1)
            result.data.push_back(MMod(prod, x.v, 0));
        else if (sign == 1)
            append_neg(result.data, MMod(prod, x.v, 0), gen_2tor_degs[prod.frontg()]);
    }
}

std::string MMod::Str() const
{
    std::string result = m.Str();
    if (m.IsUnKnown())
        return result;
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

Mod& Mod::iaddP(const Mod& rhs, Mod& tmp)
{
    data.insert(data.end(), rhs.data.begin(), rhs.data.end());
    Merge(data, tmp.data);
    return *this;
}

Mod Mod::operator-() const
{
    Mod result, tmp;
    for (auto& m : data)
        result.iaddP(-m, tmp);
    return result;
}

void Mod::iaddmulP(const Mon& m, const Mod& x, Mod& tmp, const int1d& gen_2tor_degs)
{
    for (auto& xm : x.data)
        AppendMul(*this, m, xm, gen_2tor_degs);
    Merge(data, tmp.data);
}

void Mod::isubmulP(const Mon& m, const Mod& x, Mod& tmp, const int1d& gen_2tor_degs)
{
    for (auto& xm : x.data)
        AppendNegMul(*this, m, xm, gen_2tor_degs);
    Merge(data, tmp.data);
}

Mod Mod::negP(Mod& tmp, const int1d& gen_2tor_degs) const
{
    Mod result;
    for (const MMod& m : data)
        append_neg(result.data, m, gen_2tor_degs[m.m.frontg()]);
    Merge(result.data, tmp.data);
    return result;
}

void Mod::ReduceSizeByChangingSignP(Mod& tmp1, Mod& tmp2)
{
    tmp1.data.clear();
    for (auto& m : data)
        tmp1.iaddP(-m, tmp2);
    if (tmp1.data.size() < data.size())
        ut::copy(tmp1.data, data);
}

std::string Mod::Str() const
{
    return myio::TplStrCont("", "+", "", "0", data.begin(), data.end(), [](const MMod& m) { return m.Str(); });
}

Poly subsMod(const Mod& x, const std::vector<Poly>& map, const AdamsDeg1d& v_degs)
{
    Poly result, tmp;
    for (const MMod& m : x.data) {
        if (m.IsUnKnown())
            result += Mon::O(m.fil());
        else {
            Mon m1 = m.m;
            m1.SetFil(m1.fil() - v_degs[m.v].s);
            result.iaddP(Poly(m1) * map[m.v], tmp);
        }
    }
    return result;
}

int1d Mod2Indices(const Mod& x, const MMod1d& basis)
{
    int1d result;
    for (const MMod& mon : x.data) {
        auto p = std::lower_bound(basis.begin(), basis.end(), mon);
#ifndef NDEBUG
        if (p == basis.end() || mon < (*p)) {
            std::cerr << "MyException(0x62885502U): Index not found\n";
            throw MyException(0x62885502U, "Index not found");
        }
#endif
        result.push_back(uint32_t(p - basis.begin()));
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

}  // namespace algZ
