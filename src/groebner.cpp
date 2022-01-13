#include "groebner.h"
#include "myio.h"
#include <iterator>


namespace alg {

void detail::mul(const Mon& mon1, const MonOnStack& mon2, MonOnStack& result)
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
            k++;
            l++;
        }
    }
    if (k != mon1.end())
        result.push_back(k, mon1.end());
    else
        result.push_back(l, mon2.end());
}

void detail::mul(const Mon& mon1, const Mon& mon2, MonOnStack& result)
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
            k++;
            l++;
        }
    }
    if (k != mon1.end())
        result.push_back(k, mon1.end());
    else
        result.push_back(l, mon2.end());
}

void detail::div(const Mon& mon1, const Mon& mon2, MonOnStack& result)
{
    result.clear();
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
            k++;
            l++;
        }
#ifdef NDEBUG
        else {
            k++;
            l++;
        }
#else /* DEBUG */
        else if (k->exp == l->exp) {
            k++;
            l++;
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
    result.push_back(k, mon1.end());
}

bool detail::HasGCD(const Mon& mon1, const Mon& mon2)
{
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			k++;
		else if (k->gen > l->gen)
			l++;
		else
			return true;
	}
	return false;
}

bool detail::GCD(const Mon& mon1, const Mon& mon2, MonOnStack& result, const array& gen_degs, int& deg, int d_min, int d_max)
{
    deg = 0;
    result.clear();
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            k++;
        else if (k->gen > l->gen)
            l++;
        else {
            result.emplace_back(k->gen, std::min(k->exp, l->exp));
            deg += gen_degs[k->gen] * result.back().exp;
            if (deg > d_max)
                return false;
            k++;
            l++;
        }
    }
    return deg >= d_min;
}


} /* namespace alg */