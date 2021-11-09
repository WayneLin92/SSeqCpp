#include "groebner.h"
#include "myio.h"
#include <iterator>


namespace alg {

void detail::mul(const Mon& mon1, const Mon& mon2, Mon& mon_out)
{
    mon_out.erase(mon_out.begin(), mon_out.end()); /* No deallocation here */
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            mon_out.push_back(*k++);
        else if (k->gen > l->gen)
            mon_out.push_back(*l++);
        else {
            mon_out.emplace_back(k->gen, k->exp + l->exp);
            k++;
            l++;
        }
    }
    if (k != mon1.end())
        mon_out.insert(mon_out.end(), k, mon1.end());
    else
        mon_out.insert(mon_out.end(), l, mon2.end());
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

bool detail::GCD(const Mon& mon1, const Mon& mon2, Mon& mon_out, const array& gen_degs, int& deg, int d_min, int d_max)
{
    deg = 0;
    mon_out.erase(mon_out.begin(), mon_out.end()); /* No deallocation here */
    MonInd k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->gen < l->gen)
            k++;
        else if (k->gen > l->gen)
            l++;
        else {
            mon_out.emplace_back(k->gen, std::min(k->exp, l->exp));
            deg += gen_degs[k->gen] * mon_out.back().exp;
            if (deg > d_max)
                return false;
            k++;
            l++;
        }
    }
    return deg >= d_min;
}


} /* namespace alg */