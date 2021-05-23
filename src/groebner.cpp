#include "groebner.h"
#include <iterator>


namespace grbn {

/**********************************************************
* Small functions
**********************************************************/

alg::Poly Reduce(alg::Poly poly, const alg::Poly1d& gb)
{
	alg::Poly result;
	auto pMon = poly.begin(); auto pEnd = poly.end();
	while (pMon != pEnd) {
		auto pGb = gb.begin();
		for (; pGb != gb.end(); ++pGb)
			if (divisible(pGb->front(), *pMon))
				break;
		if (pGb == gb.end())
			result.push_back(std::move(*pMon++));
		else {
			alg::Mon q = div(*pMon, pGb->front());
			alg::Poly rel1 = mul(*pGb, q);
			alg::Poly poly1;
			std::set_symmetric_difference(pMon, pEnd, rel1.begin(), rel1.end(),
				std::back_inserter(poly1));
			poly = std::move(poly1);
			pMon = poly.begin(); pEnd = poly.end();
		}
	}
	return result;
}

alg::Poly Reduce(alg::Poly poly, const GbWithCache& gb)
{
	alg::Poly result;
	auto pMon = poly.begin(); auto pEnd = poly.end();
	while (pMon != pEnd) {
		int index = -1;
		for (auto pGE = pMon->begin(); pGE != pMon->end() && index == -1; ++pGE)
			if (gb.lc.find(pGE->gen) != gb.lc.end())
				for (auto pIndex = gb.lc.at(pGE->gen).begin(); pIndex != gb.lc.at(pGE->gen).end() && index == -1; ++pIndex)
					if (divisible(gb.gb[*pIndex][0], *pMon))
						index = *pIndex;
		if (index == -1)
			result.push_back(std::move(*pMon++));
		else {
			alg::Mon q = div(*pMon, gb.gb[index][0]);
			alg::Poly rel1 = mul(gb.gb[index], q);
			alg::Poly poly1;
			std::set_symmetric_difference(pMon, pEnd, rel1.begin(), rel1.end(),
				std::back_inserter(poly1));
			poly = std::move(poly1);
			pMon = poly.begin(); pEnd = poly.end();
		}
	}
	return result;
}

bool gcd_nonzero(const alg::Mon& mon1, const alg::Mon& mon2)
{
	alg::MonInd k = mon1.begin(), l = mon2.begin();
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

GbBuffer GenerateBuffer(const alg::Poly1d& gb, const alg::array& gen_degs, const alg::array& gen_degs1, int t_min, int t_max)
{
	GbBuffer buffer;
	for (auto pg1 = gb.begin(); pg1 != gb.end(); ++pg1) {
		for (auto pg2 = pg1 + 1; pg2 != gb.end(); ++pg2) {
			if (gcd_nonzero(pg1->front(), pg2->front())) {
				alg::Mon lcm = LCM(pg1->front(), pg2->front());
				int deg_new_rel = get_deg(lcm, gen_degs, gen_degs1);
				if (t_min <= deg_new_rel && deg_new_rel <= t_max) {
					alg::Poly new_rel = (*pg1) * div(lcm, pg1->front()) + (*pg2) * div(lcm, pg2->front());
					buffer[deg_new_rel].push_back(std::move(new_rel));
				}
			}
		}
	}
	return buffer;
}

GbBufferV2 GenerateBufferV2(const alg::Poly1d& gb, const alg::array& gen_degs, const alg::array& gen_degs1, int t_min, int t_max)
{
	GbBufferV2 buffer;
	for (auto pg1 = gb.begin(); pg1 != gb.end(); ++pg1) {
		for (auto pg2 = pg1 + 1; pg2 != gb.end(); ++pg2) {
			alg::Mon gcd = GCD(pg1->front(), pg2->front());
			if (!gcd.empty()) {
				int deg_new_rel = get_deg(pg1->front(), gen_degs, gen_degs1) + get_deg(pg2->front(), gen_degs, gen_degs1) - get_deg(gcd, gen_degs, gen_degs1);
				if (t_min <= deg_new_rel && deg_new_rel <= t_max)
					buffer[deg_new_rel].push_back(std::make_unique<GcdBufferEle>(std::move(gcd), (int)(pg1 - gb.begin()), (int)(pg2 - gb.begin())));
			}
		}
	}
	return buffer;
}


} /* namespace grbn */