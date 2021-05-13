#include "algebras.h"
#include <iterator>

#ifdef _DEBUG
#include <iostream>
#endif

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
			k++; l++;
		}
	}
	if (k != mon1.end())
		result.insert(result.end(), k, mon1.end());
	else
		result.insert(result.end(), l, mon2.end());
	return result;
}

Mon div(const Mon& mon1, const Mon& mon2)
{
	Mon result;
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			result.push_back(*k++);
#ifdef _DEBUG
		else if (k->gen > l->gen){
			std::cout << "mon1/mon2 not divisible!\n";
			throw "1227de8e";
		}
#endif
		else if (k->exp > l->exp) {
			result.emplace_back(k->gen, k->exp - l->exp);
			k++; l++;
		}
		else if (k->exp == l->exp) {
			k++;
			l++;
		}
#ifdef _DEBUG
		else {
			std::cout << "mon1/mon2 not divisible!\n";
			throw "a9c74ef9";
		}
#endif
	}
#ifdef _DEBUG
	if (l != mon2.end()) {
		std::cout << "mon1/mon2 not divisible!\n";
		throw "6cdd66bd";
	}
	else
#endif
		result.insert(result.end(), k, mon1.end());
	return result;
}

Poly add(const Poly& poly1, const Poly& poly2) {
	Poly result;
	std::set_symmetric_difference(poly1.begin(), poly1.end(), poly2.begin(), poly2.end(),
		std::back_inserter(result));
	return result;
}

Poly mul(const Poly& poly, const Mon& mon) {
	Poly result;
	for (const Mon& m : poly)
		result.push_back(mul(m, mon));
	return result;
}

Poly mul(const Poly& poly1, const Poly& poly2) {
	Poly result;
	for (const Mon& mon2 : poly2)
		result += mul(poly1, mon2);
	return result;
}

Poly pow(const Poly& poly, int n)
{
	Poly result = { {} };
	if (n == 0)
		return result;
	Poly power = poly;
	while (n) {
		if (n & 1)
			result = mul(result, power);
		n >>= 1;
		if (n)
			power = mul(power, power);
	}
	return result;
}

bool divides(const Mon& mon1, const Mon& mon2)
{
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			return false;
		else if (k->gen > l->gen)
			l++;
		else if (k->exp > l->exp)
			return false;
		else {
			k++;
			l++;
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
			k++;
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
			k++;
			l++;
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
			k++;
		else if (k->gen > l->gen)
			l++;
		else {
			result.emplace_back(k->gen, std::min(k->exp, l->exp));
			k++; l++;
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
			k++; l++;
		}
	}
	if (k != mon1.end())
		result.insert(result.end(), k, mon1.end());
	else
		result.insert(result.end(), l, mon2.end());
	return result;
}

Poly get_diff(const Mon& mon, const Poly1d& diffs)
{
	Poly result;
	for (MonInd k = mon.begin(); k != mon.end(); ++k) {
		if (k->exp % 2)
			result += mul(diffs[k->gen], div(mon, { { k->gen, 1 } }));
	}
	return result;
}

Poly get_diff(const Poly& poly, const Poly1d& diffs)
{
	Poly result;
	for (const Mon& mon : poly) {
		for (MonInd k = mon.begin(); k != mon.end(); ++k) {
			if (k->exp % 2)
				result += result, mul(diffs[k->gen], div(mon, { { k->gen, 1 } }));
		}
	}
	return result;
}
