/*
** Functions for Monomials and Polynomials
** A monomial m is an array of i1, e1, i2, e2, ...
** Monomials are ordered lexicographically by exponents after filling in zero exponents.
** A polynomials is an increasing sequence of monomials.
*/

#ifndef ALGEBRAS_H
#define ALGEBRAS_H

#include <queue>

/********** STRUCTS AND CLASSES **********/
using array = std::vector<int>;
using array2d = std::vector<array>;
using array3d = std::vector<array2d>;
using array4d = std::vector<array3d>;

struct GenPow {
	int gen, exp;
	GenPow(int g, int e) : gen(g), exp(e) {}
	bool operator<(const GenPow& rhs) const { return gen > rhs.gen || (gen == rhs.gen && exp < rhs.exp); }
	bool operator==(const GenPow& rhs) const { return gen == rhs.gen && exp == rhs.exp; }
};
using Mon = std::vector<GenPow>;
using Poly = std::vector<Mon>;
using Poly1d = std::vector<Poly>;
using Poly2d = std::vector<Poly1d>;
using Mon1d = std::vector<Mon>;
using Mon2d = std::vector<Mon1d>;
using MonInd = Mon::const_iterator;

struct Deg
{
	int s, t, v;
	Deg() = default;
	Deg(int s_, int t_, int v_) : s(s_), t(t_), v(v_) {};
	bool operator<(const Deg& rhs) const {
		if (t < rhs.t)
			return true;
		else if (t == rhs.t) {
			if (s < rhs.s)
				return true;
			else if (s == rhs.s)
				if (v < rhs.v)
					return true;
		}
		return false;
	};
	Deg operator+(const Deg& rhs) const {
		return Deg({ s + rhs.s, t + rhs.t, v + rhs.v });
	};
	Deg operator-(const Deg& rhs) const {
		return Deg({ s - rhs.s, t - rhs.t, v - rhs.v });
	};
	Deg& operator+=(const Deg& rhs) {
		s += rhs.s; t += rhs.t; v += rhs.v;
		return *this;
	};
	Deg& operator-=(const Deg& rhs) {
		s -= rhs.s; t -= rhs.t; v -= rhs.v;
		return *this;
	};
	bool operator==(const Deg& rhs) const {
		return v == rhs.v && s == rhs.s && t == rhs.t;
	};
	bool operator!=(const Deg& rhs) const {
		return t != rhs.t || s != rhs.s || v != rhs.v;
	};
	Deg operator*(int rhs) const {
		return Deg({ s * rhs, t * rhs, v * rhs });
	};
};

/********** FUNCTIONS **********/
Mon mul(const Mon& mon1, const Mon& mon2);
Mon div(const Mon& mon1, const Mon& mon2);
Poly add(const Poly& poly1, const Poly& poly2);
Poly mul(const Poly& poly, const Mon& mon);
inline Poly mul(const Mon& mon, const Poly& poly) { return mul(poly, mon); }
Poly mul(const Poly& poly1, const Poly& poly2);
Mon pow(const Mon& m, int e);
Poly pow(const Poly& poly, int n);
bool divides(const Mon& m1, const Mon& m2);
/* return the largest integer e where m1 = m2^e * r */
int log(const Mon& m1, const Mon& m2);

Mon GCD(const Mon& m1, const Mon& m2);
Mon LCM(const Mon& m1, const Mon& m2);

/* Operator overloadings */
inline Poly operator+(const Poly& lhs, const Poly& rhs) { return add(lhs, rhs); }
inline Poly& operator+=(Poly& lhs, const Poly& rhs) { lhs = add(lhs, rhs); return lhs; }
inline Poly operator*(const Poly& lhs, const Poly& rhs) { return mul(lhs, rhs); }
inline Poly operator*(const Poly& lhs, const Mon& rhs) { return mul(lhs, rhs); }
inline Poly operator*(const Mon& lhs, const Poly& rhs) { return mul(lhs, rhs); }
inline Mon operator/(const Mon& lhs, const Mon& rhs) { return div(lhs, rhs); }

inline int get_deg(const Mon& mon) { int result = 0; for (MonInd p = mon.begin(); p != mon.end(); ++p) result += p->exp; return result; };
inline int get_deg(const Mon& mon, const array& gen_degs) { int result = 0; for (MonInd p = mon.begin(); p != mon.end(); ++p) result += gen_degs[p->gen] * p->exp; return result; };
inline int get_deg(const Mon& mon, const array& gen_degs, const array& gen_degs1) {
	int result = 0;
	for (MonInd p = mon.begin(); p != mon.end(); ++p)
		result += (p->gen >= 0 ? gen_degs[p->gen] : gen_degs1[size_t(-p->gen) - 1]) * p->exp;
	return result;
}
inline Deg get_deg(const Mon& mon, const std::vector<Deg>& gen_degs) { Deg result({ 0, 0, 0 }); for (MonInd p = mon.begin(); p != mon.end(); ++p) result += gen_degs[p->gen] * p->exp; return result; };
inline int get_deg_t(const Mon& mon, const std::vector<Deg>& gen_degs) { int result = 0; for (MonInd p = mon.begin(); p != mon.end(); ++p) result += gen_degs[p->gen].t * p->exp; return result; };

inline int get_deg(const Poly& poly) { return poly.empty() ? -10000 : get_deg(poly[0]); };
inline int get_deg(const Poly& poly, const array& gen_degs) { return poly.empty() ? -10000 : get_deg(poly[0], gen_degs); }
inline int get_deg(const Poly& poly, const array& gen_degs, const array& gen_degs1) { return poly.empty() ? -10000 : get_deg(poly[0], gen_degs, gen_degs1); }
inline Deg get_deg(const Poly& poly, const std::vector<Deg>& gen_degs) { return poly.empty() ? Deg({ -10000, -10000, -10000 }) : get_deg(poly[0], gen_degs); }
inline int get_deg_t(const Poly& poly, const std::vector<Deg>& gen_degs) { return poly.empty() ? -10000 : get_deg_t(poly[0], gen_degs); }

/* Functors which compute the degree of a monomial */
struct FnGetDeg {
	const array& gen_degs;
	int operator()(const Mon& mon) const { return get_deg(mon, gen_degs); }
};
struct FnGetDegV2 {
	const array& gen_degs;
	const array& gen_degs1;
	int operator()(const Mon& mon) const { return get_deg(mon, gen_degs, gen_degs1); }
};

/* Compute the differential */
Poly get_diff(const Mon& mon, const Poly1d& diffs);
Poly get_diff(const Poly& poly, const Poly1d& diffs);

#endif /* ALGEBRAS_H */