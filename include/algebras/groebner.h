/** \file groebner.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
#include <iostream> ////
#include <map>
#include <memory>

/**
 * The macro enables the multithreading for `AddRels` and functions it depends on.
 * It can be toggled off.
 */
#define GROEBNER_MULTITHREAD
#ifdef GROEBNER_MULTITHREAD
#include <future>
/**
 * The number of threads used in multithreading.
 */
constexpr int kNumThreads = 5;
#endif

// /**
//  * The macro enables the testing of some templates at compile time.
//  * It can be toggled off.
//  */
// #define GROEBNER_TEMPLATE_INSTANTIATIONS

/**
 * The namespace `grbn` encapsulates classes and functions for 
 * Groebner bases 
 */
namespace grbn {

/**
 * The auxiliary class `GbBuffer` is a set of
 * polynomials to be added to a groebner basis.
 * The polynomials are grouped by degrees.
 */
using GbBuffer = std::map<int, alg::Poly1d>;

/**
 * The class `GbLeadCache` is used to speed up
 * the incremtation of Groebner bases.
 * `lc[i]` is the array of indices of elements of `gb`
 * whose leading term's last generator id is `i`.
 */
using GbLeadCache = std::map<int, alg::array>;

/**
 * The class `GbLeadCache` consists of a Groebner basis
 * and its cache.
 */
class GbWithCache
{
public:
	GbWithCache() = default;
	GbWithCache(alg::Poly1d gb1) : gb(std::move(gb1)) { 
		for (int i = 0; i < (int)gb.size(); ++i)
			lc[gb[i][0].back().gen].push_back(i);
	}
public:
	auto begin() const { return gb.begin(); }
	auto end() const { return gb.end(); }
	auto size() const { return gb.size(); }
	void push_back(alg::Poly g) { lc[g[0].back().gen].push_back((int)gb.size()); gb.push_back(std::move(g)); }
public:
	alg::Poly1d gb;
	GbLeadCache lc;
};

/**
 * A virtual type of elements in `GbBufferV2`,
 * which is designed to save the memory.
 * @see GbBufferV2
 * 
 * The subclass `GcdBufferEle` stores the indices of
 * two elements of `gb` and cache the GCD of them and
 * it represents the new polynomials produced by those
 * two elements.
 * 
 * The subclass `PolyBufferEle` directly stores the polynomial
 * to be added to `gb`.
 */
struct BaseBufferEle {
	virtual ~BaseBufferEle() {};
	virtual alg::Poly GetPoly(const alg::Poly1d& gb) = 0;
};
/**
 * @see BaseBufferEle
 */
struct GcdBufferEle : BaseBufferEle
{
	alg::Mon gcd_; int i1_, i2_;
	GcdBufferEle(alg::Mon gcd, int i1, int i2) : gcd_(std::move(gcd)), i1_(i1), i2_(i2) {};
	alg::Poly GetPoly(const alg::Poly1d& gb) override {
		alg::Poly result = gb[i1_] * div(gb[i2_][0], gcd_) + gb[i2_] * div(gb[i1_][0], gcd_);
		alg::Mon{}.swap(gcd_); /* Deallocate */
		return result;
	}
};
/**
 * @see BaseBufferEle
 */
struct PolyBufferEle : BaseBufferEle
{
	alg::Poly p_;
	PolyBufferEle(alg::Poly p) : p_(std::move(p)) {};
	alg::Poly GetPoly(const alg::Poly1d& gb) override {
		return std::move(p_);
	}
};
/**
 * This groebner buffer type is designed to save the memory.
 */
using GbBufferV2 = std::map<int, std::vector<std::unique_ptr<BaseBufferEle>>>;

/**
 * Reduce a polynomial by a Groebner basis
 */
alg::Poly Reduce(alg::Poly poly, const alg::Poly1d& gb);
/**
 * Reduce a polynomial by a Groebner basis (with cache)
 */
alg::Poly Reduce(alg::Poly poly, const GbWithCache& gb);

/**
 * Remove elements of `cont` which are empty containers
 */
template <typename Container1d>
inline void RemoveEmptyElements(Container1d& cont)
{
	cont.erase(std::remove_if(cont.begin(), cont.end(), [](const typename Container1d::value_type& g) {return g.empty(); }), cont.end());
}

/**
 * Create the int array 1, 2, ..., n
 */
inline alg::array range(int n) {
	alg::array result;
	for (int i = 0; i < n; ++i)
		result.push_back(i);
	return result;
};

/**
 * A fast algorithm that computes
 * `poly ** n` modulo `gb`
 */
template <typename GbType>
alg::Poly pow(const alg::Poly& poly, int n, const GbType& gb)
{
	alg::Poly result = { {} };
	if (n == 0)
		return result;
	alg::Poly power = poly;
	while (n) {
		if (n & 1)
			result = Reduce(mul(result, power), gb);
		n >>= 1;
		if (n)
			power = Reduce(mul(power, power), gb);
	}
	return result;
}

/**
 * Replace the generators in `poly` with elements given in `map`
 * and evaluate modulo `gb`.
 * @param poly The polynomial to be substituted.
 * @param map `map[i]` is the polynomial that substitutes the generator of id `i`.
 * @param gb A Groebner basis.
 */
template <typename FnType, typename GbType>
alg::Poly subs(const alg::Poly& poly, FnType map, const GbType& gb)
{
	alg::Poly result;
	for (const alg::Mon& m : poly) {
		alg::Poly fm = { {} };
		for (alg::MonInd p = m.begin(); p != m.end(); ++p)
			fm = Reduce(mul(fm, pow(map(p->gen), p->exp, gb)), gb);
		result = add(result, fm);
	}
	return result;
}

bool gcd_nonzero(const alg::Mon& mon1, const alg::Mon& mon2);

/**
 * generate buffer in degree `t_min <= t <= t_max`
 */
GbBuffer GenerateBuffer(const alg::Poly1d& gb, const alg::array& gen_degs, const alg::array& gen_degs1, int t_min, int t_max);
/**
 * generate buffer in degree `t_min <= t <= t_max`
 */
GbBufferV2 GenerateBufferV2(const alg::Poly1d& gb, const alg::array& gen_degs, const alg::array& gen_degs1, int t_min, int t_max);


#ifdef GROEBNER_MULTITHREAD
/*
 * Auxiliary function used in multi-threading computing.
 */
template<typename GbType, typename InType>
void _BatchReduce(const GbType& gb, const InType& rels_in, alg::Poly1d& rels_out, int i_start)
{
	for (int i = i_start; i < (int)rels_in.size(); i += kNumThreads) {
		if constexpr (std::is_same<InType, alg::Poly1d>::value)
			rels_out[i] = Reduce(rels_in[i], gb);
		else if constexpr (std::is_same<GbType, alg::Poly1d>::value) /* InType == std::vector<std::unique_ptr<BaseBufferEle>> */
			rels_out[i] = Reduce(rels_in[i]->GetPoly(gb), gb);
		else
			rels_out[i] = Reduce(rels_in[i]->GetPoly(gb.gb), gb);
	}
}

/*
 * Auxiliary function used in multi-threading computing.
 */
template<typename GbType, typename FnType>
std::vector<std::pair<int, alg::Poly>> _BatchNewBufferElements(const GbType& gb, const alg::Poly& rel, FnType _get_deg, int deg_max, int i_start)
{
	std::vector<std::pair<int, alg::Poly>> result;
	for (auto pg = gb.begin() + i_start; pg < gb.end(); pg += kNumThreads) {
		if (gcd_nonzero(rel[0], pg->front()) && !(rel[0][0].gen < 0 && pg->front()[0].gen < 0 && rel[0][0].gen != pg->front()[0].gen)) {
			alg::Mon lcm = LCM(rel[0], pg->front());
			int deg_new_rel = _get_deg(lcm);
			if (deg_max == -1 || deg_new_rel <= deg_max) {
				alg::Poly new_rel = rel * div(lcm, rel[0]) + (*pg) * div(lcm, pg->front());
				if (!new_rel.empty())
					result.push_back(std::make_pair(deg_new_rel, std::move(new_rel)));
			}
		}
	}
	return result;
}
/*
 * Auxiliary function used in multi-threading computing.
 */
template<typename GbType, typename FnType>
std::vector<std::pair<int, GcdBufferEle>> _BatchNewBufferElementsV2(const GbType& gb, const alg::Poly& rel, FnType _get_deg, int deg, int deg_max, int i_start)
{
	std::vector<std::pair<int, GcdBufferEle>> result;
	for (auto pg = gb.begin() + i_start; pg < gb.end(); pg += kNumThreads) {
		if (gcd_nonzero(rel[0], pg->front()) && !(rel[0][0].gen < 0 && pg->front()[0].gen < 0 && rel[0][0].gen != pg->front()[0].gen)) {
			alg::Mon gcd = GCD(rel[0], pg->front());
			int deg_new_rel = deg + _get_deg(pg->front()) - _get_deg(gcd);
			if (deg_max == -1 || deg_new_rel <= deg_max)
				result.push_back(std::make_pair(deg_new_rel, GcdBufferEle(std::move(gcd), (int)(pg - gb.begin()), (int)gb.size())));
		}
	}
	return result;
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _TestBatchReduce(alg::Poly1d& gb, alg::Poly1d& rels_in, alg::Poly1d& rels_out, int i_start) { _BatchReduce(gb, rels_in, rels_out, i_start); }
inline void _TestBatchReduce(GbWithCache& gb, std::vector<std::unique_ptr<BaseBufferEle>>& rels_in, alg::Poly1d& rels_out, int i_start) { _BatchReduce(gb, rels_in, rels_out, i_start); }
inline void _TestBatchNewBufferElements(const alg::Poly1d& gb, const alg::Poly& rel, alg::FnGetDeg f, int i_start)   { _BatchNewBufferElements  (gb, rel, f, -1,     i_start); }
inline void _TestBatchNewBufferElementsV2(const alg::Poly1d& gb, const alg::Poly& rel, alg::FnGetDeg f, int i_start) { _BatchNewBufferElementsV2(gb, rel, f, -1, -1, i_start); }
#endif
#endif

/**
 * Comsume relations from `buffer` in degree `<= deg`
 * while adding new relations back to `buffer` in degree `<= deg_max`.
 * `deg=-1` or `deg_max=-1` means infinity.
 */
template <typename GbType, typename BufferType, typename FnType>
void AddRelsB(GbType& gb, BufferType& buffer, FnType _get_deg, int deg, int deg_max)
{
	auto p_buffer = buffer.begin();
	for (; p_buffer != buffer.end() && (deg == -1 || p_buffer->first <= deg); ++p_buffer) {
		/* Reduce relations from buffer in degree `p_buffer->first` */
		std::cout << "t=" << p_buffer->first << '\n'; ////
		alg::Poly1d rels;
#ifndef GROEBNER_MULTITHREAD /* Singlethreading */
		for (auto& poly : p_buffer->second) {
			alg::Poly rel = [&]() {
				if constexpr (std::is_same<BufferType, GbBuffer>::value)
					return Reduce(std::move(poly), gb);
				else if constexpr (std::is_same<GbType, alg::Poly1d>::value)
					return Reduce(std::move(poly)->GetPoly(gb), gb);
				else
					return Reduce(std::move(poly)->GetPoly(gb.gb), gb);
			}();
			for (alg::Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#else /* Multithreading */
		alg::Poly1d rels_tmp(p_buffer->second.size());
		std::vector<std::future<void>> futures;
		for (int i = 0; i < kNumThreads; ++i)
			futures.push_back(std::async(std::launch::async, [&gb, &p_buffer, &rels_tmp, i]() { return _BatchReduce(gb, p_buffer->second, rels_tmp, i); }));
		for (int i = 0; i < kNumThreads; ++i)
			futures[i].wait();
		for (auto& rel : rels_tmp) {
			for (alg::Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#endif

		/* Add these relations */
		for (auto& rel : rels) {
			if (!rel.empty()) {
#ifndef GROEBNER_MULTITHREAD /* Singlethreading */
				for (auto pg = gb.begin(); pg != gb.end(); ++pg) {
					if (gcd_nonzero(rel[0], pg->front()) && !(rel[0][0].gen < 0 && pg->front()[0].gen < 0 && rel[0][0].gen != pg->front()[0].gen)) {
						if constexpr (std::is_same<BufferType, GbBuffer>::value) {
							alg::Mon lcm = LCM(rel[0], pg->front());
							int deg_new_rel = _get_deg(lcm);
							if (deg_max == -1 || deg_new_rel <= deg_max) {
								alg::Poly new_rel = rel * div(lcm, rel[0]) + (*pg) * div(lcm, pg->front());
								if (!new_rel.empty())
									buffer[deg_new_rel].push_back(std::move(new_rel));
							}
						}
						else { /* BufferType == GbBufferV2 */
							alg::Mon gcd = GCD(rel[0], pg->front());
							int deg_new_rel = p_buffer->first + _get_deg(pg->front()) - _get_deg(gcd);
							if (deg_max == -1 || deg_new_rel <= deg_max)
								buffer[deg_new_rel].push_back(std::make_unique<GcdBufferEle>(std::move(gcd), (int)(pg - gb.begin()), (int)gb.size()));
						}
					}
				}
#else /* Multithreading */
				if constexpr (std::is_same<BufferType, GbBuffer>::value) { /* GbBuffer */
					std::vector<std::future<std::vector<std::pair<int, alg::Poly>>>> futures1;
					for (int i = 0; i < kNumThreads; ++i)
						futures1.push_back(std::async(std::launch::async, [&gb, &rel, _get_deg, deg_max, i]() { return _BatchNewBufferElements(gb, rel, _get_deg, deg_max, i); }));
					for (int i = 0; i < kNumThreads; ++i) {
						futures1[i].wait();
						for (auto& [d, p] : futures1[i].get())
							buffer[d].push_back(std::move(p));
					}
				}
				else { /* GbBufferV2 */
					std::vector<std::future<std::vector<std::pair<int, GcdBufferEle>>>> futures1;
					for (int i = 0; i < kNumThreads; ++i)
						futures1.push_back(std::async(std::launch::async, [&gb, &rel, _get_deg, &p_buffer, deg_max, i]() { return _BatchNewBufferElementsV2(gb, rel, _get_deg, p_buffer->first, deg_max, i); }));
					for (int i = 0; i < kNumThreads; ++i) {
						futures1[i].wait();
						for (auto& [d, p] : futures1[i].get())
							buffer[d].push_back(std::make_unique<GcdBufferEle>(std::move(p.gcd_), p.i1_, p.i2_));
					}
				}
#endif
				gb.push_back(std::move(rel));
			}
		}
	}
	buffer.erase(buffer.begin(), p_buffer);
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _TestAddRelsB(alg::Poly1d& gb, GbBuffer& buffer, alg::FnGetDeg f, int deg_max) { AddRelsB(gb, buffer, f, -1, deg_max); }
inline void _TestAddRelsB(GbWithCache& gb, GbBufferV2& buffer, alg::FnGetDeg f, int deg_max) { AddRelsB(gb, buffer, f, -1, deg_max); }
#endif

/**
 * @see AddRelsB(GbType&, BufferType&, FnType, int, int)
 */
template <typename GbType, typename BufferType>
void AddRelsB(GbType& gb, BufferType& buffer, const alg::array& gen_degs, int deg, int deg_max) { AddRelsB(gb, buffer, alg::FnGetDeg{ gen_degs }, deg, deg_max); }

/**
 * Convert `rels` into a groebner basis in degree `<= deg`
 * `deg_max=-1` means infinity.
 */
template <typename GbType, typename FnType, typename BufferType = GbBuffer>
void AddRels(GbType& gb, alg::Poly1d rels, FnType _get_deg, int deg_max)
{
	BufferType buffer;
	for (alg::Poly& rel : rels) {
		if (!rel.empty()) {
			int deg = _get_deg(rel[0]);
			if constexpr (std::is_same<BufferType, GbBufferV2>::value)
				buffer[deg].push_back(std::make_unique<PolyBufferEle>(std::move(rel)));
			else
				buffer[deg].push_back(std::move(rel));
		}
	}
	AddRelsB(gb, buffer, _get_deg, -1, deg_max);
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _TestAddRels(alg::Poly1d& gb, alg::Poly1d rels, alg::FnGetDeg f, int deg_max) { AddRels(gb, std::move(rels), f, deg_max); }
inline void _TestAddRels(GbWithCache& gb, alg::Poly1d rels, alg::FnGetDeg f, int deg_max) { AddRels<GbWithCache, alg::FnGetDeg, GbBufferV2>(gb, std::move(rels), f, deg_max); }
#endif

/**
 * @see AddRels(GbType&, alg::Poly1d, FnType, int)
 */
template <typename GbType>
void AddRels(GbType& gb, alg::Poly1d rels, const alg::array& gen_degs, int deg_max) { AddRels(gb, std::move(rels), alg::FnGetDeg{ gen_degs }, deg_max); }
/**
 * @see AddRels(GbType&, alg::Poly1d, FnType, int)
 */
template <typename GbType, typename FnType>
void AddRelsV2(GbType& gb, alg::Poly1d rels, FnType _get_deg, int deg_max) { AddRels<GbType, FnType, GbBufferV2>(gb, std::move(rels), _get_deg, deg_max); }
/**
 * @see AddRels(GbType&, alg::Poly1d, FnType, int)
 */
template <typename GbType>
void AddRelsV2(GbType& gb, alg::Poly1d rels, const alg::array& gen_degs, int deg_max) { AddRels<GbType, alg::FnGetDeg, GbBufferV2>(gb, std::move(rels), alg::FnGetDeg{ gen_degs }, deg_max); }


/**********************************************************
* Algorithms that use Groebner basis
**********************************************************/

/**
 * Compute the minimal generating set of `vectors` inplace.
 * 
 * Each element of vectors is considered as an element of $R^n$ where $R$ is the algebra
 * determined by the Groebner basis `gb`.
 * 
 * The degree of the basis of $R^n$ is determined by `basis_degs`.
 */
template<typename GbType>
alg::Poly2d& indecomposables(const GbType& gb, alg::Poly2d& vectors, const alg::array& gen_degs, const alg::array& basis_degs)
{
	if (vectors.empty())
		return vectors;
	GbWithCache gb1 = gb;

	/* Convert each vector v into a relation \\sum vi x_{-i-1} */
	alg::Poly1d rels;
	alg::array degs;
	for (const alg::Poly1d& v : vectors) {
		alg::Poly rel;
		for (int i = 0; i < basis_degs.size(); ++i)
			if (!v[i].empty())
				rel += v[i] * alg::Mon{ {-i - 1, 1} };
		degs.push_back(get_deg(rel[0], gen_degs, basis_degs));
		rels.push_back(std::move(rel));
	}
	alg::array indices = range((int)vectors.size());
	std::sort(indices.begin(), indices.end(), [&degs](int i, int j) {return degs[i] < degs[j]; });

	/* Add relations ordered by degree to gb1 */
	GbBuffer buffer1;
	int deg_max = degs[indices.back()];
	for (int i : indices) {
		AddRelsB(gb1, buffer1, alg::FnGetDegV2{ gen_degs, basis_degs }, degs[i], deg_max);
		alg::Poly rel = Reduce(rels[i], gb1);
		if (!rel.empty())
			buffer1[degs[i]].push_back(std::move(rel));
		else
			vectors[i].clear();
	}

	/* Keep only the indecomposables in `vectors` */
	RemoveEmptyElements(vectors);
	return vectors;
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _Test_indecomposables(const GbWithCache& gb, alg::Poly2d& vectors, const alg::array& gen_degs, const alg::array& basis_degs) { indecomposables(gb, vectors, gen_degs, basis_degs); }
#endif

/**
 * Compute the generating set of linear relations among `polys`.
 * 
 * The result is truncated by `deg<=deg_max`.
 */
template<typename GbType>
alg::Poly2d ann_seq(const GbType& gb, const alg::Poly1d& polys, const alg::array& gen_degs, int deg_max)
{
	alg::Poly2d result;
	if (polys.empty())
		return result;
	alg::Poly1d rels;
	alg::array gen_degs1;
	int N = (int)polys.size();

	/* Add relations Xi=polys[i] to gb to obtain gb1 */
	for (int i = 0; i < N; ++i) {
		alg::Poly p = polys[i];
		gen_degs1.push_back(get_deg(p, gen_degs));
		p.push_back({ {-i - 1, 1} });
		rels.push_back(std::move(p));
	}
	GbWithCache gb1 = gb;
	AddRels(gb1, rels, alg::FnGetDegV2{ gen_degs, gen_degs1 }, deg_max);

	/* Extract linear relations from gb1 */
	for (const alg::Poly& g : gb1) {
		if (g[0][0].gen < 0) {
			alg::Poly1d ann;
			ann.resize(N);
			for (const alg::Mon& m : g) {
				alg::MonInd p = m.begin();
				for (; p != m.end() && p->gen < 0; ++p);
				alg::Mon m1(m.begin(), p), m2(p, m.end());
				ann[size_t(-m1[0].gen) - 1] += Reduce(mul(subs({ div(m1, { {m1[0].gen, 1} }) }, [&polys](int i) {return polys[size_t(-i) - 1]; }, gb), m2), gb);
			}
			result.push_back(std::move(ann));
		}
	}

	/* Add commutators to linear relations */
	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; ++j) {
			if (deg_max == -1 || gen_degs1[i] + gen_degs1[j] <= deg_max) {
				alg::Poly1d result_i;
				result_i.resize(N);
				result_i[i] = polys[j];
				result_i[j] = polys[i];
				result.push_back(std::move(result_i));
			}
		}
	}

	indecomposables(gb, result, gen_degs, gen_degs1);
	return result;
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _Test_ann_seq(const GbWithCache& gb, const alg::Poly1d& polys, const alg::array& gen_degs, int deg_max) { ann_seq(gb, polys, gen_degs, deg_max); }
#endif

} /* namespace grbn */

#endif /* GROEBNER_H */
