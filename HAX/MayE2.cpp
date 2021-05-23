#include "MayE2.h"
#include "algebras/linalg.h"
#include "algebras/benchmark.h"
#include "algebras/myexception.h"

#define MULTITHREAD
#define BENCHMARK

#ifdef MULTITHREAD
#include <future>
#endif

alg::array ReduceToIndices(alg::Poly poly, const grbn::GbWithCache& gb, const alg::Mon1d& basis)
{
	alg::array result;
	auto pbegin = poly.begin(); auto pend = poly.end();
	while (pbegin != pend) {
		auto first = std::lower_bound(basis.begin(), basis.end(), *pbegin);
		if (first != basis.end() && !(*pbegin < *first)) {
			result.push_back(int(first - basis.begin()));
			++pbegin;
		}
		else {
			auto pGb = gb.begin();
			for (; pGb != gb.end(); ++pGb)
				if (divisible(pGb->front(), *pbegin))
					break;
			alg::Mon q = div(*pbegin, pGb->front());
			alg::Poly rel1 = mul(*pGb, q);
			alg::Poly poly1;
			std::set_symmetric_difference(pbegin, pend, rel1.begin(), rel1.end(),
				std::back_inserter(poly1));
			poly = std::move(poly1);
			pbegin = poly.begin(); pend = poly.end();
		}
	}
	return result;
}

/* Add rels from gb in degree `t` to gb1 */
void AddRelsFromGb(const grbn::GbWithCache& gb, const alg::array& gen_degs, grbn::GbWithCache& gb1, grbn::GbBufferV2& buffer1, int t, int t_max)
{
	/* Add relations from gb to gb1 */
	auto p1 = std::lower_bound(gb.begin(), gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t; });
	auto p2 = std::lower_bound(p1, gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t + 1; });
	for (auto pg = p1; pg != p2; ++pg)
		buffer1[t].push_back(std::make_unique<grbn::PolyBufferEle>(*pg));
	grbn::AddRelsB(gb1, buffer1, gen_degs, t, t_max);
}

/* Add reindexed rels from gb in degree `t` to gb1 */
void AddRelsFromGb(const grbn::GbWithCache& gb, const alg::array& gen_degs, grbn::GbWithCache& gb1, grbn::GbBufferV2& buffer1, const alg::array& gen_degs1, const alg::array& map_gen_id, int t, int t_max)
{
	/* Add relations from gb to gb1 */
	auto p1 = std::lower_bound(gb.begin(), gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t; });
	auto p2 = std::lower_bound(p1, gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t + 1; });
	for (auto pg = p1; pg != p2; ++pg)
		buffer1[t].push_back(std::make_unique<grbn::PolyBufferEle>(reindex(*pg, map_gen_id)));
	grbn::AddRelsB(gb1, buffer1, gen_degs1, t, t_max);
}

/* compute ann * polys = 0 in degree t */
alg::Poly2d ExtendAnn(const grbn::GbWithCache& gb, const alg::array& gen_degs, grbn::GbWithCache& gb1, grbn::GbBufferV2& buffer1, const alg::Poly1d& polys, const alg::array& deg_polys, int t, int t_max) // Copy polys by value
{
	/* Add relations from gb to gb1 */
	auto p1 = std::lower_bound(gb.begin(), gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t; });
	auto p2 = std::lower_bound(p1, gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t + 1; });
	for (auto pg = p1; pg != p2; ++pg)
		buffer1[t].push_back(std::make_unique<grbn::PolyBufferEle>(*pg));

	/* Add relations Xi=polys[i] to gb1 */
	int N = (int)polys.size();
	for (int i = N - 1; i >= 0; --i) {
		if (deg_polys[i] != t)
			break;
		alg::Poly p = polys[i];
		p.push_back({ {-i - 1, 1} });
		buffer1[t].push_back(std::make_unique<grbn::PolyBufferEle>(std::move(p)));
	}
	grbn::AddRelsB(gb1, buffer1, alg::FnGetDegV2{ gen_degs, deg_polys }, t, t_max);

	alg::Poly2d result;
	if (polys.empty())
		return result;

	/* Extract linear relations from gb1 */
	for (auto pg = gb1.gb.rbegin(); pg != gb1.gb.rend(); ++pg) {
		if (get_deg(pg->front(), gen_degs, deg_polys) != t)
			break;
		if (pg->front()[0].gen < 0) {
			alg::Poly1d ann;
			ann.resize(N);
			for (const alg::Mon& m : *pg) {
				alg::MonInd p = m.begin();
				for (; p != m.end() && p->gen < 0; ++p);
				alg::Mon m1(m.begin(), p), m2(p, m.end());
				ann[size_t(-m1[0].gen) - 1] += grbn::Reduce(mul(grbn::subs({ div(m1, { {m1[0].gen, 1} }) }, [&polys](int i) {return polys[size_t(-i) - 1]; }, gb), m2), gb);
			}
			result.push_back(std::move(ann));
		}
	}

	/* Add commutators to linear relations */
	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; ++j) {
			if (deg_polys[i] + deg_polys[j] == t) {
				alg::Poly1d result_i;
				result_i.resize(N);
				result_i[i] = polys[j];
				result_i[j] = polys[i];
				result.push_back(std::move(result_i));
			}
		}
	}

	return result;
}

/* Remove decomposables */
void Indecomposables(const grbn::GbWithCache& gb, const alg::array& gen_degs, grbn::GbWithCache& gb1, grbn::GbBufferV2& buffer1, alg::Poly2d& vectors, const alg::array& basis_degs, int t, int t_max)
{
	/* Add relations from gb to gb1 */
	auto p1 = std::lower_bound(gb.begin(), gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t; });
	auto p2 = std::lower_bound(p1, gb.end(), t, [&gen_degs](const alg::Poly& g, int t) {return get_deg(g, gen_degs) < t + 1; });
	for (auto pg = p1; pg != p2; ++pg)
		buffer1[t].push_back(std::make_unique<grbn::PolyBufferEle>(*pg));
	grbn::AddRelsB(gb1, buffer1, alg::FnGetDegV2{ gen_degs, basis_degs }, t, t_max);

	/* Convert each vector v to a relation \\sum vi x_{-i-1} */
	for (size_t j = 0; j < vectors.size(); ++j) {
		alg::Poly rel;
		for (int i = 0; i < basis_degs.size(); ++i)
			if (!vectors[j][i].empty())
				rel += vectors[j][i] * alg::Mon{ {-i - 1, 1} };
		rel = grbn::Reduce(std::move(rel), gb1);
		if (!rel.empty())
			buffer1[t].push_back(std::make_unique<grbn::PolyBufferEle>(std::move(rel)));
		else
			vectors[j].clear();
	}
	grbn::AddRelsB(gb1, buffer1, alg::FnGetDegV2{ gen_degs, basis_degs }, t, t_max);

	/* Keep only the indecomposables in `vectors` */
	grbn::RemoveEmptyElements(vectors);
}

std::map<alg::Deg, alg::Mon1d> ExtendBasis(const std::vector<alg::Deg>& gen_degs, const alg::Mon2d& leadings, std::map<alg::Deg, alg::Mon1d>& basis, int t)
{
	std::map<alg::Deg, alg::Mon1d> basis_t;
	if (t == 0) {
		basis_t[alg::Deg{ 0, 0, 0 }].push_back({});
		return basis_t;
	}
	for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 0; --gen_id) {
		int t1 = t - gen_degs[gen_id].t;
		if (t1 >= 0) {
			auto p1_basis = basis.lower_bound(alg::Deg{ 0, t1, 0 });
			auto p2_basis = basis.lower_bound(alg::Deg{ 0, t1 + 1, 0 });
			for (auto p_basis = p1_basis; p_basis != p2_basis; ++p_basis) {
				for (auto p_m = p_basis->second.begin(); p_m != p_basis->second.end(); ++p_m) {
					if (p_m->empty() || gen_id <= p_m->front().gen) {
						alg::Mon mon(mul(*p_m, { {gen_id, 1} }));
						if ((size_t)gen_id >= leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
							[&mon](const alg::Mon& _m) { return divisible(_m, mon); }))
							basis_t[p_basis->first + gen_degs[gen_id]].push_back(std::move(mon));
					}
				}
			}
		}
	}
	return basis_t;
}

/* Compute relations in HB where B=A[X] in degrees (s, t, v) with fixed t and v2s=v+2*s */
alg::Poly1d FindRels(std::map<alg::Deg, alg::Mon1d>& basis_HB, const std::map<alg::Deg, DgaBasis1>& basis_A, const std::map<alg::Deg, DgaBasis1>& basis_X, const grbn::GbWithCache& gb_A, const alg::Poly1d& gen_reprs_B,
	int t, int v2s, const alg::array& v_s)
{
	alg::Mon1d basis_B_s; /* in deg s */
	alg::Poly1d diffs_B_s; /* in deg s + 1 */
	alg::Mon1d basis_B_sm1; /* in deg s - 1 */
	alg::Poly1d diffs_B_sm1; /* in deg s */
	alg::Poly1d result;
	for (size_t i = 0; i < v_s.size(); ++i) {
		int s = v_s[i];
		alg::Deg d = { s, t, v2s - 2 * s };
		if (i > 0 && v_s[i - 1] == s - 1) {
			std::swap(basis_B_sm1, basis_B_s); basis_B_s.clear();
			std::swap(diffs_B_sm1, diffs_B_s); diffs_B_s.clear();
		}
		else {
			get_basis_with_diff_B(basis_A, basis_X, basis_B_sm1, diffs_B_sm1, d - alg::Deg{ 1, 0, -2 });
		}
		if (i < v_s.size() - 1 && v_s[i + 1] == s + 1)
			get_basis_with_diff_B(basis_A, basis_X, basis_B_s, diffs_B_s, d);
		else
			get_basis_B(basis_A, basis_X, basis_B_s, d);
		std::sort(basis_B_s.begin(), basis_B_s.end());

		alg::array2d map_diff;
		for (alg::Poly& diff : diffs_B_sm1)
			map_diff.push_back(ReduceToIndices(diff, gb_A, basis_B_s)); // Profiler: t=100, 5650
		alg::array2d image_diff = lina::GetSpace(map_diff);

		alg::array2d map_repr;
		for (const alg::Mon& mon : basis_HB[d]) {
			alg::array repr = lina::Residue(image_diff, Poly_to_indices(get_image({ mon }, gen_reprs_B, gb_A), basis_B_s));
			map_repr.push_back(std::move(repr));
		}
		alg::array2d image_repr, kernel_repr, g_repr;
		lina::SetLinearMap(map_repr, image_repr, kernel_repr, g_repr);

		for (const alg::array& rel_indices : kernel_repr)
			result.push_back(indices_to_Poly(rel_indices, basis_HB[d]));
		for (const alg::array& rel_indices : kernel_repr)
			basis_HB[d][rel_indices[0]].clear();

		grbn::RemoveEmptyElements(basis_HB[d]);
	}
	return result;
}

void generate_HB(const Database& db, int t_max, int t_max_compute=-1, bool drop_existing=false)
{
	/*# Load data */
	size_t index_x = (size_t)db.get_int("SELECT COUNT(*) FROM A_generators;") - 1; /* the gen_id of Xn is index_x + n */
	size_t n = (size_t)db.get_int("SELECT COUNT(*) FROM B_generators;") - index_x - 1;
#ifdef BENCHMARK
	n = 6;
#endif
	int t_min;
	if (drop_existing)
		t_min = 0;
	else {
		try { t_min = db.get_int("SELECT MAX(t) FROM HA1_basis;") + 1; } ////
		catch (MyException&) { t_min = 0; }
	}
	std::vector<alg::Deg> gen_degs_B = db.load_gen_degs("B_generators");
	alg::Poly1d gen_diffs_B = db.load_gen_diffs("B_generators");
	grbn::GbWithCache gb_A0 = db.load_gb("A_relations", t_max);
	const std::map<alg::Deg, DgaBasis1> basis_A0 = load_dga_basis(db, "A_basis", 2, t_max);
	std::vector<std::map<alg::Deg, DgaBasis1>> basis_X;

	/* Data that needs to be updated */
	std::vector<std::vector<alg::Deg>> gen_degs_HA;
	std::vector<alg::array> gen_degs_t_HA; /* derived from gen_degs_HA */
	std::vector<alg::Poly1d> gen_reprs_HA;
	std::vector<grbn::GbWithCache> gb_HA;
	std::vector<grbn::GbLeadCache> lc_HA(n + 1);
	std::vector<alg::Mon2d> leadings_HA; /* derived from gb_HA */
	std::vector<grbn::GbBufferV2> buffer_HA(n + 1); /* derived from gb_HA */
	std::vector<std::map<alg::Deg, alg::Mon1d>> basis_HA;
	std::vector<alg::Poly1d> y;
	std::vector<alg::array> t_y;
	std::vector<alg::array> map_gen_id; /* new gen_id of HA[i] in HA[i+1] */

	std::vector<grbn::GbWithCache> gb_HA_ann_c;
	std::vector<grbn::GbWithCache> gb_HA_ind_y;
	std::vector<grbn::GbWithCache> gb_HA_ann_y;
	std::vector<grbn::GbWithCache> gb_HA_ind_a;
	std::vector<grbn::GbBufferV2> buffer_HA_ann_c(n); /* derived from gb_HA_ann_c */
	std::vector<grbn::GbBufferV2> buffer_HA_ind_y(n); /* derived from gb_HA_ind_y */
	std::vector<grbn::GbBufferV2> buffer_HA_ann_y(n); /* derived from gb_HA_ann_y */
	std::vector<grbn::GbBufferV2> buffer_HA_ind_a(n); /* derived from gb_HA_ind_a */
	for (int i = 0; i <= n; ++i) {
		std::string table_prefix = i == 0 ? "HA" : "HA" + std::to_string(i);

		/* Load gen_degs_HA */
		try { db.execute_cmd("CREATE TABLE " + table_prefix + "_generators (gen_id INTEGER PRIMARY KEY, gen_name TEXT UNIQUE, gen_diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
		catch (MyException&) {}
		if (i > 0 && drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_generators;");
		gen_degs_HA.push_back(db.load_gen_degs(table_prefix + "_generators"));

		/* Initialize gen_degs_t_HA */
		gen_degs_t_HA.push_back({});
		for (auto p = gen_degs_HA.back().begin(); p < gen_degs_HA.back().end(); ++p)
			gen_degs_t_HA.back().push_back(p->t);

		/* Load gen_reprs_HA */
		gen_reprs_HA.push_back(db.load_gen_reprs(table_prefix + "_generators"));

		/* Load gb_HA */
		try { db.execute_cmd("CREATE TABLE " + table_prefix + "_relations (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
		catch (MyException&) {}
		if (i > 0 && drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_relations;");
		gb_HA.push_back(db.load_gb(table_prefix + "_relations", t_max));

		/* Initialize leadings_HA */
		leadings_HA.push_back({});

		/* Generate buffer_HA */
		buffer_HA[i] = i == 0 ? grbn::GbBufferV2{} : grbn::GenerateBufferV2(gb_HA[i].gb, gen_degs_t_HA[i], {}, t_min, t_max);

		/* Load basis_HA */
		try { db.execute_cmd("CREATE TABLE " + table_prefix + "_basis  (mon_id INTEGER PRIMARY KEY, mon TEXT NOT NULL UNIQUE, diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
		catch (MyException&) {}
		if (drop_existing && i > 0) db.execute_cmd("DELETE FROM " + table_prefix + "_basis;");
		basis_HA.push_back(db.load_basis(table_prefix + "_basis", t_max));

		/* Load basis_X */
		basis_X.push_back(load_basis_X(db, "X" + std::to_string(i) + "_basis", t_max, 2));

		if (i < n){
			/* Load y and t_y */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_y  (y TEXT, t SMALLINT);"); }
			catch (MyException&) {}
			y.push_back({}); t_y.push_back({});
			if (drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_y;");
			load_y(db, table_prefix + "_y", y.back(), t_y.back());

			/* Load map_gen_id */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_map_gen_id  (gen_id INT);"); }
			catch (MyException&) {}
			if (drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_map_gen_id;");
			map_gen_id.push_back(db.get_ints(table_prefix + "_map_gen_id", "gen_id"));

			/* Load gb_HA_ann_c */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ann_c_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (MyException&) {}
			if (drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_ann_c_gb;");
			gb_HA_ann_c.push_back(db.load_gb(table_prefix + "_ann_c_gb", t_max));

			/* Load gb_HA_ann_y */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ann_y_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (MyException&) {}
			if (drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_ann_y_gb;");
			gb_HA_ann_y.push_back(db.load_gb(table_prefix + "_ann_y_gb", t_max));

			/* Load gb_HA_ind_y */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ind_y_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (MyException&) {}
			if (drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_ind_y_gb;");
			gb_HA_ind_y.push_back(db.load_gb(table_prefix + "_ind_y_gb", t_max));

			/* Load gb_HA_ind_a */
			try { db.execute_cmd("CREATE TABLE " + table_prefix + "_ind_a_gb (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);"); }
			catch (MyException&) {}
			if (drop_existing) db.execute_cmd("DELETE FROM " + table_prefix + "_ind_a_gb;");
			gb_HA_ind_a.push_back(db.load_gb(table_prefix + "_ind_a_gb", t_max));

			int t_x = gen_degs_B[index_x + i + 1].t; /* x = x_{i + 1} */
			buffer_HA_ann_c[i] = grbn::GenerateBufferV2(gb_HA_ann_c[i].gb, gen_degs_t_HA[i], { gen_degs_B[index_x + i + 1].t }, t_min, t_max);
			buffer_HA_ind_y[i] = grbn::GenerateBufferV2(gb_HA_ind_y[i].gb, gen_degs_t_HA[i], { gen_degs_B[index_x + i + 1].t }, t_min, t_max);
			buffer_HA_ann_y[i] = grbn::GenerateBufferV2(gb_HA_ann_y[i].gb, gen_degs_t_HA[i], t_y[i], t_min - t_x, t_max - t_x);
			buffer_HA_ind_a[i] = grbn::GenerateBufferV2(gb_HA_ind_a[i].gb, gen_degs_t_HA[i], t_y[i], t_min - t_x, t_max - t_x);
		}
	}

	/*# Compute */
	alg::Poly1d c(n); /* dx[i+1] = c[i] */
	for (size_t i = 0; i < (size_t)n; ++i) {
		int t_x = gen_degs_B[index_x + i + 1].t; /* x = x_{i + 1} */
		if (t_min > t_x)
			c[i] = proj(gen_diffs_B[index_x + i + 1], gen_degs_B, gen_diffs_B, gb_A0, basis_A0, basis_X[i], gen_reprs_HA[i], basis_HA[i]);
	}
	int t_m = t_max_compute == -1 ? t_max : t_max_compute;
	for (int t = t_min; t <= t_m; ++t) {
		std::cout << "\033[0;32m" << "t=" << t << "\033[0m\n";
		db.begin_transaction();
		for (size_t i = 0; i < (size_t)n; ++i) {
			/*## Compute ann(c) and ann(y) */
			/* Add x_{-1}=c to gb_HA_ann_c in order to compute ann_HA(c) */

			std::string table_prefix = i == 0 ? "HA" : "HA" + std::to_string(i);
			std::string table_H_prefix = "HA" + std::to_string(i + 1);
			alg::Deg deg_x = gen_degs_B[index_x + i + 1];
			int t_x = deg_x.t; /* x = x_{i + 1} */
			if (t == t_x) {
				c[i] = proj(gen_diffs_B[index_x + i + 1], gen_degs_B, gen_diffs_B, gb_A0, basis_A0, basis_X[i], gen_reprs_HA[i], basis_HA[i]);
				buffer_HA[i + 1][t].push_back(std::make_unique<grbn::PolyBufferEle>(c[i]));
				grbn::AddRelsB(gb_HA[i + 1], buffer_HA[i + 1], gen_degs_t_HA[i + 1], t, t_max); /* Add relation dx=0 */
			}

			alg::Poly2d y_new; /* annilators of c in gb_HA[i] in degree t - t_x */
			std::map<alg::Deg, alg::Mon1d> basis_t;

			int gen_degs_HA_t_start = (int)gen_degs_HA[i + 1].size(); /* starting index of generators of HA[i+1] in t */
			int gb_HA_t_start = (int)gb_HA[i + 1].size(); /* starting index of relations in t */
			int map_gen_id_t_start = (int)map_gen_id[i].size(); /* starting index of gen_id in HA[i] in t */
			auto p1_degs = std::lower_bound(gen_degs_HA[i].begin(), gen_degs_HA[i].end(), t, [](alg::Deg d, int t_) {return d.t < t_; });
			auto p2_degs = std::lower_bound(p1_degs, gen_degs_HA[i].end(), t + 1, [](alg::Deg d, int t_) {return d.t < t_; });
			for (auto p_degs = p1_degs; p_degs != p2_degs; ++p_degs) { /* Add generators from HA[i] to HA[i + 1] */
				map_gen_id[i].push_back((int)gen_degs_HA[i + 1].size());
				gen_degs_HA[i + 1].push_back(*p_degs);
				gen_degs_t_HA[i + 1].push_back(gen_degs_HA[i + 1].back().t);
				gen_reprs_HA[i + 1].push_back(gen_reprs_HA[i][p_degs - gen_degs_HA[i].begin()]);
			}
			AddRelsFromGb(gb_HA[i], gen_degs_t_HA[i], gb_HA[i + 1], buffer_HA[i + 1], gen_degs_t_HA[i + 1], map_gen_id[i], t, t_max); /* Add relations of HA[i] to HA[i+1] */
			if (t < t_x) {
				AddRelsFromGb(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_c[i], buffer_HA_ann_c[i], t, t_max);
				AddRelsFromGb(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_y[i], buffer_HA_ann_y[i], t - t_x, t_max);
				leadings_HA[i + 1].clear();
				leadings_HA[i + 1].resize(gen_degs_HA[i + 1].size());
				for (const alg::Poly& g : gb_HA[i + 1])
					leadings_HA[i + 1][g[0][0].gen].push_back(g[0]);
				basis_t = ExtendBasis(gen_degs_HA[i + 1], leadings_HA[i + 1], basis_HA[i + 1], t);
			}
			else {
				y_new = ExtendAnn(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_c[i], buffer_HA_ann_c[i], { c[i] }, { t_x }, t, t_max);
				Indecomposables(gb_HA[i], gen_degs_t_HA[i], gb_HA_ind_y[i], buffer_HA_ind_y[i], y_new, { t_x }, t, t_max);
				for (auto& yi : y_new) {
					y[i].push_back(yi[0]);
					t_y[i].push_back(t - t_x);
				}
				save_y(db, table_prefix + "_y", y_new, t - t_x);

				alg::Poly2d a = ExtendAnn(gb_HA[i], gen_degs_t_HA[i], gb_HA_ann_y[i], buffer_HA_ann_y[i], y[i], t_y[i], t - t_x, t_max - t_x);
				Indecomposables(gb_HA[i], gen_degs_t_HA[i], gb_HA_ind_a[i], buffer_HA_ind_a[i], a, t_y[i], t - t_x, t_max - t_x);

				/* Add new generators to HA[i+1] */
				alg::Poly1d cy_inv;
				for (const alg::Poly1d& yk : y_new) {
					alg::Poly cyk = c[i] * yk[0];
					alg::Poly cyk_repr = get_image(cyk, gen_reprs_HA[i], gb_A0);
					cy_inv.push_back(d_inv(cyk_repr, gen_degs_B, gen_diffs_B, gb_A0, basis_A0, basis_X[i]));
				}
				if (t == t_x * 2) { /* Add [x^2] */
					gen_degs_HA[i + 1].push_back(deg_x * 2);
					gen_degs_t_HA[i + 1].push_back(gen_degs_HA[i + 1].back().t);
					gen_reprs_HA[i + 1].push_back({ {{int(index_x + i + 1), 2}} });
				}
				for (size_t j = 0; j < y_new.size(); ++j) { /* Add [xy+d^{-1}(cy)] */
					gen_degs_HA[i + 1].push_back(deg_x + get_deg(y_new[j][0], gen_degs_HA[i]));
					gen_degs_t_HA[i + 1].push_back(gen_degs_HA[i + 1].back().t);
					gen_reprs_HA[i + 1].push_back(alg::Mon{ {int(index_x + i + 1), 1} } * get_image(y_new[j][0], gen_reprs_HA[i], gb_A0) + cy_inv[j]);
				}

				/* Add new relations to HA[i + 1]. Compute degrees of relations first. */
				std::map<int, alg::array> rel_degs; /* Degrees of relations */
				alg::array range_gen_degs_HA = grbn::range((int)gen_degs_HA[i + 1].size());
				alg::array range_g = lina::AddVectors(range_gen_degs_HA, map_gen_id[i]);
				for (size_t j = 0; j < range_g.size(); ++j) {
					alg::Poly check = gen_reprs_HA[i + 1][range_g[j]] + alg::Poly{ {{int(index_x + i + 1), 2}} };
					if (check.empty()) {
						range_g.erase(range_g.begin() + j); /* Want indices of gj. Remove the index of [x^2]. */
						break;
					}
				}
				for (size_t j = 0; j < range_g.size(); ++j)
					for (size_t k = j; k < range_g.size(); ++k) {
						alg::Deg deg_gjgk = gen_degs_HA[i + 1][range_g[j]] + gen_degs_HA[i + 1][range_g[k]]; /* alg::Deg of [xy_j][xy_k] */
						if (deg_gjgk.t == t)
							rel_degs[deg_gjgk.v + 2 * deg_gjgk.s].push_back(deg_gjgk.s);
					}
				for (const alg::Poly1d& ak : a)
					for (size_t j = 0; j < range_g.size(); ++j)
						if (!ak[j].empty()) {
							alg::Deg deg_akgj = get_deg(ak[j], gen_degs_HA[i]) + gen_degs_HA[i + 1][range_g[j]]; /* alg::Deg of a_{kj}[xy_j] */
							rel_degs[deg_akgj.v + 2 * deg_akgj.s].push_back(deg_akgj.s);
							break;
						}
				for (auto& [v1, list_s] : rel_degs) {
					std::sort(list_s.begin(), list_s.end());
					list_s.erase(std::unique(list_s.begin(), list_s.end()), list_s.end());
				}

				/* Compute relations */
				leadings_HA[i + 1].clear();
				leadings_HA[i + 1].resize(gen_degs_HA[i + 1].size());
				for (const alg::Poly& g : gb_HA[i + 1])
					leadings_HA[i + 1][g[0][0].gen].push_back(g[0]);
				basis_t = ExtendBasis(gen_degs_HA[i + 1], leadings_HA[i + 1], basis_HA[i + 1], t);
#ifndef MULTITHREAD
				for (auto& [v2s, v_s] : rel_degs) {
					alg::Poly1d rels = FindRels(basis_t, basis_A0, basis_X[i + 1], gb_A0, gen_reprs_HA[i + 1], t, v2s, v_s);
					for (alg::Poly& rel : rels) // Change return type
						buffer_HA[i + 1][t].push_back(std::make_unique<grbn::PolyBufferEle>(std::move(rel)));
				}
#else
				std::vector<std::future<alg::Poly1d>> futures;
				for (auto& [v2s, v_s] : rel_degs) {
					futures.push_back(std::async(std::launch::async, FindRels, std::ref(basis_t), std::ref(basis_A0), std::ref(basis_X[i + 1]), std::ref(gb_A0), 
						std::ref(gen_reprs_HA[i + 1]), t, v2s, std::ref(v_s)));
				}
				for (size_t j = 0; j < futures.size(); ++j) {
					futures[j].wait();
					std::cout << "t=" << t << ", i=" << i << ", completed thread=" << j + 1 << '/' << futures.size() << "          \r";
					for (auto& rel : futures[j].get())
						buffer_HA[i + 1][t].push_back(std::make_unique<grbn::PolyBufferEle>(std::move(rel)));
				}
#endif
				grbn::AddRelsB(gb_HA[i + 1], buffer_HA[i + 1], gen_degs_t_HA[i + 1], t, t_max);
			}

			/* Save data */
			db.save_generators(table_H_prefix + "_generators", gen_degs_HA[i + 1], gen_reprs_HA[i + 1], gen_degs_HA_t_start);
			db.save_gb(table_H_prefix + "_relations", gb_HA[i + 1].gb, gen_degs_HA[i + 1], gb_HA_t_start);
			db.save_basis(table_H_prefix + "_basis", basis_t); basis_HA[i + 1].merge(basis_t);
			save_map_gen_id(db, table_prefix + "_map_gen_id", map_gen_id[i], map_gen_id_t_start);
			SaveGb(db, table_prefix + "_ann_c_gb", gb_HA_ann_c[i].gb, gen_degs_t_HA[i], { t_x }, t);
			SaveGb(db, table_prefix + "_ind_y_gb", gb_HA_ind_y[i].gb, gen_degs_t_HA[i], { t_x }, t);
			SaveGb(db, table_prefix + "_ann_y_gb", gb_HA_ann_y[i].gb, gen_degs_t_HA[i], t_y[i], t - t_x);
			SaveGb(db, table_prefix + "_ind_a_gb", gb_HA_ind_a[i].gb, gen_degs_t_HA[i], t_y[i], t - t_x);
		}
		db.end_transaction();
	}
}

int main_test_51e8aa0a(int argc, char** argv)
{
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t_benchmark.db)");
	 
	return 0;
}

int main_benchmark_E4t100(int argc, char** argv)
{
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t_benchmark.db)");
	Timer timer;
	int t_max = 100;
	generate_HB(db, t_max, -1, true);
	return 0;
}

int main(int argc, char** argv)
{
#ifdef BENCHMARK
	return main_benchmark_E4t100(argc, argv);
#else
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	Timer timer;
	int t_max = 189;
	generate_HB(db, t_max, -1, false);
	return 0;
#endif
}