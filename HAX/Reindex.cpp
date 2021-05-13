#include "main_E4t.h"
#include "benchmark.h"

Poly reindex_v2(const Poly& poly, array map_gen_id)
{
	Poly result;
	for (const Mon& m : poly) {
		Mon m1;
		for (GenPow ge : m)
			m1.push_back({ map_gen_id[ge.gen], ge.exp });
		std::sort(m1.begin(), m1.end(), [](const GenPow& lhs, const GenPow& rhs) {return lhs.gen < rhs.gen; });
		result.push_back(m1);
	}
	std::sort(result.begin(), result.end());
	return result;
}

/* This function reorders the generator by t and reproduce the Groebner basis */
std::pair<array, Poly1d> ReorderGens(const std::vector<Deg>& gen_degs, const Poly1d& gb, int t_max)
{
	array map_gen_id_inv = grbn::range((int)gen_degs.size()); /* the i`th new generator is the old map_gen_id_inv[i]`th generator */
	std::sort(map_gen_id_inv.begin(), map_gen_id_inv.end(), [&gen_degs](int i, int j) {return gen_degs[i].t < gen_degs[j].t; });
	array map_gen_id; map_gen_id.resize(gen_degs.size()); /* the i`th old generator becomes the map_gen_id[i]`th generator */
	array gen_degs_new;
	for (int i = 0; i < (int)gen_degs.size(); ++i) {
		map_gen_id[map_gen_id_inv[i]] = i;
		gen_degs_new.push_back(gen_degs[map_gen_id_inv[i]].t);
	}

	grbn::GbBuffer buffer;
	for (const Poly& g : gb)
		buffer[get_deg_t(g, gen_degs)].push_back(reindex_v2(g, map_gen_id));
	Poly1d gb_new;
	grbn::AddRelsB(gb_new, buffer, gen_degs_new, -1, t_max);
	return std::make_pair(std::move(map_gen_id_inv), std::move(gb_new));
}

void ReorderHA()
{
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	int t_max = 200;

	std::vector<Deg> gen_degs_HA = db.load_gen_degs("HA_generators");
	Poly1d gen_reprs_HA = db.load_gen_reprs("HA_generators");
	Poly1d gb_HA = db.load_gb("HA_relations", t_max);
	auto [map_gen_id_inv, gb_new] = ReorderGens(gen_degs_HA, gb_HA, t_max);

	std::vector<Deg> gen_degs_HA_new;
	Poly1d gen_reprs_HA_new;
	for (int i = 0; i < (int)gen_degs_HA.size(); ++i) {
		gen_degs_HA_new.push_back(gen_degs_HA[map_gen_id_inv[i]]);
		gen_reprs_HA_new.push_back(gen_reprs_HA[map_gen_id_inv[i]]);
	}

	db.execute_cmd("DELETE FROM HA_generators_ordered;");
	db.execute_cmd("DELETE FROM HA_relations_ordered;");

	db.begin_transaction();
	db.save_generators("HA_generators_ordered", gen_degs_HA_new, gen_reprs_HA_new);
	db.save_gb("HA_relations_ordered", gb_new, gen_degs_HA_new);
	db.end_transaction();
}

void DeleteDecomposableGenerators()
{
	Timer timer;
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	std::string table_name = "HA3";
	int t_max = 189;

	std::vector<Deg> gen_degs = db.load_gen_degs(table_name + "_generators");
	Poly1d gen_reprs = db.load_gen_reprs(table_name + "_generators");
	grbn::GbWithCache gb;
	gb = db.load_gb(table_name + "_relations", 189);

	std::vector<Deg> gen_degs1;
	array gen_degs1_t;
	Poly1d gen_reprs1;
	Poly1d gb1;

	array indices_decomposables;
	array map_gen_id(gen_degs.size());
	int count = 0;
	for (int i = 0; i < (int)gen_degs.size(); ++i) {
		if (grbn::Reduce({ {{i, 1}} }, gb) != Poly{ {{i, 1}} }) {
			indices_decomposables.push_back(i);
			++count;
			map_gen_id[i] = -1;
		}
		else {
			map_gen_id[i] = i - count;
			gen_degs1.push_back(gen_degs[i]);
			gen_degs1_t.push_back(gen_degs[i].t);
			gen_reprs1.push_back(gen_reprs[i]);
		}
	}
	std::cout << "num of decomposables = " << indices_decomposables.size() << '\n';

	for (const Poly& g : gb)
		if (!std::binary_search(indices_decomposables.begin(), indices_decomposables.end(), g[0][0].gen))
			gb1.push_back(reindex(g, map_gen_id));

	std::string table_name1 = "E4t";

	db.execute_cmd("DELETE FROM " + table_name1 + "_generators;");
	db.execute_cmd("DELETE FROM " + table_name1 + "_relations;");

	db.begin_transaction();
	db.save_generators(table_name1 + "_generators", gen_degs1, gen_reprs1);
	db.save_gb(table_name1 + "_relations", gb1, gen_degs1);
	db.end_transaction();
}

int main_generate_Reindex(int argc, char** argv)
{
	DeleteDecomposableGenerators();
	return 0;
}