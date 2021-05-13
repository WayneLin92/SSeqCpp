#include "main_E4t.h"
#include "myexception.h"

/* return a basis for polynomails of xi for t<=t_max */
std::map<Deg, DgaBasis1> get_basis_X(const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, int start_x, int num_x, int t_max)
{
	std::map<Deg, DgaBasis1> result;
	result[Deg{ 0, 0, 0 }].basis.push_back({});
	result[Deg{ 0, 0, 0 }].diffs.push_back({});

	for (int t = 1; t <= t_max; ++t) {
		std::map<Deg, DgaBasis1> basis_new;
		std::cout << "Computing basis_x, t=" << t << "          \r";
		for (int gen_id = start_x + num_x - 1; gen_id >= start_x; --gen_id) { /* 58 is the gen_id of b1 */
			int t1 = t - gen_degs[gen_id].t;
			if (t1 >= 0) {
				auto p1 = result.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = result.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					auto p_m = p->second.basis.begin();
					auto p_diff = p->second.diffs.begin();
					for (; p_m != p->second.basis.end(); ++p_m, ++p_diff) {
						if (p_m->empty() || gen_id <= p_m->front().gen) {
							Mon mon(mul(*p_m, { {gen_id, 1} }));

							int s = p->first.s + gen_degs[gen_id].s;
							int v = p->first.v + gen_degs[gen_id].v;
							basis_new[Deg{ s, t, v }].basis.push_back(std::move(mon));
							Poly diff = add(mul(gen_diffs[gen_id], *p_m), mul({ {gen_id, 1} }, *p_diff));
							basis_new[Deg{ s, t, v }].diffs.push_back(std::move(diff));
						}
					}
				}
			}
		}
		result.merge(basis_new);
	}
	return result;
}

void save_basis_X(Database& db, int start_x, int num_x)
{
	int t_max = 200;
	std::vector<Deg> gen_degs_B = db.load_gen_degs("B_generators");
	Poly1d diffs_B = db.load_gen_diffs("B_generators");
	std::map<Deg, DgaBasis1> basis_X = get_basis_X(gen_degs_B, diffs_B, start_x, num_x, t_max);

	std::string table_name = "X" + std::to_string(num_x) + "_basis";
	try { db.execute_cmd("DROP TABLE " + table_name + ';'); }
	catch (MyException&) {}
	db.execute_cmd("CREATE TABLE " + table_name + " (mon_id INTEGER PRIMARY KEY, mon TEXT NOT NULL UNIQUE, diff TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
	Statement stmt(db, "INSERT INTO " + table_name + " (s, t, v, mon, diff) VALUES (?1, ?2, ?3, ?4, ?5);");

	db.begin_transaction();
	for (auto& [deg, basis_d] : basis_X) {
		for (size_t i = 0; i < basis_d.basis.size(); ++i) {
			stmt.bind_int(1, deg.s);
			stmt.bind_int(2, deg.t);
			stmt.bind_int(3, deg.v);
			stmt.bind_str(4, Mon_to_str(basis_d.basis[i]));
			stmt.bind_str(5, Poly_to_str(basis_X[deg].diffs[i]));
			stmt.step_and_reset();
		}
	}
	db.end_transaction();
	std::cout << table_name << " is created, number of degrees=" << basis_X.size() << '\n';
}


int main_generate_X_basis(int argc, char** argv)
{
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\HB.db)");

	for (int n = 1; n <= 7; n++) /* Create Xi_basis */
		save_basis_X(db, 58, n);
	return 0;
}