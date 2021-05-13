#ifndef _MAY_E2_H_
#define _MAY_E2_H_

#include "algebras/groebner.h"
#include "algebras/database.h"
#include <iostream>
#include <future>

//#define GENERATE_E4T_1


/********** STRUCTS AND CLASSES **********/

struct DgaBasis1
{
	Mon1d basis;
	Poly1d diffs;
};

/********** FUNCTIONS **********/
template<typename GbType>
Poly get_image(const Poly& poly, const Poly1d& gen_reprs, const GbType& gb)
{
	return grbn::evaluate(poly, [&gen_reprs](int i) {return gen_reprs[i]; }, gb);
}

std::map<Deg, DgaBasis1> load_dga_basis(const Database& db, const std::string& table_name, int r, int t_max);
std::map<Deg, DgaBasis1> load_basis_X(const Database& db, const std::string& table_name, int t_max, int r);
void load_y(const Database& db, const std::string& table_name, Poly1d& y, array& t_y);
void save_y(const Database& db, const std::string& table_name, const Poly2d& y_t, int t);
void save_map_gen_id(const Database& db, const std::string& table_name, const array& map_gen_id, int i_start);
void SaveGb(const Database& db, const std::string& table_name, const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t);
void get_basis_B(const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, Mon1d& basis_B, const Deg& deg);
void get_basis_with_diff_B(const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, Mon1d& basis_B, Poly1d& mon_diffs_B, const Deg& deg);
/* Assume poly is a boundary. Return the chain with it as boundary */
Poly d_inv(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& diffs, const grbn::GbWithCache& gb, const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X);
/* Assume poly is a cycle. Return the homology class */
Poly proj(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, const grbn::GbWithCache& gb, const std::map<Deg, DgaBasis1>& basis_A,
	const std::map<Deg, DgaBasis1>& basis_X, const Poly1d& gen_reprs, std::map<Deg, Mon1d>& basis_H);

/* assume that the sequence map_gen_id is increasing */
Poly reindex(const Poly& poly, const array& map_gen_id);

#endif /* _MAY_E2_H_ */