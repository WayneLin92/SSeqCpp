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
	alg::Mon1d basis;
	alg::Poly1d diffs;
};

/********** FUNCTIONS **********/
template<typename GbType>
alg::Poly get_image(const alg::Poly& poly, const alg::Poly1d& gen_reprs, const GbType& gb)
{
	return alg::subs(poly, [&gen_reprs](int i) {return gen_reprs[i]; }, gb);
}

std::map<alg::Deg, DgaBasis1> load_dga_basis(const Database& db, const std::string& table_name, int r, int t_max);
std::map<alg::Deg, DgaBasis1> load_basis_X(const Database& db, const std::string& table_name, int t_max, int r);
void load_y(const Database& db, const std::string& table_name, alg::Poly1d& y, alg::array& t_y);
void save_y(const Database& db, const std::string& table_name, const alg::Poly2d& y_t, int t);
void save_map_gen_id(const Database& db, const std::string& table_name, const alg::array& map_gen_id, int i_start);
void SaveGb(const Database& db, const std::string& table_name, const alg::Poly1d& gb, const alg::array& gen_degs, const alg::array& gen_degs1, int t);
void get_basis_B(const std::map<alg::Deg, DgaBasis1>& basis_A, const std::map<alg::Deg, DgaBasis1>& basis_X, alg::Mon1d& basis_B, const alg::Deg& deg);
void get_basis_with_diff_B(const std::map<alg::Deg, DgaBasis1>& basis_A, const std::map<alg::Deg, DgaBasis1>& basis_X, alg::Mon1d& basis_B, alg::Poly1d& mon_diffs_B, const alg::Deg& deg);
/* Assume poly is a boundary. Return the chain with it as boundary */
alg::Poly d_inv(const alg::Poly& poly, const std::vector<alg::Deg>& gen_degs, const alg::Poly1d& diffs, const alg::Groebner& gb, const std::map<alg::Deg, DgaBasis1>& basis_A, const std::map<alg::Deg, DgaBasis1>& basis_X);
/* Assume poly is a cycle. Return the homology class */
alg::Poly proj(const alg::Poly& poly, const std::vector<alg::Deg>& gen_degs, const alg::Poly1d& gen_diffs, const alg::Groebner& gb, const std::map<alg::Deg, DgaBasis1>& basis_A,
	const std::map<alg::Deg, DgaBasis1>& basis_X, const alg::Poly1d& gen_reprs, std::map<alg::Deg, alg::Mon1d>& basis_H);

/* assume that the sequence map_gen_id is increasing */
alg::Poly reindex(const alg::Poly& poly, const alg::array& map_gen_id);

#endif /* _MAY_E2_H_ */