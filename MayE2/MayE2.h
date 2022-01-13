#ifndef _MAY_E2_H_
#define _MAY_E2_H_

#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/myexception.h"
#include <iostream>

std::vector<std::pair<alg::array, alg::array>> H(int a, int b);
std::string get_name(const alg::array& S);

void ReorderHX(int n, int s_max);
int bit_length(int t);

void GetD2(int n);

#endif /* _MAY_E2_H_ */