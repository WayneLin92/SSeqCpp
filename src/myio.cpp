#include "myio.h"

/*********** FUNCTIONS **********/

void dump_MonPow(std::ostream& sout, const GenPow& p)
{
	sout << "x_";
	if (0 <= p.gen && p.gen < 10)
		sout << p.gen;
	else
		sout << '{' << p.gen << '}';
	if (p.exp > 1) {
		sout << '^';
		if (0 <= p.exp && p.exp < 10)
			sout << p.exp;
		else
			sout << '{' << p.exp << '}';
	}
};