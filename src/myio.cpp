#include "myio.h"

/*********** FUNCTIONS **********/

void dump_MonPow(std::ostream& sout, const alg::GenPow& p)
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
}

void dump_MonPowV2(std::ostream& sout, const alg::GenPow& p, const std::vector<std::string>& gen_names)
{
    sout << gen_names[p.gen];
    if (p.exp > 1) {
        sout << '^';
        if (0 <= p.exp && p.exp < 10)
            sout << p.exp;
        else
            sout << '{' << p.exp << '}';
    }
}

void dump_MonV2(std::ostream& sout, const alg::Mon& mon, const std::vector<std::string>& gen_names)
{
    for (auto i = mon.begin(); i != mon.end(); ++i) {
        dump_MonPowV2(sout, *i, gen_names);
    }
}

void dump_PolyV2(std::ostream& sout, const alg::Poly& poly, const std::vector<std::string>& gen_names)
{
    for (auto i = poly.begin(); i != poly.end(); ++i) {
        if (i != poly.begin())
            sout << '+';
        dump_MonV2(sout, *i, gen_names);
    }
}