#ifndef MYIO_H
#define MYIO_H

#include "algebras.h"
#include "myparser.h"

// clang-format off
/* arrays */
inline void load_array(std::istream& sin, alg::array& mon) { load_vector(sin, mon, "(", ",", ")"); }
inline void load_array2d(std::istream& sin, alg::array2d& poly) { load_vector(sin, poly, "[", ",", "]", load_array); }
inline void load_array3d(std::istream& sin, alg::array3d& polys) { load_vector(sin, polys, "{", ",", "}", load_array2d); }
inline void dump_array(std::ostream& sout, const alg::array& mon) { dump_vector(sout, mon, "(", ", ", ")"); }
inline void dump_array2d(std::ostream& sout, const alg::array2d& poly) { dump_vector(sout, poly, "[", ", ", "]", dump_array); }
inline void dump_array3d(std::ostream& sout, const alg::array3d& polys) { dump_vector(sout, polys, "{", ", ", "}", dump_array2d); }

inline std::istream& operator>>(std::istream& sin, alg::array& mon) { load_array(sin, mon); return sin; }
inline std::istream& operator>>(std::istream& sin, alg::array2d& poly) { load_array2d(sin, poly); return sin; }
inline std::istream& operator>>(std::istream& sin, alg::array3d& polys) { load_array3d(sin, polys); return sin; }

inline std::ostream& operator<<(std::ostream& sout, const alg::array& mon) { dump_array(sout, mon); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const alg::array2d& poly) { dump_array2d(sout, poly); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const alg::array3d& polys) { dump_array3d(sout, polys); return sout; }

/* algebras */
void dump_MonPow(std::ostream& sout, const alg::GenPow& p);
inline void dump_Mon(std::ostream& sout, const alg::Mon& mon) { dump_vector(sout, mon, "", "", "", dump_MonPow); }
inline void dump_Poly(std::ostream& sout, const alg::Poly& poly) { if (poly.empty()) sout << '0'; else dump_vector(sout, poly, "", "+", "", dump_Mon); }
inline void dump_Poly1d(std::ostream& sout, const alg::Poly1d& polys) { dump_vector(sout, polys, "(", ", ", ")", dump_Poly); }
inline void dump_Poly2d(std::ostream& sout, const alg::Poly2d& polys) { dump_vector(sout, polys, "[", ", ", "]", dump_Poly1d); }


void dump_MonPowV2(std::ostream& sout, const alg::GenPow& p, const std::vector<std::string>& gen_names);
void dump_MonV2(std::ostream& sout, const alg::Mon& mon, const std::vector<std::string>& gen_names);
void dump_PolyV2(std::ostream& sout, const alg::Poly& poly, const std::vector<std::string>& gen_names);

inline std::ostream& operator<<(std::ostream& sout, const alg::Deg& d) { sout << '(' << d.s << ',' << d.t << ',' << d.v << ')'; return sout; }
inline std::ostream& operator<<(std::ostream& sout, const alg::GenPow& p) { dump_MonPow(sout, p); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const alg::Mon& mon) { dump_Mon(sout, mon); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const alg::Poly& poly) { dump_Poly(sout, poly); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const alg::Poly1d& polys) { dump_Poly1d(sout, polys); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const alg::Poly2d& polys) { dump_Poly2d(sout, polys); return sout; }

inline alg::Poly GetPolyByName(const std::vector<std::string>& gen_names, const std::string& gn)
{
    return alg::Poly{{{(int)(std::find(gen_names.begin(), gen_names.end(), gn) - gen_names.begin()), 1}}};
}

// clang-format on

#endif /* MYIO_H */