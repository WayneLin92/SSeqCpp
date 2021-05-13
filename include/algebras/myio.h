#ifndef MYIO_H
#define MYIO_H

#include "algebras.h"
#include "myparser.h"

/********** STRUCTS AND CLASSES **********/

/********** FUNCTIONS **********/

/* arrays */
inline void load_array(std::istream& sin, array& mon) { load_vector(sin, mon, "(", ",", ")"); }
inline void load_array2d(std::istream& sin, array2d& poly) { load_vector(sin, poly, "[", ",", "]", load_array); }
inline void load_array3d(std::istream& sin, array3d& polys) { load_vector(sin, polys, "{", ",", "}", load_array2d); }
inline void dump_array(std::ostream& sout, const array& mon) { dump_vector(sout, mon, "(", ", ", ")"); }
inline void dump_array2d(std::ostream& sout, const array2d& poly) { dump_vector(sout, poly, "[", ", ", "]", dump_array); }
inline void dump_array3d(std::ostream& sout, const array3d& polys) { dump_vector(sout, polys, "{", ", ", "}", dump_array2d); }

inline std::istream& operator>>(std::istream& sin, array& mon) { load_array(sin, mon); return sin; }
inline std::istream& operator>>(std::istream& sin, array2d& poly) { load_array2d(sin, poly); return sin; }
inline std::istream& operator>>(std::istream& sin, array3d& polys) { load_array3d(sin, polys); return sin; }

inline std::ostream& operator<<(std::ostream& sout, const array& mon) { dump_array(sout, mon); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array2d& poly) { dump_array2d(sout, poly); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array3d& polys) { dump_array3d(sout, polys); return sout; }

/* algebras */
void dump_MonPow(std::ostream& sout, const GenPow& p);
inline void dump_Mon(std::ostream& sout, const Mon& mon) { dump_vector(sout, mon, "", "", "", dump_MonPow); }
inline void dump_Poly(std::ostream& sout, const Poly& poly) { if (poly.empty()) sout << '0'; else dump_vector(sout, poly, "", "+", "", dump_Mon); }
inline void dump_Poly1d(std::ostream& sout, const Poly1d& polys) { dump_vector(sout, polys, "(", ", ", ")", dump_Poly); }
inline void dump_Poly2d(std::ostream& sout, const Poly2d& polys) { dump_vector(sout, polys, "[", ", ", "]", dump_Poly1d); }

inline std::ostream& operator<<(std::ostream& sout, const Deg& d) { sout << '(' << d.s << ',' << d.t << ',' << d.v << ')'; return sout; }
inline std::ostream& operator<<(std::ostream& sout, const GenPow& p) { dump_MonPow(sout, p); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Mon& mon) { dump_Mon(sout, mon); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly& poly) { dump_Poly(sout, poly); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly1d& polys) { dump_Poly1d(sout, polys); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly2d& polys) { dump_Poly2d(sout, polys); return sout; }

#endif /* MYIO_H */