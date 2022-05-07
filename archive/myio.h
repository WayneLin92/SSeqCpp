

/* algebras */
template <typename FnNames>
void DumpGenPow(std::ostream& sout, const alg::GenPow& p, FnNames gen_names)
{
    sout << gen_names(p.gen);
    if (p.exp > 1) {
        sout << '^';
        if (0 <= p.exp && p.exp < 10)
            sout << p.exp;
        else
            sout << '{' << p.exp << '}';
    }
}
template <typename FnNames>
void DumpMon(std::ostream& sout, const alg::Mon& mon, FnNames gen_names)
{
    for (auto i = mon.begin(); i != mon.end(); ++i) {
        DumpGenPow(sout, *i, gen_names);
    }
}
template <typename FnNames>
void DumpPoly(std::ostream& sout, const alg::Mon1d& poly, FnNames gen_names)
{
    if (poly.empty()) {
        sout << '0';
        return;
    }
    for (auto i = poly.begin(); i != poly.end(); ++i) {
        if (i != poly.begin())
            sout << '+';
        DumpMon(sout, *i, gen_names);
    }
}

inline void DumpGenPow(std::ostream& sout, const alg::GenPow& p, const std::vector<std::string>& gen_names)
{
    DumpGenPow(sout, p, [&gen_names](int i) { return gen_names[i]; });
}
inline void DumpMon(std::ostream& sout, const alg::Mon& mon, const std::vector<std::string>& gen_names)
{
    DumpMon(sout, mon, [&gen_names](int i) { return gen_names[i]; });
}
inline void DumpPoly(std::ostream& sout, const alg::Mon1d& poly, const std::vector<std::string>& gen_names)
{
    DumpPoly(sout, poly, [&gen_names](int i) { return gen_names[i]; });
}

void DumpGenPowV2(std::ostream& sout, const alg::GenPow& p);
inline void DumpMonV2(std::ostream& sout, const alg::Mon& mon)
{
    dump_vector(sout, mon, "", "", "", DumpGenPowV2);
}
inline void DumpPolyV2(std::ostream& sout, const alg::Mon1d& poly)
{
    if (poly.empty())
        sout << '0';
    else
        dump_vector(sout, poly, "", "+", "", DumpMonV2);
}

template <typename FnNames>
std::string StrPoly(const alg::Mon1d& poly, FnNames gen_names)
{
    std::ostringstream s;
    DumpPoly(s, poly, gen_names);
    return s.str();
}
inline std::string StrPoly(const alg::Mon1d& poly, const std::vector<std::string>& gen_names)
{
    return StrPoly(poly, [&gen_names](int i) { return gen_names[i]; });
}
inline constexpr int GEN_IDEAL = 0x4000000;
inline constexpr int GEN_SUB = 0x2000000;
inline std::string StrPolyH(const alg::Mon1d& poly, const std::vector<std::string>& gen_names)
{
    return StrPoly(poly, [&gen_names](int i) { return ((i & GEN_IDEAL) ? "b_{" + std::to_string(i - GEN_IDEAL - GEN_SUB) + "}" : ((i & GEN_SUB) ? "x_{" + std::to_string(i - GEN_SUB) + "}" : gen_names[i])); });
}

inline std::ostream& operator<<(std::ostream& sout, const alg::MayDeg& d)
{
    sout << '(' << d.s << ',' << d.t << ',' << d.v << ')';
    return sout;
}
inline std::ostream& operator<<(std::ostream& sout, const alg::GenPow& p)
{
    DumpGenPowV2(sout, p);
    return sout;
}
inline std::ostream& operator<<(std::ostream& sout, const alg::Mon& mon)
{
    DumpMonV2(sout, mon);
    return sout;
}
inline std::ostream& operator<<(std::ostream& sout, const alg::Mon1d& poly)
{
    DumpPolyV2(sout, poly);
    return sout;
}

template <typename FnCmp>
inline alg::Polynomial<FnCmp> GetPolyByName(const std::vector<std::string>& gen_names, const std::string& gn)
{
#ifndef NDEBUG
    if (std::find(gen_names.begin(), gen_names.end(), gn) == gen_names.end())
        throw MyException(0xa98af59aU, "Generator not found");
#endif
    return alg::Polynomial<FnCmp>::Gen((int)(std::find(gen_names.begin(), gen_names.end(), gn) - gen_names.begin()));
}

void DumpGenPowV2(std::ostream& sout, const alg::GenPow& p)
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