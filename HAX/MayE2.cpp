#include "MayE2.h"
#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/linalg.h"
#include "algebras/myio.h"

/**
 * The algebra X can have infinitely many gradings with values in Z[Z].
 * We call an element of Z[Z] a signature.
 */
struct Signature
{
    alg::array deg;
    Signature() = default;
    Signature(int s, int t)
    {
        for (int i = 0; i < s; ++i)
            deg.push_back(0);
        for (int i = 0; i < t - s; ++i)
            deg.push_back(1);
    }
    auto size() const
    {
        return deg.size();
    }
    auto& operator[](size_t index) const
    {
        return deg[index];
    }
    Signature operator+(const Signature& rhs) const
    {
        Signature result;
        if (this->size() > rhs.size()) {
            result.deg = this->deg;
            for (size_t i = 0; i < rhs.deg.size(); ++i) {
                result.deg[i] += rhs[i];
            }
        }
        else {
            result.deg = rhs.deg;
            for (size_t i = 0; i < this->size(); ++i) {
                result.deg[i] += this->deg[i];
            }
        }
        return result;
    }
    Signature operator-(const Signature& rhs) const
    {
        Signature result;
        if (this->size() >= rhs.size()) {
            result.deg = this->deg;
            for (size_t i = 0; i < rhs.deg.size(); ++i) {
                result.deg[i] -= rhs[i];
            }
        }
        else {
            std::cerr << "Error: "
                      << "a=" << deg << " "
                      << "b=" << rhs.deg << '\n';
            throw MyException(0x15205197U, "a-b=negative signature.");
        }
        while (!result.deg.empty() && !result.deg.back())
            result.deg.pop_back();
        return result;
    }
    Signature operator*(int rhs) const
    {
        Signature result;
        if (rhs < 0)
            throw MyException(0x822f8715U, "Negative signature.");
        else if (rhs == 0)
            return result;
        result.deg = this->deg;
        for (size_t i = 0; i < this->size(); ++i)
            result.deg[i] *= rhs;
        return result;
    }
    int operator/(const Signature& rhs) const
    {
        if (this->size() < rhs.size()) {
            return 0;
        }
        else {
            int min = INT_MAX;
            for (size_t i = 0; i < rhs.size(); ++i) {
                if (rhs[i] && this->deg[i] / rhs[i] < min)
                    min = this->deg[i] / rhs[i];
            }
            return min;
        }
    }
    bool operator==(const Signature& rhs) const
    {
        return this->deg == rhs.deg;
    }
};

std::ostream& operator<<(std::ostream& sout, const Signature& sig)
{
    sout << sig.deg;
    return sout;
}

/* m is a sparse vector while [p1, p2) is a dense vector starting with index `index` */
bool divisible(const alg::Mon& m, alg::array::const_iterator p1, alg::array::const_iterator p2, size_t index)
{
    for (alg::MonInd p = m.begin(); p != m.end(); ++p)
        if (p->exp > *(p1 + (size_t(p->gen) - index)))
            return false;
    return true;
}

/**
 * Return a basis in degree `deg`.
 * @param leadings Leadings[i] is all monomials starting with generator g_i.
 * @param gen_sigs This gives the signatures of generators and also the number of generators considered in this function.
 * @param sig The signature of the basis.
 */
alg::Mon1d get_basis(const alg::Mon2d& leadings, const std::vector<Signature>& gen_sigs, Signature sig)
{
    alg::Mon1d result;
    alg::array mon_dense;
    mon_dense.resize(gen_sigs.size());
    size_t index = mon_dense.size() - 1;

    while (true) {
        mon_dense[index] = INT_MAX;  
        int e_max = INT_MAX;
        for (const alg::Mon& lead : leadings[index])
            if (divisible(lead, mon_dense.begin() + index, mon_dense.end(), index) && lead[0].exp - 1 < e_max)
                e_max = lead[0].exp - 1;
        e_max = std::min({e_max, sig / gen_sigs[index]});
        mon_dense[index] = e_max;
        sig = sig - gen_sigs[index] * e_max;

        bool move_right = false;
        if (!sig.deg.empty()) {
            if (index > 0)
                index--;
            else
                move_right = true;
        }
        else {
            alg::Mon mon;
            for (size_t i = index; i < mon_dense.size(); ++i)
                if (mon_dense[i])
                    mon.emplace_back(int(i), mon_dense[i]);
            result.push_back(std::move(mon));
            if (index > 0 && e_max > 0) {
                sig = sig + gen_sigs[index];  //
                mon_dense[index--]--;
            }
            else
                move_right = true;
        }
        if (move_right) {
            for (sig = sig + gen_sigs[index] * mon_dense[index], ++index; index < mon_dense.size() && mon_dense[index] == 0; ++index)
                ;
            if (index == mon_dense.size()) {
                break;
            }
            else {
                sig = sig + gen_sigs[index];
                mon_dense[index--]--;
            }
        }
    }
    return result;
}

std::vector<std::pair<alg::array, alg::array>> Hp(int a, int b);

/**
 * Iterator of (S, T) such that $h_{ST}\in \sH_{a,b-1}$.
 */
std::vector<std::pair<alg::array, alg::array>> H(int a, int b)
{
    std::vector<std::pair<alg::array, alg::array>> result;
    if (b - a == 0 || (b - a) % 2)
        return result;
    for (auto [S, T] : Hp(a + 1, b - 1)) {
        S.insert(S.begin(), a);
        T.insert(T.end(), b - 1);
        result.push_back(std::make_pair(std::move(S), std::move(T)));
    }
    return result;
}

/**
 * Iterator of (S, T) such that $h_{ST}\in \sH^\prime_{a,b-1}$.
 */
std::vector<std::pair<alg::array, alg::array>> Hp(int a, int b)
{
    std::vector<std::pair<alg::array, alg::array>> result;
    if ((b - a) % 2)
        return result;
    if (b - a == 0) {
        result.push_back(std::make_pair(alg::array{}, alg::array{}));
        return result;
    }
    for (int i = a + 2; i < b + 2; i += 2) {
        for (auto [S1, T1] : H(a, i)) {
            for (auto [S2, T2] : Hp(i, b)) {
                alg::array S, T;
                std::merge(S1.begin(), S1.end(), S2.begin(), S2.end(), std::back_inserter(S));
                std::merge(T1.begin(), T1.end(), T2.begin(), T2.end(), std::back_inserter(T));
                result.push_back(std::make_pair(S, T));
            }
        }
    }
    return result;
}

/**
 * Return the short name for $h_{S, T}$.
 *
 * Here we assume that len(S) >= 1;
 */
std::string get_name(const alg::array& S)
{
    if (S.size() < 1)
        throw MyException(0xc56fc560, "get_name(S): S.size() < 1)");
    std::string result = "h_" + std::to_string(S[0]);
    if (S.size() > 1) {
        result += "(";
        for (size_t i = 1; i < S.size(); ++i) {
            if (i < S.size() - 1)
                result += std::to_string(S[i] - S[0]) + ",";
            else
                result += std::to_string(S[i] - S[0]) + ")";
        }
    }
    return result;
}

/**
 * Return the short name for $h_{S, T - {max(T)} + {n}}$.
 *
 * Here we assume that len(S) >= 2;
 */
std::string get_name(const alg::array& S, int n)
{
    std::string result = "r_{" + std::to_string(S[0]) + std::to_string(n) + "}";
    if (S.size() > 1) {
        result += "(";
        for (size_t i = 1; i < S.size(); ++i) {
            if (i < S.size() - 1)
                result += std::to_string(S[i] - S[0]) + ",";
            else
                result += std::to_string(S[i] - S[0]) + ")";
        }
    }
    return result;
}

Signature get_sig(const alg::array& S, const alg::array& T)
{
    Signature result;
    for (size_t i = 0; i < S.size(); ++i)
        result = result + Signature(S[i], T[i]);
    return result;
}

int get_id(int i, int j)
{
    return (19 - i) * i / 2 + (j - i - 1);
}

alg::PolyLex get_repr(const alg::array& S, alg::array T)
{
    alg::PolyLex result;
    do {
        bool nonzero = true;
        for (size_t i = 0; i < S.size(); ++i)
            if (S[i] > T[i])
                nonzero = false;
        if (nonzero) {
            alg::Mon m;
            for (size_t i = 0; i < S.size(); ++i)
                m.push_back({get_id(S[i], T[i]), 1});
            result += alg::PolyLex{{std::move(m)}};
        }
    } while (std::next_permutation(T.begin(), T.end()));
    return result;
}

alg::PolyLex get_repr(const alg::array& S, const alg::array& T, int n)
{
    alg::PolyLex result;
    int m = T.back();
    for (int k = m; k < n; ++k) {
        alg::array T1(T);
        T1.back() = k;
        do {
            bool nonzero = true;
            for (size_t i = 0; i < S.size(); ++i)
                if (S[i] > T1[i])
                    nonzero = false;
            if (nonzero) {
                alg::Mon m;
                for (size_t i = 0; i < S.size(); ++i)
                    m.push_back({get_id(S[i], T1[i]), 1});
                m.push_back({get_id(k, n), 1});
                result += alg::PolyLex{{m}};
            }
        } while (std::next_permutation(T1.begin(), T1.end()));
    }
    return result;
}

/**
 * Return the bit length of t.
 */
int bit_length(int t)
{
    int bits;
    unsigned var = (t < 0) ? -t : t;
    for (bits = 0; var != 0; ++bits)
        var >>= 1;
    return bits;
}

/**
 * Compute the signature of repr which is a polynomials in X9.
 * @param poly A polynomial in X9.
 * @param gen_degs_X9 Degrees of generators of X9.
 */
Signature get_sig(const alg::PolyLex& poly, const std::vector<alg::MayDeg>& gen_degs_X9)
{
    Signature result;
    for (alg::MonInd p = poly.GetLead().begin(); p != poly.GetLead().end(); ++p) {
        alg::MayDeg deg = gen_degs_X9[p->gen];
        int j = bit_length(deg.t);
        int i = j - (deg.v + 1);
        result = result + Signature(i, j) * p->exp;
    }
    return result;
}

/**
 * Compute the signature of repr which is a polynomials.
 * @param poly A polynomial.
 * @param gen_sigs Signatures of generators in poly.
 */
Signature get_sig(const alg::PolyRevlex& poly, const std::vector<Signature>& gen_sigs)
{
    Signature result;
    for (alg::MonInd p = poly.GetLead().begin(); p != poly.GetLead().end(); ++p)
        result = result + gen_sigs[p->gen] * p->exp;
    return result;
}

/**
 * Solve the extension problem for `rel`.
 *
 * Search for all basis in A in the right degree and determine for which p such that rel+p=0.
 * If such p does not exist, return false.
 * If such p is not unique, return false.
 * Otherwise, modify rel and return true.
 */
bool solve_extension(alg::PolyRevlex& rel, const std::vector<std::string>& gen_names, const alg::PolyLex1d& gen_repr, const std::vector<Signature>& gen_sigs, const alg::Mon2d& leadings,
                     const std::vector<alg::MayDeg>& gen_degs_X9, const alg::GroebnerLex& gb_B9)
{
    std::cout << "rel=" << StrPoly(rel.data, gen_names) << '\n';

    alg::PolyLex repr_rel = gb_B9.Reduce(alg::subs(rel.data, gen_repr));
    if (!repr_rel)
        return true;
    Signature sig_rel = get_sig(repr_rel, gen_degs_X9);

    // std::cout << "leadings=" << leadings << '\n';
    // std::cout << "gen_sigs=";
    // dump_vector(std::cout, gen_sigs, "[", ", ", "]");
    // std::cout << '\n';
    // std::cout << "sig_rel=" << sig_rel << '\n';

    alg::Mon1d basis = get_basis(leadings, gen_sigs, sig_rel);

    /*std::cout << "basis=(";
    for (const auto& b : basis) {
        dump_MonV2(std::cout, b, gen_names);
        std::cout << ", ";
    }
    std::cout << ")\n";*/

    // std::cout << "repr_basis=("; //
    lina::array2d fx;
    for (const auto& m : basis) {
        alg::PolyLex repr_m = gb_B9.Reduce(alg::subs({m}, gen_repr));
        // std::cout << repr_m.data << "="; //
        fx.push_back(alg::hash1d(repr_m.data));
        // std::cout << alg::hash1d(repr_m.data) << ", "; //
    }
    // std::cout << ")\n"; //
    // std::cout << "repr_rel=" << repr_rel.data << '=' << alg::hash1d(repr_rel.data) << '\n'; //
    fx.push_back(alg::hash1d(repr_rel.data));

    lina::array2d image, kernel, g;
    lina::SetLinearMap(fx, image, kernel, g);
    // std::cout << "kernel[0]=" << kernel[0] << ' ' << kernel[0].back() << ' ' << basis.size() << '\n';
    if (kernel.size() != 1 || kernel[0].back() != (int)basis.size()) {
        if (kernel.empty())
            std::cerr << "Bug: no solution to the extension problem.\n";
        else {
            if (kernel[0].back() != (int)basis.size()) {
                alg::PolyRevlex k;
                for (size_t i = 0; i < kernel[0].size(); ++i)
                    k += alg::PolyRevlex{{basis[kernel[0][i]]}};

                std::cerr << "k=" << StrPoly(k.data, gen_names) << '\n';
                std::cerr << "repr_k=" << gb_B9.Reduce(alg::subs(k.data, gen_repr)).data << '\n';
                std::cerr << "deg(k)=" << alg::GetDegT(alg::subs(k.data, gen_repr).data.front(), gen_degs_X9) << '\n';
                std::cerr << "Unexpected solution.\n";
            }
            else {
                std::cerr << "basis.size()=" << basis.size() << '\n';
                std::cerr << "kernel[0]=" << kernel[0] << '\n';
                std::cerr << "kernel[1]=" << kernel[1] << '\n';
                alg::PolyRevlex ext1 = rel;
                alg::PolyRevlex ext2 = rel;
                for (size_t i = 0; i < kernel[0].size() - 1; ++i)
                    ext1 += alg::PolyRevlex{{basis[kernel[0][i]]}};
                for (size_t i = 0; i < kernel[1].size() - 1; ++i)
                    ext2 += alg::PolyRevlex{{basis[kernel[1][i]]}};

                std::cerr << "ext1=" << StrPoly(ext1.data, gen_names) << '\n';

                auto ext1_repr = gb_B9.Reduce(alg::subs(ext1.data, gen_repr));
                std::cerr << "Reduce(repr(ext1))=" << ext1_repr.data << '\n';

                std::cerr << "ext2=" << StrPoly(ext2.data, gen_names) << '\n';

                std::cerr << "Bug: nonunique solutions to the extension problem.\n";
            }
        }
        return false;
    }
    for (size_t i = 0; i < kernel[0].size() - 1; ++i)
        rel += alg::PolyRevlex{{basis[kernel[0][i]]}};
    return true;
}

/**
 * Compute HX_{nm} where X_{nm} is the polynomial DGA generated by R<=R_{n-m, n}.
 */
int generate_Xnm(int n, int m, int t_max)
{
    myio::DbAlg db("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/HX9.db");
    std::string table_prefix;
    if (m - 1 == 0)
        table_prefix = std::string("HX") + std::to_string(n - 1);
    else
        table_prefix = std::string("HX") + std::to_string(n) + std::to_string(m - 1);
    std::string table1_prefix = std::string("HX") + std::to_string(n) + std::to_string(m);  //

    std::cout << "A=" << table_prefix << '\n';
    std::cout << "B=" << table1_prefix << '\n';

    db.execute_cmd("CREATE TABLE IF NOT EXISTS " + table1_prefix + "_generators (gen_id INTEGER PRIMARY KEY, gen_name TEXT UNIQUE, gen_diff TEXT, repr TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
    db.execute_cmd("CREATE TABLE IF NOT EXISTS " + table1_prefix + "_relations (leading_term TEXT, basis TEXT, s SMALLINT, t SMALLINT, v SMALLINT);");
    db.execute_cmd("DELETE FROM " + table1_prefix + "_generators");
    db.execute_cmd("DELETE FROM " + table1_prefix + "_relations");

    alg::MayDeg deg_x = {1, (1 << n) - (1 << (n - m)), m - 1};
    int id_x = get_id(n - m, n);
    Signature sig_x = Signature(n - m, n);

    std::cout << "id_x=" << id_x << '\n';

    std::vector<alg::MayDeg> gen_degs = db.load_gen_maydegs(table_prefix + "_generators");
    alg::array gen_degs_t;
    for (auto p = gen_degs.begin(); p < gen_degs.end(); ++p)
        gen_degs_t.push_back(p->t);
    std::vector<std::string> gen_names = db.load_gen_names(table_prefix + "_generators");
    alg::PolyLex1d gen_reprs = db.load_gen_reprs<alg::CmpLex>(table_prefix + "_generators");
    alg::GroebnerRevlex gb = db.load_gb<alg::CmpRevlex>(table_prefix + "_relations", t_max);

    std::vector<alg::MayDeg> gen_degs_X9 = db.load_gen_maydegs("X9_generators");
    alg::PolyLex1d gen_diffs_E1 = db.load_gen_diffs<alg::CmpLex>("X9_generators");

    alg::Groebner gb1 = gb;
    std::vector<Signature> gen_sigs;
    for (size_t i = 0; i < gen_reprs.size(); ++i)
        gen_sigs.push_back(get_sig(gen_reprs[i], gen_degs_X9));

    myio::DbAlg db_B9("C:/Users/lwnpk/Documents/MyData/Math_AlgTop/databases/B9.db");
    alg::GroebnerLex gb_B9 = db_B9.load_gb<alg::CmpLex>("B" + std::to_string(n) + std::to_string(m) + "_relations", t_max);

    if (m == 1) {
        /* Add h_{n-1} */
        gen_degs.push_back(gen_degs_X9[id_x]);
        gen_names.push_back("h_" + std::to_string(n - 1));
        gen_reprs.push_back(alg::PolyLex::Gen(id_x));
    }
    else {
        /* # Compute new generators */
        alg::PolyRevlex c;
        if (m == 2) {
            std::string gen_name = "h_" + std::to_string(n - 2);
            auto p = std::find(gen_names.begin(), gen_names.end(), gen_name);
            if (p == gen_names.end()) {
                std::cerr << "Bug: do not find the class $h_{n-2}$\n";
                std::cerr << "gen_name=" << gen_name << '\n';
                std::cerr << "gen_names=";
                dump_vector(std::cerr, gen_names, "(", ", ", ")");
                std::cerr << '\n';
                return 1;
            }
            int gen_id = int(p - gen_names.begin());
            c = alg::PolyRevlex::Gen(gen_id) * alg::PolyRevlex::Gen((int)gen_degs.size() - 1); /* dx = c = h_{n-2}h_{n-1} */
        }
        else {
            std::string gen_name = "r_{" + std::to_string(n - m) + std::to_string(n) + "}";
            auto p = std::find(gen_names.begin(), gen_names.end(), gen_name);
            if (p == gen_names.end()) {
                std::cerr << "Bug: do not find the boundary class $r_{n-m, n}$\n";
                std::cerr << "gen_name=" << gen_name << '\n';
                std::cerr << "gen_names=";
                dump_vector(std::cerr, gen_names, "(", ", ", ")");
                std::cerr << '\n';
                return 1;
            }
            int gen_id = int(p - gen_names.begin());
            c = alg::PolyRevlex::Gen(gen_id); /* dx = c */
        }

        std::vector<alg::PolyRevlex1d> y2d = alg::AnnSeq(gb, {c}, gen_degs_t, t_max); /* y=ann(c) */
        alg::PolyRevlex1d y;
        for (alg::PolyRevlex1d& v : y2d)
            y.push_back(std::move(v[0]));
        // std::cout << "y=" << y << '\n';
        std::cout << "y.size()=" << y.size() << '\n';

        /* Create the list of expected new indecomposables */
        std::vector<std::string> new_gen_names;
        std::vector<Signature> new_gen_sigs;
        std::vector<alg::PolyLex> new_gen_reprs;

        for (auto [S, T] : H(n - m, n + 1)) {
            new_gen_names.push_back(get_name(S));
            new_gen_sigs.push_back(get_sig(S, T));
            new_gen_reprs.push_back(get_repr(S, T));
        }
        for (int i = n - m - 1; i >= 0; i -= 2) {
            for (auto [S, T] : H(i, n - m + 1)) {
                new_gen_names.push_back(get_name(S, n));
                new_gen_sigs.push_back(get_sig(S, T) + sig_x);
                new_gen_reprs.push_back(get_repr(S, T, n));
            }
        }

        std::cout << "Expected new generators like h_i(S): ";
        dump_vector(std::cout, new_gen_names, "(", ", ", ")");
        std::cout << '\n';

        if (y.size() != new_gen_names.size()) {
            std::cerr << "Bug: unexpected Indecomposables among y.\n";
            for (size_t i = 0; i < y.size(); ++i) {
                std::cerr << "y[" << i << "]=";
                dump_PolyV2(std::cerr, y[i].data, gen_names);
                std::cerr << '\n';
            }
            return 2;
        }

        /* Add [x^2] */
        gen_degs.push_back(gen_degs_X9[id_x] * 2);
        gen_degs_t.push_back(gen_degs.back().t);
        gen_names.push_back("b_{" + std::to_string(n - m) + std::to_string(n) + "}");
        gen_reprs.push_back(alg::PolyLex::GenExp(id_x, 2));

        std::cout << "generator " << gen_names.back() << " added\n";

        /* Add [xy+d^{-1}(cy)] */
        alg::array indices_y; /* Indices in new_gen_names corresponding to y */
        for (auto& yi : y) {
            Signature sig_gi = get_sig(yi, gen_sigs) + sig_x;
            bool found = false;
            for (size_t i = 0; i < new_gen_sigs.size(); ++i) {
                if (sig_gi == new_gen_sigs[i]) {
                    indices_y.push_back((int)i);
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Surprise: unexpected indecomposables:" << yi.data << '\n';
                std::cerr << "sig_gi=" << sig_gi << '\n';
                std::cerr << "new_gen_sigs=";
                dump_vector(std::cerr, new_gen_sigs, "(", ", ", ")");
                std::cerr << '\n';
                std::cerr << "new_gen_names=";
                dump_vector(std::cerr, new_gen_names, "(", ", ", ")");
                std::cerr << '\n';
                return 3;
            }
        }

        for (size_t i = 0; i < indices_y.size(); ++i) {
            gen_degs.push_back(y[i].GetMayDeg(gen_degs) + deg_x);
            gen_degs_t.push_back(gen_degs.back().t);
            gen_names.push_back(new_gen_names[indices_y[i]]);
            gen_reprs.push_back(new_gen_reprs[indices_y[i]]);
        }

        std::cout << indices_y.size() << " generators added\n";

        /* # Compute the relations */

        /* Add relation c=0 */
        std::cout << "Adding c=" << StrPoly(c.data, gen_names) << "=0\n";
        alg::AddRels(gb, {c}, gen_degs_t, t_max);

        alg::Mon2d leadings(gen_sigs.size());  // TODO: group by the last generator
        for (size_t i = 0; i < gb.data.size(); ++i)
            leadings[gb.data[i].GetLead()[0].gen].push_back(gb.data[i].GetLead());

        /* Add relations gi * gj = b * yi * yj + [...] */
        std::cout << "Adding relations gi * gj = b * yi * yj + [...]\n";
        alg::PolyRevlex1d new_relations;
        int id_b = (int)gen_degs.size() - (int)y.size() - 1;
        alg::PolyRevlex b = alg::PolyRevlex::Gen(id_b);
        for (size_t i = 0; i < y.size(); ++i) {
            int id_gi = (int)gen_degs.size() - (int)y.size() + (int)i;
            alg::PolyRevlex gi = alg::PolyRevlex::Gen(id_gi);
            for (size_t j = i; j < y.size(); ++j) {
                int id_gj = (int)gen_degs.size() - (int)y.size() + (int)j;
                if (gen_degs_t[id_gi] + gen_degs_t[id_gj] <= t_max) {
                    alg::PolyRevlex gj = alg::PolyRevlex::Gen(id_gj);
                    alg::PolyRevlex rel = gi * gj + b * y[i] * y[j];
                    if (!solve_extension(rel, gen_names, gen_reprs, gen_sigs, leadings, gen_degs_X9, gb_B9)) {
                        return 4;
                    }
                    else {
                        std::cout << "Extended to rel=" << StrPoly(rel.data, gen_names) << "\n\n";
                        new_relations.push_back(std::move(rel));
                    }
                }
            }
        }

        /* Add relations a1 g1 + ... + ak gk = [...] */
        std::cout << "Adding relations a1 g1 + ... + ak gk = [...]\n";
        std::vector<alg::PolyRevlex1d> a = alg::AnnSeq(gb, y, gen_degs_t, t_max - deg_x.t);
        std::cout << "a.size()=" << a.size() << '\n';

        for (const auto& ai : a) {
            alg::PolyRevlex rel;
            for (size_t j = 0; j < ai.size(); ++j) {
                int id_gj = (int)gen_degs.size() - (int)y.size() + (int)j;
                alg::PolyRevlex gj = alg::PolyRevlex::Gen(id_gj);
                rel += ai[j] * gj;
            }
            if (!solve_extension(rel, gen_names, gen_reprs, gen_sigs, leadings, gen_degs_X9, gb_B9)) {
                return 5;
            }
            else {
                std::cout << "Extended to rel=";
                dump_PolyV2(std::cout, rel.data, gen_names);
                std::cout << "\n\n";
                new_relations.push_back(std::move(rel));
            }
        }

        std::sort(new_relations.begin(), new_relations.end(), [&gen_degs](const alg::PolyRevlex& rel1, const alg::PolyRevlex& rel2) { return rel1.GetMayDeg(gen_degs) < rel2.GetMayDeg(gen_degs); });
        alg::AddRels(gb, new_relations, gen_degs_t, t_max);
    }

    db.begin_transaction();

    db.save_gen_maydegs(table1_prefix + "_generators", gen_degs);
    db.save_gen_names(table1_prefix + "_generators", gen_names);
    db.save_gen_reprs(table1_prefix + "_generators", gen_reprs, 0);
    db.save_gb(table1_prefix + "_relations", gb, gen_degs);

    db.end_transaction();

    return 0;
}

int main(int argc, char** argv)
{
    Timer timer;

    for (int n = 1; n <= 9; ++n) {
        for (int m = 1; m <= n; ++m) {
            if (generate_Xnm(n, m, 1024))
                return 1;
            std::cout << "\n-----------\n";
        }
        ReorderHX(n, 1024);
        std::cout << "\n-----------\n";
        std::cout << "-----------\n";
    }

    //generate_Xnm(9, 3, 1024);
    // ReorderHX(2, 1024);

    return 0;
}
