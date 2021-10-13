/** \file groebner.h
 * A Component for Groebner bases.
 */

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
//#include <iostream> ////
#include <map>
#include <memory>

/**
 * The macro enables the multithreading for `AddRels` and functions it depends on.
 * It can be toggled off.
 */
#define GROEBNER_MULTITHREAD
#ifdef GROEBNER_MULTITHREAD
#include <future>
/**
 * The number of threads used in multithreading.
 */
constexpr int kNumThreads = 4;
#endif

// /**
//  * The macro enables the testing of some templates at compile time.
//  * It can be toggled off.
//  */
// #define GROEBNER_TEMPLATE_INSTANTIATIONS

/* extension of namespace alg in algebras.h */
namespace alg {

/**
 * The auxiliary class `GbBuffer` is a set of
 * polynomials to be added to a groebner basis.
 *
 * The polynomials are grouped by degrees.
 */
using GbBuffer = std::map<int, Poly1d>;

/**
 * A virtual type of elements used in Groebner buffer V2,
 * which is designed to save the memory.
 *
 * The subclass `GcdBufferEle` stores the indices of
 * two elements of `gb` and cache the GCD of them and
 * it represents the new polynomials produced by those
 * two elements.
 *
 * The subclass `PolyBufferEle` directly stores the polynomial
 * to be added to `gb`.
 *
 * @see GbBuffer::BufferV2
 */
template <typename GbType> struct BaseBufferEle
{
    virtual ~BaseBufferEle(){};
    virtual Poly GetPoly(const GbType& gb) = 0;
};
/**
 * @see BaseBufferEle
 */
template <typename GbType> struct GcdBufferEle : BaseBufferEle<GbType>
{
    Mon gcd_;
    int i1_, i2_;
    GcdBufferEle(Mon gcd, int i1, int i2) : gcd_(std::move(gcd)), i1_(i1), i2_(i2){};
    Poly GetPoly(const GbType& gb) override
    {
        Poly result = gb[i1_] * div(GbType::GetLead(gb[i2_]), gcd_) + gb[i2_] * div(GbType::GetLead(gb[i1_]), gcd_);
        Mon{}.swap(gcd_); /* Deallocate */
        return result;
    }
};
/**
 * @see BaseBufferEle
 */
template <typename GbType> struct PolyBufferEle : BaseBufferEle<GbType>
{
    Poly p_;
    PolyBufferEle(Poly p) : p_(std::move(p)){};
    Poly GetPoly(const GbType& gb) override
    {
        return std::move(p_);
    }
};
/**
 * The class `Groebner` consists of a Groebner basis and its cache.
 *
 * The monomial ordering is the reversed lexicographical ordering.
 *
 * The member `lc` is used to speed up the calculation of `Reduce(p, gb)`.
 * `lc[i]` is the array of indices of elements of `gb`
 * whose leading term's last generator id is `i`.
 * This is the most effective if the generators are ordered by degree.
 *
 * TODO: refactor gb with gb_ and make it private
 */
class Groebner
{
public:
    Groebner() = default;
    Groebner(Poly1d gb1) : gb(std::move(gb1))
    {
        for (int i = 0; i < (int)gb.size(); ++i)
            lc[GetLead(gb[i]).back().gen].push_back(i);
    }
    static const Mon& GetLead(const Poly& poly)
    {
        return poly.front();
    }
    using BufferV2 = std::map<int, std::vector<std::unique_ptr<BaseBufferEle<Groebner>>>>;

public: /* Convenient interfaces for the member `gb` */
    auto begin() const
    {
        return gb.begin();
    }
    auto end() const
    {
        return gb.end();
    }
    auto rbegin() const
    {
        return gb.rbegin();
    }
    auto rend() const
    {
        return gb.rend();
    }
    auto size() const
    {
        return gb.size();
    }
    auto& operator[](size_t index) const
    {
        return gb[index];
    }
    void push_back(Poly g)
    {
        lc[GetLead(g).back().gen].push_back((int)gb.size());
        gb.push_back(std::move(g));
    }

public:
    Poly1d gb;
    std::map<int, array> lc;
};

/**
 * The class `GroebnerLex` consists of a Groebner basis and its cache.
 *
 * The monomial ordering is the lexicographical ordering.
 */
class GroebnerLex : public Groebner
{
public:
    GroebnerLex() = default;
    GroebnerLex(Poly1d gb1)
    {
        gb = std::move(gb1);
        for (int i = 0; i < (int)gb.size(); ++i)
            lc[GetLead(gb[i]).back().gen].push_back(i);
    }
    static const Mon& GetLead(const Poly& poly)
    {
        return poly.back();
    };
    using BufferV2 = std::map<int, std::vector<std::unique_ptr<BaseBufferEle<GroebnerLex>>>>;

public:
    void push_back(Poly g)
    {
        lc[GetLead(g).back().gen].push_back((int)gb.size());
        gb.push_back(std::move(g));
    }
};

/**
 * Reduce a polynomial by a Groebner basis
 */
Poly Reduce(Poly poly, const Poly1d& gb);
/**
 * Reduce a polynomial by a Groebner basis (with cache)
 */
Poly Reduce(Poly poly, const Groebner& gb);
/**
 * Reduce a polynomial by a Groebner basis with lexicographical ordering (with cache)
 */
Poly Reduce(Poly poly, const GroebnerLex& gb);

/**
 * Remove elements of `cont` which are empty containers
 */
template <typename Container1d> inline void RemoveEmptyElements(Container1d& cont)
{
    cont.erase(std::remove_if(cont.begin(), cont.end(), [](const typename Container1d::value_type& g) { return g.empty(); }), cont.end());
}

/**
 * Create the int array 1, 2, ..., n
 */
inline array range(int n)
{
    array result;
    for (int i = 0; i < n; ++i)
        result.push_back(i);
    return result;
};

/**
 * A fast algorithm that computes
 * `poly ** n` modulo `gb`
 */
template <typename GbType> Poly pow(const Poly& poly, int n, const GbType& gb)
{
    Poly result = { {} };
    if (n == 0)
        return result;
    Poly power = poly;
    while (n) {
        if (n & 1)
            result = Reduce(mul(result, power), gb);
        n >>= 1;
        if (n)
            power = Reduce(mul(power, power), gb);
    }
    return result;
}

/**
 * Replace the generators in `poly` with elements given in `map`
 * and evaluate modulo `gb`.
 * @param poly The polynomial to be substituted.
 * @param map `map[i]` is the polynomial that substitutes the generator of id `i`.
 * @param gb A Groebner basis.
 */
template <typename FnType, typename GbType> Poly subs(const Poly& poly, FnType map, const GbType& gb)
{
    Poly result;
    for (const Mon& m : poly) {
        Poly fm = { {} };
        for (MonInd p = m.begin(); p != m.end(); ++p)
            fm = Reduce(mul(fm, pow(map(p->gen), p->exp, gb)), gb);
        result = add(result, fm);
    }
    return result;
}

/**
 * Specialization.
 * @see subs(const Poly&, FnType, const GbType&)
 */
template <typename GbType> Poly subs(const Poly& poly, const Poly1d& map, const GbType& gb)
{
    Poly result;
    for (const Mon& m : poly) {
        Poly fm = { {} };
        for (MonInd p = m.begin(); p != m.end(); ++p)
            fm = Reduce(mul(fm, pow(map[p->gen], p->exp, gb)), gb);
        result = add(result, fm);
    }
    return result;
}

/**
 * Return if `mon1` and `mon2` have a nontrivial common factor.
 */
bool gcd_nonzero(const Mon& mon1, const Mon& mon2);

/**
 * Generate buffer in degree `t_min <= t <= t_max`
 */
GbBuffer GenerateBuffer(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t_min, int t_max);
/**
 * Generate buffer in degree `t_min <= t <= t_max`
 */
Groebner::BufferV2 GenerateBufferV2(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t_min, int t_max);

#ifdef GROEBNER_MULTITHREAD
/*
 * Auxiliary function used in multi-threading computing.
 */
template <typename GbType, typename InType> void _BatchReduce(const GbType& gb, const InType& rels_in, Poly1d& rels_out, int i_start)
{
    for (size_t i = i_start; i < rels_in.size(); i += kNumThreads) {
        if constexpr (std::is_same<InType, Poly1d>::value)
            rels_out[i] = Reduce(rels_in[i], gb);
        else /* InType == std::vector<std::unique_ptr<BaseBufferEle>> */
            rels_out[i] = Reduce(rels_in[i]->GetPoly(gb), gb);
    }
}

/*
 * Auxiliary function used in multi-threading computing.
 */
template <typename GbType, typename FnType> std::vector<std::pair<int, Poly>> _BatchNewBufferElements(const GbType& gb, const Poly& rel, FnType _get_deg, int deg_max, int i_start)
{
    std::vector<std::pair<int, Poly>> result;
    for (size_t i = i_start; i < gb.size(); i += kNumThreads) {
        if (gcd_nonzero(GbType::GetLead(rel), GbType::GetLead(gb[i]))
            && !(GbType::GetLead(rel)[0].gen < 0 && GbType::GetLead(gb[i])[0].gen < 0 && GbType::GetLead(rel)[0].gen != GbType::GetLead(gb[i])[0].gen)) {
            Mon lcm = LCM(GbType::GetLead(rel), GbType::GetLead(gb[i]));
            int deg_new_rel = _get_deg(lcm);
            if (deg_max == -1 || deg_new_rel <= deg_max) {
                Poly new_rel = rel * div(lcm, GbType::GetLead(rel)) + (gb[i]) * div(lcm, GbType::GetLead(gb[i]));
                if (!new_rel.empty())
                    result.push_back(std::make_pair(deg_new_rel, std::move(new_rel)));
            }
        }
    }
    return result;
}
/*
 * Auxiliary function used in multi-threading computing.
 * GbType must be ordered with revlex.
 */
template <typename GbType, typename FnType>
std::vector<std::pair<int, GcdBufferEle<GbType>>> _BatchNewBufferElementsV2(const GbType& gb, const Poly& rel, FnType _get_deg, int deg, int deg_max, int i_start)
{
    std::vector<std::pair<int, GcdBufferEle<GbType>>> result;
    for (size_t i = i_start; i < gb.size(); i += kNumThreads) {
        if (gcd_nonzero(GbType::GetLead(rel), GbType::GetLead(gb[i]))
            && !(GbType::GetLead(rel)[0].gen < 0 && GbType::GetLead(gb[i])[0].gen < 0 && GbType::GetLead(rel)[0].gen != GbType::GetLead(gb[i])[0].gen)) {
            Mon gcd = GCD(GbType::GetLead(rel), GbType::GetLead(gb[i]));
            int deg_new_rel = deg + _get_deg(GbType::GetLead(gb[i])) - _get_deg(gcd);
            if (deg_max == -1 || deg_new_rel <= deg_max)
                result.push_back(std::make_pair(deg_new_rel, GcdBufferEle<GbType>(std::move(gcd), (int)i, (int)gb.size())));
        }
    }
    return result;
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _TestBatchReduce(Poly1d& gb, Poly1d& rels_in, Poly1d& rels_out, int i_start)
{
    _BatchReduce(gb, rels_in, rels_out, i_start);
}
inline void _TestBatchReduce(Groebner& gb, std::vector<std::unique_ptr<BaseBufferEle>>& rels_in, Poly1d& rels_out, int i_start)
{
    _BatchReduce(gb, rels_in, rels_out, i_start);
}
inline void _TestBatchNewBufferElements(const Poly1d& gb, const Poly& rel, FnGetDeg f, int i_start)
{
    _BatchNewBufferElements(gb, rel, f, -1, i_start);
}
inline void _TestBatchNewBufferElementsV2(const Poly1d& gb, const Poly& rel, FnGetDeg f, int i_start)
{
    _BatchNewBufferElementsV2(gb, rel, f, -1, -1, i_start);
}
#endif
#endif

/**
 * Comsume relations from `buffer` in degree `<= deg`
 * while adding new relations back to `buffer` in degree `<= deg_max`.
 * `deg=-1` or `deg_max=-1` means infinity.
 */
template <typename GbType, typename BufferType, typename FnType> void AddRelsB(GbType& gb, BufferType& buffer, FnType _get_deg, int deg, int deg_max)
{
    auto p_buffer = buffer.begin();
    for (; p_buffer != buffer.end() && (deg == -1 || p_buffer->first <= deg); ++p_buffer) {
        /* Reduce relations from buffer in degree `p_buffer->first` */
        // std::cout << "t=" << p_buffer->first << '\n'; ////
        Poly1d rels;
#ifndef GROEBNER_MULTITHREAD /* Singlethreading */
        for (auto& poly : p_buffer->second) {
            Poly rel = [&]() {
                if constexpr (std::is_same<BufferType, GbBuffer>::value)
                    return Reduce(std::move(poly), gb);
                else
                    return Reduce(std::move(poly)->GetPoly(gb), gb);
            }();
            for (Poly& rel1 : rels)
                if (std::binary_search(rel.begin(), rel.end(), GbType::GetLead(rel1)))
                    rel += rel1;
            if (!rel.empty())
                rels.push_back(std::move(rel));
        }
#else /* Multithreading */
        Poly1d rels_tmp(p_buffer->second.size());
        std::vector<std::future<void>> futures;
        for (int i = 0; i < kNumThreads; ++i)
            futures.push_back(std::async(std::launch::async, [&gb, &p_buffer, &rels_tmp, i]() { return _BatchReduce(gb, p_buffer->second, rels_tmp, i); }));
        for (int i = 0; i < kNumThreads; ++i)
            futures[i].wait();
        for (auto& rel : rels_tmp) {
            for (Poly& rel1 : rels)
                if (std::binary_search(rel.begin(), rel.end(), GbType::GetLead(rel1)))
                    rel += rel1;
            if (!rel.empty())
                rels.push_back(std::move(rel));
        }
#endif

        /* Add these relations */
        for (auto& rel : rels) {
            if (!rel.empty()) {
#ifndef GROEBNER_MULTITHREAD /* Singlethreading */
                for (auto pg = gb.begin(); pg != gb.end(); ++pg) {
                    if (gcd_nonzero(GbType::GetLead(rel), GbType::GetLead(*pg))
                        && !(GbType::GetLead(rel)[0].gen < 0 && GbType::GetLead(*pg)[0].gen < 0 && GbType::GetLead(rel)[0].gen != GbType::GetLead(*pg)[0].gen)) {
                        if constexpr (std::is_same<BufferType, GbBuffer>::value) {
                            Mon lcm = LCM(GbType::GetLead(rel), GbType::GetLead(*pg));
                            int deg_new_rel = _get_deg(lcm);
                            if (deg_max == -1 || deg_new_rel <= deg_max) {
                                Poly new_rel = rel * div(lcm, GbType::GetLead(rel)) + (*pg) * div(lcm, GbType::GetLead(*pg));
                                if (!new_rel.empty())
                                    buffer[deg_new_rel].push_back(std::move(new_rel));
                            }
                        }
                        else { /* BufferType == GbBufferV2 */
                            Mon gcd = GCD(GbType::GetLead(rel), GbType::GetLead(*pg));
                            int deg_new_rel = p_buffer->first + _get_deg(GbType::GetLead(*pg)) - _get_deg(gcd);
                            if (deg_max == -1 || deg_new_rel <= deg_max)
                                buffer[deg_new_rel].push_back(std::make_unique<GcdBufferEle<GbType>>(std::move(gcd), (int)(pg - gb.begin()), (int)gb.size()));
                        }
                    }
                }
#else /* Multithreading */
                if constexpr (std::is_same<BufferType, GbBuffer>::value) { /* GbBuffer */
                    std::vector<std::future<std::vector<std::pair<int, Poly>>>> futures1;
                    for (int i = 0; i < kNumThreads; ++i)
                        futures1.push_back(std::async(std::launch::async, [&gb, &rel, _get_deg, deg_max, i]() { return _BatchNewBufferElements(gb, rel, _get_deg, deg_max, i); }));
                    for (int i = 0; i < kNumThreads; ++i) {
                        futures1[i].wait();
                        for (auto& [d, p] : futures1[i].get())
                            buffer[d].push_back(std::move(p));
                    }
                }
                else { /* GbBufferV2 */
                    std::vector<std::future<std::vector<std::pair<int, GcdBufferEle<GbType>>>>> futures1;
                    for (int i = 0; i < kNumThreads; ++i)
                        futures1.push_back(
                            std::async(std::launch::async, [&gb, &rel, _get_deg, &p_buffer, deg_max, i]() { return _BatchNewBufferElementsV2(gb, rel, _get_deg, p_buffer->first, deg_max, i); }));
                    for (int i = 0; i < kNumThreads; ++i) {
                        futures1[i].wait();
                        for (auto& [d, p] : futures1[i].get())
                            buffer[d].push_back(std::make_unique<GcdBufferEle<GbType>>(std::move(p.gcd_), p.i1_, p.i2_));
                    }
                }
#endif
                gb.push_back(std::move(rel));
            }
        }
    }
    buffer.erase(buffer.begin(), p_buffer);
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _TestAddRelsB(Poly1d& gb, GbBuffer& buffer, FnGetDeg f, int deg_max)
{
    AddRelsB(gb, buffer, f, -1, deg_max);
}
inline void _TestAddRelsB(Groebner& gb, GbBufferV2& buffer, FnGetDeg f, int deg_max)
{
    AddRelsB(gb, buffer, f, -1, deg_max);
}
#endif

/**
 * @see AddRelsB(GbType&, BufferType&, FnType, int, int)
 */
template <typename GbType, typename BufferType> void AddRelsB(GbType& gb, BufferType& buffer, const array& gen_degs, int deg, int deg_max)
{
    AddRelsB(gb, buffer, FnGetDeg{ gen_degs }, deg, deg_max);
}

/**
 * Convert `rels` into a groebner basis in degree `<= deg`
 * `deg_max=-1` means infinity.
 */
template <typename GbType, typename FnType, typename BufferType = GbBuffer> void AddRels(GbType& gb, Poly1d rels, FnType _get_deg, int deg_max)
{
    BufferType buffer;
    for (Poly& rel : rels) {
        if (!rel.empty()) {
            int deg = _get_deg(GbType::GetLead(rel));
            if constexpr (std::is_same<BufferType, typename GbType::BufferV2>::value)
                buffer[deg].push_back(std::make_unique<PolyBufferEle<GbType>>(std::move(rel)));
            else
                buffer[deg].push_back(std::move(rel));
        }
    }
    AddRelsB(gb, buffer, _get_deg, -1, deg_max);
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _TestAddRels(Poly1d& gb, Poly1d rels, FnGetDeg f, int deg_max)
{
    AddRels(gb, std::move(rels), f, deg_max);
}
inline void _TestAddRels(Groebner& gb, Poly1d rels, FnGetDeg f, int deg_max)
{
    AddRels<Groebner, FnGetDeg, GbBufferV2>(gb, std::move(rels), f, deg_max);
}
#endif

/**
 * @see AddRels(GbType&, Poly1d, FnType, int)
 */
template <typename GbType> void AddRels(GbType& gb, Poly1d rels, const array& gen_degs, int deg_max)
{
    AddRels(gb, std::move(rels), FnGetDeg{ gen_degs }, deg_max);
}
/**
 * @see AddRels(GbType&, Poly1d, FnType, int)
 */
template <typename GbType, typename FnType> void AddRelsV2(GbType& gb, Poly1d rels, FnType _get_deg, int deg_max)
{
    AddRels<GbType, FnType, typename GbType::BufferV2>(gb, std::move(rels), _get_deg, deg_max);
}
/**
 * @see AddRels(GbType&, Poly1d, FnType, int)
 */
template <typename GbType> void AddRelsV2(GbType& gb, Poly1d rels, const array& gen_degs, int deg_max)
{
    AddRels<GbType, FnGetDeg, typename GbType::BufferV2>(gb, std::move(rels), FnGetDeg{ gen_degs }, deg_max);
}

/**********************************************************
 * Algorithms that use Groebner basis
 **********************************************************/

/**
 * Compute the minimal generating set of `vectors` inplace.
 *
 * Each element of vectors is considered as an element of $R^n$ where $R$ is the algebra
 * determined by the Groebner basis `gb`.
 *
 * The degree of the basis of $R^n$ is determined by `basis_degs`.
 */
template <typename GbType> Poly2d& indecomposables(const GbType& gb, Poly2d& vectors, const array& gen_degs, const array& basis_degs)
{
    if (vectors.empty())
        return vectors;
    Groebner gb1 = gb;

    /* Convert each vector v into a relation \\sum vi x_{-i-1} */
    Poly1d rels;
    array degs;
    for (const Poly1d& v : vectors) {
        Poly rel;
        for (int i = 0; i < (int)basis_degs.size(); ++i)
            if (!v[i].empty())
                rel += v[i] * Mon{ { -i - 1, 1 } };
        degs.push_back(get_deg(GbType::GetLead(rel), gen_degs, basis_degs));
        rels.push_back(std::move(rel));
    }
    array indices = range((int)vectors.size());
    std::sort(indices.begin(), indices.end(), [&degs](int i, int j) { return degs[i] < degs[j]; });

    /* Add relations ordered by degree to gb1 */
    GbBuffer buffer1;
    int deg_max = degs[indices.back()];
    for (int i : indices) {
        AddRelsB(gb1, buffer1, FnGetDegV2{ gen_degs, basis_degs }, degs[i], deg_max);
        Poly rel = Reduce(rels[i], gb1);
        if (!rel.empty())
            buffer1[degs[i]].push_back(std::move(rel));
        else
            vectors[i].clear();
    }

    /* Keep only the indecomposables in `vectors` */
    RemoveEmptyElements(vectors);
    return vectors;
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _Test_indecomposables(const Groebner& gb, Poly2d& vectors, const array& gen_degs, const array& basis_degs)
{
    indecomposables(gb, vectors, gen_degs, basis_degs);
}
#endif

/**
 * Compute the generating set of linear relations among `polys`.
 *
 * The result is truncated by `deg<=deg_max`.
 */
template <typename GbType> Poly2d ann_seq(const GbType& gb, const Poly1d& polys, const array& gen_degs, int deg_max)
{
    Poly2d result;
    if (polys.empty())
        return result;
    Poly1d rels;
    array gen_degs1;
    int N = (int)polys.size();

    /* Add relations Xi=polys[i] to gb to obtain gb1 */
    for (int i = 0; i < N; ++i) {
        Poly p = polys[i];
        gen_degs1.push_back(get_deg(p, gen_degs));
        p.push_back({ { -i - 1, 1 } });
        rels.push_back(std::move(p));
    }
    Groebner gb1 = gb;
    AddRels(gb1, rels, FnGetDegV2{ gen_degs, gen_degs1 }, deg_max);

    /* Extract linear relations from gb1 */
    for (const Poly& g : gb1) {
        if (GbType::GetLead(g)[0].gen < 0) {
            Poly1d ann;
            ann.resize(N);
            for (const Mon& m : g) {
                MonInd p = m.begin();
                for (; p != m.end() && p->gen < 0; ++p)
                    ;
                Mon m1(m.begin(), p), m2(p, m.end());
                ann[size_t(-m1[0].gen) - 1] += Reduce(mul(subs(
                                                              { div(m1, { { m1[0].gen, 1 } }) }, [&polys](int i) { return polys[size_t(-i) - 1]; }, gb),
                                                          m2),
                                                      gb);
            }
            result.push_back(std::move(ann));
        }
    }

    /* Add commutators to linear relations */
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            if (deg_max == -1 || gen_degs1[i] + gen_degs1[j] <= deg_max) {
                Poly1d result_i;
                result_i.resize(N);
                result_i[i] = polys[j];
                result_i[j] = polys[i];
                result.push_back(std::move(result_i));
            }
        }
    }

    indecomposables(gb, result, gen_degs, gen_degs1);
    return result;
}
#ifdef GROEBNER_TEMPLATE_INSTANTIATIONS
inline void _Test_ann_seq(const Groebner& gb, const Poly1d& polys, const array& gen_degs, int deg_max)
{
    ann_seq(gb, polys, gen_degs, deg_max);
}
#endif

} /* namespace alg */

#endif /* GROEBNER_H */
