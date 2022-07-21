#include "linalg.h"
#include "myexception.h"

#ifndef NDEBUG
#include <iostream>
#endif

namespace lina {

int2d GetSpace(const int2d& vectors)
{
    int2d result;
    for (int1d v : vectors) { /* Create copy on purpose */
        for (const auto& vi : result)
            if (std::binary_search(v.begin(), v.end(), vi[0]))
                v = AddVectors(v, vi);
        if (!v.empty())
            result.push_back(std::move(v));
    }
    return result;
}

int2d& SimplifySpace(int2d& spaceV)
{
    for (size_t i = spaceV.size() - 1; i != -1; i--)
        for (size_t j = 0; j < i; j++)
            if (std::binary_search(spaceV[j].begin(), spaceV[j].end(), spaceV[i][0]))
                spaceV[j] = AddVectors(spaceV[j], spaceV[i]);
    return spaceV;
}

int1d Residue(int2dIt spaceV_first, int2dIt spaceV_last, int1d v)
{
    for (auto p_vi = spaceV_first; p_vi != spaceV_last; ++p_vi)
        if (std::binary_search(v.begin(), v.end(), p_vi->front()))
            v = AddVectors(v, *p_vi);
    return v;
}

inline void AddToSpace(int2d& spaceV, const int1d& v)
{
    int1d v1 = Residue(spaceV, v);
    if (!v1.empty())
        spaceV.push_back(std::move(v1));
}

void GetInvMap(const int2d& fx, int2d& image, int2d& g)
{
    /* f(g[i]) = image[i] */
    for (size_t i = 0; i < fx.size(); ++i) {
        int1d src = {int(i)};
        int1d tgt = fx[i];
        for (size_t j = 0; j < image.size(); j++) {
            if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
                tgt = AddVectors(tgt, image[j]);
                src = AddVectors(src, g[j]);
            }
        }
        if (!tgt.empty()) {
            image.push_back(std::move(tgt));
            g.push_back(std::move(src));
        }
    }
}

void SetLinearMap(const int2d& fx, int2d& image, int2d& kernel, int2d& g)
{
    /* f(g[i]) = image[i] */
    for (size_t i = 0; i < fx.size(); ++i) {
        int1d src = {int(i)};
        int1d tgt = fx[i];
        for (size_t j = 0; j < image.size(); j++) {
            if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
                tgt = AddVectors(tgt, image[j]);
                src = AddVectors(src, g[j]);
            }
        }
        if (tgt.empty())
            AddToSpace(kernel, src);
        else {
            image.push_back(std::move(tgt));
            g.push_back(std::move(src));
        }
    }
}

void SetLinearMapV2(int2dIt x_first, int2dIt x_last, int2dIt fx_first, int2dIt fx_last, int2d& image, int2d& kernel, int2d& g)
{
    /* f(g[i]) = image[i] */
    for (auto p_x = x_first, p_fx = fx_first; p_x != x_last; ++p_x, ++p_fx) {
        int1d src = *p_x;
        int1d tgt = *p_fx;
        for (size_t j = 0; j < image.size(); j++) {
            if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
                tgt = AddVectors(tgt, image[j]);
                src = AddVectors(src, g[j]);
            }
        }
        if (tgt.empty())
            AddToSpace(kernel, src);
        else {
            image.push_back(std::move(tgt));
            g.push_back(std::move(src));
        }
    }
}

void SetLinearMapV2(const int1d& x, const int2d& fx, int2d& image, int2d& kernel, int2d& g)
{
    /* f(g[i]) = image[i] */
    for (size_t i = 0; i < fx.size(); ++i) {
        int1d src = {x[i]};
        int1d tgt = fx[i];
        for (size_t j = 0; j < image.size(); j++) {
            if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
                tgt = AddVectors(tgt, image[j]);
                src = AddVectors(src, g[j]);
            }
        }
        if (tgt.empty())
            AddToSpace(kernel, src);
        else {
            image.push_back(std::move(tgt));
            g.push_back(std::move(src));
        }
    }
}

void SetLinearMapV3(const int2d& x, const int2d& fx, int2d& domain, int2d& f, int2d& image, int2d& g, int2d& kernel)
{
    /* f(g[i]) = image[i] */
    for (size_t i = 0; i < fx.size(); ++i) {
        int1d src = x[i];
        int1d tgt = fx[i];
        for (size_t k = 0; k < domain.size(); ++k) {
            if (std::binary_search(src.begin(), src.end(), domain[k][0])) {
                src = AddVectors(src, domain[k]);
                tgt = AddVectors(tgt, f[k]);
            }
        }
        if (src.empty()) {
            if (!tgt.empty())
                throw MyException(0x67b8b67dU, "conflicting linear map definition");
        }
        else {
            domain.push_back(src);
            f.push_back(tgt);
            for (size_t j = 0; j < image.size(); j++) {
                if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
                    tgt = AddVectors(tgt, image[j]);
                    src = AddVectors(src, g[j]);
                }
            }
            if (tgt.empty())
                AddToSpace(kernel, src);
            else {
                image.push_back(std::move(tgt));
                g.push_back(std::move(src));
            }
        }
    }
}

int1d GetImage(int2dIt spaceV_first, int2dIt spaceV_last, int2dIt f_first, int2dIt f_last, int1d v)
{
    int1d result;
    for (auto p_Vi = spaceV_first, p_fi = f_first; p_Vi != spaceV_last && !v.empty(); ++p_Vi, ++p_fi)
        if (std::binary_search(v.begin(), v.end(), p_Vi->front())) {
            v = AddVectors(v, *p_Vi);
            result = AddVectors(result, *p_fi);
        }
#ifndef NDEBUG
    if (!v.empty()) {
        throw MyException(0x6a4fe8a1U, "v is not in spaceV");
    }
#endif
    return result;
}

int1d GetInvImage(const int2d& spaceV, int1d w)
{
    int1d result;
    for (size_t j = 0; j < spaceV.size(); j++) {
        if (std::binary_search(w.begin(), w.end(), spaceV[j][0])) {
            w = AddVectors(w, spaceV[j]);
            result.push_back((int)j);
        }
    }
#ifndef NDEBUG
    if (!w.empty())
        throw MyException(0x1a4ef6d8U, "GetInvImage not well defined");
#endif
    return result;
}

int2d QuotientSpace(const int2d& spaceV, const int2d& spaceW)
{
    int2d quotient;
    size_t dimQuo = spaceV.size() - spaceW.size();
#ifndef NDEBUG
    for (size_t i = 0; i < spaceV.size(); i++)
#else
    for (size_t i = 0; i < spaceV.size() && quotient.size() < dimQuo; i++)
#endif
    {
        auto v1 = Residue(quotient, Residue(spaceW, spaceV[i]));
        if (!v1.empty())
            quotient.push_back(std::move(v1));
    }
#ifndef NDEBUG
    if (quotient.size() != dimQuo) {
        std::cerr << "W is not a subspace of V!\n";
        throw "cec7f701";
    }
#endif
    return quotient;
}

} /* namespace lina */