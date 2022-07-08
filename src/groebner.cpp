#include "groebner.h"
#include "myio.h"
#include <iterator>

namespace alg {

bool detail::HasGCD(const Mon& mon1, const Mon& mon2)
{
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw())
            k++;
        else if (k->g_raw() < l->g_raw())
            l++;
        else
            return true;
    }
    return false;
}

int detail::DegLCM(const Mon& mon1, const Mon& mon2, const int1d& gen_degs)
{
    int result = 0;
    auto k = mon1.begin(), l = mon2.begin();
    while (k != mon1.end() && l != mon2.end()) {
        if (k->g_raw() > l->g_raw()) {
            result += gen_degs[k->g()] * k->e();
            ++k;
        }
        else if (k->g_raw() < l->g_raw()) {
            result += gen_degs[l->g()] * l->e();
            ++l;
        }
        else {
            if (k->e() < l->e())
                result += gen_degs[l->g()] * l->e();
            else
                result += gen_degs[k->g()] * k->e();
            ++k;
            ++l;
        }
    }
    for (; k != mon1.end(); ++k)
        result += gen_degs[k->g()] * k->e();
    for (; l != mon2.end(); ++l)
        result += gen_degs[l->g()] * l->e();
    return result;
}

void detail::ReduceP(Poly& f, const Poly& g, const size_t index, Mon& tmp_q, Mon& tmp_prod)
{
    Mon1d h;
    h.reserve(f.data.size() + g.data.size() - (size_t)index - 2);
    divP(f.data[index], g.data[0], tmp_q);
    size_t i = index + 1, j = 1;
    bool updateProd = true;
    while (j < g.data.size()) {
        if (i == f.data.size()) {
            f.data.resize(index + h.size() + g.data.size() - j);
            for (size_t k = 0; k < h.size(); ++k)
                f.data[index + k] = std::move(h[k]);
            size_t index1 = index + h.size();
            for (size_t k = j; k < g.data.size(); ++k)
                mulP(g.data[k], tmp_q, f.data[index1 + k - j]);
            return;
        }
        if (updateProd)
            mulP(g.data[j], tmp_q, tmp_prod);
        if (f.data[i] < tmp_prod) {
            updateProd = false;
            h.push_back(std::move(f.data[i++]));
        }
        else {
            updateProd = true;
            ++j;
            if (tmp_prod < f.data[i])
                h.push_back(Mon(tmp_prod));
            else
                ++i;
        }
    }
    size_t index1 = index + h.size();
    size_t new_size = index1 + f.data.size() - i;
    if (new_size > f.data.size()) {
        size_t old_size = f.data.size();
        f.data.resize(new_size);
        for (size_t k = old_size; k-- > i;)
            f.data[index1 + k - i] = std::move(f.data[k]);
    }
    else if (new_size < f.data.size()) {
        for (size_t k = i; k < f.data.size(); ++k)
            f.data[index1 + k - i] = std::move(f.data[k]);
        f.data.erase(f.data.begin() + new_size, f.data.end());
    }
    for (size_t k = 0; k < h.size(); ++k)
        f.data[index + k] = std::move(h[k]);
}

void CriPair::SetFromLM(CriPair& result, const Mon& lead1, const Mon& lead2, int i, int j)
{
    auto k = lead1.begin(), l = lead2.begin();
    while (k != lead1.end() && l != lead2.end()) {
        if (k->g_raw() > l->g_raw())
            result.m2.push_back(*k++);
        else if (k->g_raw() < l->g_raw())
            result.m1.push_back(*l++);
        else {
            if (k->e() < l->e())
                result.m1.push_back(GE(l->data - k->e()));
            else if (k->e() > l->e())
                result.m2.push_back(GE(k->data - l->e()));
            k++;
            l++;
        }
    }
    if (k != lead1.end())
        result.m2.data.insert(result.m2.end(), k, lead1.end());
    else
        result.m1.data.insert(result.m1.end(), l, lead2.end());
    result.i1 = i;
    result.i2 = j;
    result.trace_m2 = result.m2.Trace();
}

} /* namespace alg */