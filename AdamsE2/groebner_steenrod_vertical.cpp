#include "groebner_steenrod.h"

namespace steenrod {
void AddRels(Groebner& gb, const Mod1d& rels, int deg, array& min_gb, Mod1d& redundent_gb)
{
    int deg_max = gb.deg_trunc();
    if (deg > deg_max)
        throw MyException(0xb2474e19U, "deg is bigger than the truncation degree.");

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, array> rels_graded;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            int d = rels[i].GetLead().deg(gb.basis_degs());
            if (d <= deg)
                rels_graded[d].push_back((int)i);
        }
    }
    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg && (d <= deg_max_rels || !gb.gb_pairs().empty_pairs_for_gb()); ++d) {
        std::cout << "t=" << d << "            \r";
        int size_pairs_d;
        PCriticalPair1d pairs_d;
        if (gb.gb_pairs().empty_pairs_for_gb_d(d))
            size_pairs_d = 0;
        else {
            gb.AddPairsAndMinimize(d);
            pairs_d = gb.pairs(d);
            size_pairs_d = (int)pairs_d.size();
        }
        auto p_rels_d = rels_graded.find(d);
        int size_rels_d = p_rels_d != rels_graded.end() ? (int)p_rels_d->second.size() : 0;
        if (size_pairs_d + size_rels_d == 0)
            continue;
        Mod1d rels_tmp(size_pairs_d + size_rels_d);

        /* Reduce relations in degree d */
        if (size_pairs_d) {
            ut::Range range(0, size_pairs_d);
            std::for_each(std::execution::seq, range.begin(), range.end(), [&gb, &pairs_d, &rels_tmp](int i) { rels_tmp[i] = gb.Reduce(*pairs_d[i]); });
        }
        if (size_rels_d) {
            auto& rels_d = p_rels_d->second;
            ut::Range range(0, size_rels_d);
            std::for_each(range.begin(), range.end(), [&gb, &rels, &rels_d, &rels_tmp, size_pairs_d](int i) { rels_tmp[size_pairs_d + i] = gb.Reduce(rels[rels_d[i]]); });
        }

        /* Triangulate these relations */
        Mod1d rels_d;
        for (size_t i = 0; i < rels_tmp.size(); ++i) {
            for (size_t j = 0; j < rels_d.size(); ++j) {
                if (std::binary_search(rels_tmp[i].data.begin(), rels_tmp[i].data.end(), rels_d[j].GetLead())) {
                    rels_tmp[i] += rels_d[j];
                    if (i < size_pairs_d)
                        pairs_d[i]->x += MMod{MMilnor{0}, int(gb.size() + j)};
                }
            }
            if (rels_tmp[i]) {
                if (i < size_pairs_d) {
                    redundent_gb.push_back(pairs_d[i]->x);
                    pairs_d[i]->x += MMod{MMilnor{0}, int(gb.size() + rels_d.size())};
                }
                else
                    min_gb.push_back(int(gb.size() + rels_d.size()));
                rels_d.push_back(std::move(rels_tmp[i]));
            }
        }

        /* Add these relations */
        for (auto& rel : rels_d) {
            gb.push_back(std::move(rel));
        }
    }
}

void AddRelsMinRes(GroebnerMinRes& gb, const Mod1d& rels, int deg)
{
    int deg_max = gb.deg_trunc();
    if (deg > deg_max)
        throw MyException(0xb2474e19U, "deg is bigger than the truncation degree.");

    /* Calculate the degrees of `rels` and group them by degree */
    std::map<int, array> rels_graded;
    for (size_t i = 0; i < rels.size(); ++i) {
        if (rels[i]) {
            auto& lead = rels[i].GetLead();
            int d = lead.m.deg() + GetBaseDeg(gb.basis_degs(), lead.v);
            if (d <= deg)
                rels_graded[d].push_back((int)i);
        }
    }
    int deg_max_rels = rels_graded.empty() ? 0 : rels_graded.rbegin()->first;
    for (int d = 1; d <= deg && (d <= deg_max_rels || !gb.gb_pairs().empty_pairs_for_gb()); ++d) {
        std::cout << "t=" << d << "                   \n";

        PCriticalPairMinRes2d pairs_d;
        if (!gb.gb_pairs().empty_pairs_for_gb_d(d)) {
            gb.AddPairsAndMinimize(d);
            pairs_d = gb.pairs(d);
        }
        for (int s = (int)pairs_d.size(); s-- > -1;) {
            std::cout << "  s=" << s << "            \r";
            Mod1d rels_tmp;
            if (s == -1) {
                auto p_rels_d = rels_graded.find(d);
                if (p_rels_d != rels_graded.end()) {
                    rels_tmp.resize(p_rels_d->second.size());
                    auto& rels_d = p_rels_d->second;
                    ut::Range range(0, (int)rels_tmp.size());
                    std::for_each(range.begin(), range.end(), [&](int i) { rels_tmp[i] = gb.HeadReduce(rels[rels_d[i]], 0); });
                }
            }
            else {
                if (s < pairs_d.size() && pairs_d[s].size()) {
                    rels_tmp.resize(pairs_d[s].size());
                    ut::Range range(0, (int)rels_tmp.size());
                    std::for_each(std::execution::seq, range.begin(), range.end(), [&](int i) { rels_tmp[i] = gb.HeadReduce(*pairs_d[s][i]); });
                }
            }
            if (rels_tmp.empty())
                continue;

            /* Triangulate these relations */
            Mod1d rels_d;
            Mod1d kernel_splus_tmp;
            for (size_t i = 0; i < rels_tmp.size(); ++i) {
                for (size_t j = 0; j < rels_d.size(); ++j)
                    if (std::binary_search(rels_tmp[i].data.begin(), rels_tmp[i].data.end(), rels_d[j].GetLead()))
                        rels_tmp[i] += rels_d[j];
                if (rels_tmp[i]) {
                    if (GetS(rels_tmp[i].GetLead().v) == s)
                        rels_d.push_back(std::move(rels_tmp[i]));
                    else {
                        kernel_splus_tmp.push_back(gb.HeadReduce(std::move(rels_tmp[i]), s + 1));
                    }
                }
            }
            for (auto& k : kernel_splus_tmp)
                ut::RemoveIf(k.data, [s](MMod x) { return GetS(x.v) > s + 1; });
            Mod1d kernel_splus;
            for (size_t i = 0; i < kernel_splus_tmp.size(); ++i) {
                for (size_t j = 0; j < kernel_splus.size(); ++j)
                    if (std::binary_search(kernel_splus_tmp[i].data.begin(), kernel_splus_tmp[i].data.end(), kernel_splus[j].GetLead()))
                        kernel_splus_tmp[i] += kernel_splus[j];
                if (kernel_splus_tmp[i])
                    kernel_splus.push_back(std::move(kernel_splus_tmp[i]));
            }

            /* Add these relations */
            for (auto& rel : rels_d)
                gb.push_back(std::move(rel));
            for (auto& rel : kernel_splus) {
#ifndef NDEBUG
                if (rel && GetS(rel.GetLead().v) != s + 1) {
                    std::cout << "rel=";
                    print(rel);
                    std::cout << "\n\n";
                    throw MyException(0, "BUG");
                }
                if (s != -1 && TplSubs(rel, [&](int v) { return gb.kernel()[s][v & 0xffff]; })) {
                    std::cout << "rel=";
                    print(rel);
                    std::cout << "\n\n";
                    /*for (auto& g : gb.data()) {
                        print(g);
                        std::cout << '\n';
                    }*/
                    throw MyException(0, "BUG");
                }
#endif
                /*std::cout << "rel=";
                print(rel);
                std::cout << '\n';
                if (d == 5 && s == 0)
                    std::cout << "test\n";*/
                gb.push_back_kernel(std::move(rel), d, (int)s + 1);
                /*print(gb.data().back());
                std::cout << '\n';*/
            }

            if (size_t size_k = kernel_splus.size())
                std::cout << "  s=" << s + 2 << " dim=" << size_k << '\n';
        }
    }
}

}  // namespace steenrod