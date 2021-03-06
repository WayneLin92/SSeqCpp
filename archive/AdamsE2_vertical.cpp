#include "algebras/benchmark.h"
#include "algebras/dbalg.h"
#include "algebras/utility.h"
#include "groebner_steenrod.h"

void AdamsE2()
{
    using namespace steenrod;

    Timer timer;

    int t_trunc = 70;
    int s_trunc = t_trunc / 3 + 1;
    Mod1d rels;
    for (int i = 0; i < 10; ++i) {
        rels.push_back(MMod{MMilnor::P(i, i + 1), 0});
        if (MMilnor::P(i, i + 1).deg() > t_trunc)
            break;
    }
    array basis_degrees = {0};

    for (int s = 1; s <= s_trunc; ++s) {
        std::cout << "s=" << s << '\n';
        array min_gb;
        Mod1d red_gb;
        Groebner gb(t_trunc, basis_degrees);
        AddRels(gb, rels, t_trunc, min_gb, red_gb);

        std::cout << "dim=" << min_gb.size() << "\n\n";

        Mod1d map(gb.size());
        for (size_t i = 0, j = 0; i < min_gb.size(); ++i)
            map[min_gb[i]] = MMod{MMilnor{0}, (int)j++};
        for (size_t i = 0, j = 0; i < map.size(); ++i)
            if (map[i].data.empty())
                map[i] = subs(red_gb[j++], map);

        rels.clear();
        for (Mod syz : gb.MinSyzOfGb()) {
            Mod syz1 = subs(syz, map);
            if (syz1)
                rels.push_back(syz1);
        }

        basis_degrees.clear();
        for (auto i : min_gb)
            basis_degrees.push_back(gb[i].GetLead().deg(gb.basis_degs()));

        /*auto gb_data = gb.data();
        for (size_t i = 0; i < gb_data.size(); ++i)
            std::cout << "g_" << i << '=' << gb_data[i] << '\n';

        std::cout << "\nSyzygies:\n";
        for (auto& syz : syzygies)
            std::cout << syz << '\n';*/
    }
}

void AdamsE2V2()
{
    using namespace steenrod;

    Timer timer;

    int t_trunc = 50;
    Mod1d rels;
    for (int i = 0; i < 10; ++i) {
        rels.push_back(MMod{MMilnor::P(i, i + 1), GetV(0, 0)});
        if (MMilnor::P(i, i + 1).deg() > t_trunc)
            break;
    }
    array2d basis_degrees = {{0}};
    GroebnerMinRes gb(t_trunc, basis_degrees);
    AddRelsMinRes(gb, rels, t_trunc);
    // std::cout << "gb:\n";
    // for (auto& g : gb.data()) {////
    //     std::cout << "t=" << g.GetLead().m.deg() + GetBaseDeg(gb.basis_degs(), g.GetLead().v) << ' ';
    //     print(g);
    //     std::cout << '\n';
    // }
}

int main()
{
    AdamsE2V2();
    return 0;
}