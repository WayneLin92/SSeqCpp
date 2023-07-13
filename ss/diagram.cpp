/* This module deals with operations that might affect all spectra */

#include "main.h"
#include "mylog.h"

/* Add a node */
void Diagram::AddNode(DeduceFlag flag)
{
    for (auto& ring : rings_)
        ring.nodes_ss.push_back({});
    for (auto& mod : modules_)
        mod.nodes_ss.push_back({});

    if (flag & DeduceFlag::pi) {
        for (auto& ring : rings_) {
            ring.pi_gb.AddNode();
            ring.nodes_pi_basis.push_back({});
        }
        for (auto& mod : modules_) {
            mod.pi_gb.AddNode();
            mod.nodes_pi_basis.push_back({});
        }

        //// TODO: Add pi_maps
    }
}

/* Pop the lastest node */
void Diagram::PopNode(DeduceFlag flag)
{
    for (auto& ring : rings_)
        ring.nodes_ss.pop_back();
    for (auto& mod : modules_)
        mod.nodes_ss.pop_back();
    if (flag & DeduceFlag::pi) {
        for (auto& ring : rings_) {
            ring.pi_gb.PopNode();
            ring.nodes_pi_basis.pop_back();
        }
        for (auto& mod : modules_) {
            mod.pi_gb.PopNode();
            mod.nodes_pi_basis.pop_back();
        }
        //// TODO: Add pi_maps
        // UpdateAllPossEinf();
    }
}

int1d MapRing2Ring::map(const int1d& x, AdamsDeg deg_x, const RingSp1d& rings) const
{

    int1d result;
    if (!x.empty()) {
        auto x_alg = Indices2Poly(x, rings[from].basis.at(deg_x));
        auto fx_alg = rings[to].gb.Reduce(subs(x_alg, images));
        if (fx_alg)
            result = Poly2Indices(fx_alg, rings[to].basis.at(deg_x));
    }
    return result;
}

int1d MapMod2Mod::map(const int1d& x, AdamsDeg deg_x, const ModSp1d& mods) const
{
    int1d result;
    if (!x.empty()) {
        auto x_alg = Indices2Mod(x, mods[from].basis.at(deg_x));
        auto fx_alg = mods[to].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + AdamsDeg(fil, fil - sus);
            result = Mod2Indices(fx_alg, mods[to].basis.at(deg_fx));
        }
    }
    return result;
}

auto subs(const Mod& x, const std::vector<Poly>& map1, const std::vector<Mod>& map2)
{
    Mod result{}, tmp{};
    Poly tmp_prod_p, tmf_p;
    for (const MMod& m : x.data) {
        Poly fm = Poly::Unit();
        for (auto p = m.m.begin(); p != m.m.end(); ++p) {
            powP(map1[p->g()], p->e_masked(), tmp_prod_p, tmf_p);
            fm.imulP(tmp_prod_p, tmf_p);
        }
        result.iaddP(fm * map2[m.v], tmp);
    }
    return result;
}

int1d MapMod2ModV2::map(const int1d& x, AdamsDeg deg_x, const ModSp1d& mods, const Map1d& maps) const
{
    int1d result;
    auto& ring_map = std::get<MapRing2Ring>(maps[over].map);
    if (!x.empty()) {
        auto x_alg = Indices2Mod(x, mods[from].basis.at(deg_x));
        auto fx_alg = mods[to].gb.Reduce(subs(x_alg, ring_map.images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + AdamsDeg(fil, fil - sus);
            result = Mod2Indices(fx_alg, mods[to].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMod2Ring::map(const int1d& x, AdamsDeg deg_x, const ModSp1d& mods, const RingSp1d& rings) const
{
    int1d result;
    if (!x.empty()) {
        auto x_alg = Indices2Mod(x, mods[from].basis.at(deg_x));
        auto fx_alg = rings[to].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + AdamsDeg(fil, fil - sus);
            result = Poly2Indices(fx_alg, rings[to].basis.at(deg_fx));
        }
    }
    return result;
}

int Diagram::SetRingDiffGlobal(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, DeduceFlag flag)
{
    int count = 0;
    auto& ring = rings_[iRing];
    auto& nodes_ss = ring.nodes_ss;
    int depth = int(nodes_ss.size() - 2);

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
        if (x.empty())
            count += SetRingBoundaryLeibniz(iRing, deg_dx, dx, r - 1);
        else {
            int r_min = LEVEL_MIN;
            while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
                ++r_min;
            if (dx.empty())
                r = NextRTgt(nodes_ss, ring.t_max, deg_x, r + 1) - 1;
            count += SetRingDiffLeibniz(iRing, deg_x, x, dx, r, r_min, flag);
        }

        for (size_t iMap : ring.ind_maps) {
            auto& map = maps_[iMap];
            if (deg_x.t <= map.t_max) {
                auto& f = std::get<MapRing2Ring>(map.map);
                int1d fdx;
                if (deg_dx.t <= map.t_max)
                    fdx = f.map(dx, deg_dx, rings_);
                else if (!dx.empty())
                    continue;
                auto fx = f.map(x, deg_x, rings_);
                if (IsNewDiff(rings_[f.to].nodes_ss, deg_x, fx, fdx, r)) {
                    Logger::LogDiff(depth, enumReason::nat, fmt::format("({}) {}", map.display, rings_[f.to].name), deg_x, fx, fdx, r);
                    count += SetRingDiffGlobal(f.to, deg_x, fx, fdx, r, true, flag);
                }
            }
        }

        /* x^2 are cycles in E_r */
        auto deg_xx = deg_x * 2;
        if (deg_xx.t <= ring.t_max && !x.empty() && dx.empty() && r < R_PERM - 1) {
            Poly poly_x = Indices2Poly(x, ring.basis.at(deg_x));
            Poly poly_xx;
            poly_x.frobP(poly_xx);
            poly_xx = ring.gb.Reduce(std::move(poly_xx));
            if (poly_xx) {
                int1d xx = Poly2Indices(poly_xx, ring.basis.at(deg_xx));
                if (IsNewDiff(nodes_ss, deg_xx, xx, {}, r + 1)) {
                    Logger::LogDiff(depth, enumReason::deduce_v2, ring.name, deg_xx, xx, {}, r + 1);
                    count += SetRingDiffGlobal(iRing, deg_xx, xx, {}, r + 1, true, flag);
                }
            }
        }
    }
    return count;
}

int Diagram::SetModuleDiffGlobal(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, DeduceFlag flag)
{
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    int count = 0;

    auto& mod = modules_[iMod];
    auto& ring = rings_[mod.iRing];
    auto& basis = mod.basis;
    auto& nodes_ss = mod.nodes_ss;
    int t_max = mod.t_max;

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        if (x.empty())
            count += SetModuleBoundaryLeibniz(iMod, deg_dx, dx, r - 1);
        else {
            int r_min = LEVEL_MIN;
            while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
                ++r_min;
            if (dx.empty())
                r = NextRTgt(nodes_ss, t_max, deg_x, r + 1) - 1;
            count += SetModuleDiffLeibniz(iMod, deg_x, x, dx, r, r_min, flag);
        }

        for (size_t iMap : mod.ind_maps) {
            auto& map = maps_[iMap];
            if (deg_x.t <= map.t_max) {
                if (std::holds_alternative<MapMod2Ring>(map.map)) {
                    auto& f = std::get<MapMod2Ring>(map.map);
                    int1d fdx;
                    if (deg_dx.t <= map.t_max)
                        fdx = f.map(dx, deg_dx, modules_, rings_);
                    else if (!dx.empty())
                        continue;
                    auto fx = f.map(x, deg_x, modules_, rings_);
                    AdamsDeg deg_fx = deg_x + AdamsDeg(f.fil, f.fil - f.sus);
                    if (IsNewDiff(rings_[f.to].nodes_ss, deg_fx, fx, fdx, r)) {
                        Logger::LogDiff(int(nodes_ss.size() - 2), enumReason::nat, fmt::format("({}) {}", map.display, rings_[f.to].name), deg_fx, fx, fdx, r);
                        count += SetRingDiffGlobal(f.to, deg_fx, fx, fdx, r, true, flag);
                    }
                }
                else if (std::holds_alternative<MapMod2Mod>(map.map)) {
                    auto& f = std::get<MapMod2Mod>(map.map);
                    int1d fdx;
                    if (deg_dx.t <= map.t_max)
                        fdx = f.map(dx, deg_dx, modules_);
                    else if (!dx.empty())
                        continue;
                    auto fx = f.map(x, deg_x, modules_);
                    AdamsDeg deg_fx = deg_x + AdamsDeg(f.fil, f.fil - f.sus);
                    if (IsNewDiff(modules_[f.to].nodes_ss, deg_fx, fx, fdx, r)) {
                        Logger::LogDiff(int(nodes_ss.size() - 2), enumReason::nat, fmt::format("({}) {}", map.display, modules_[f.to].name), deg_fx, fx, fdx, r);
                        count += SetModuleDiffGlobal(f.to, deg_fx, fx, fdx, r, true, flag);
                    }
                }
                else {
                    auto& f = std::get<MapMod2ModV2>(map.map);
                    int1d fdx;
                    if (deg_dx.t <= map.t_max)
                        fdx = f.map(dx, deg_dx, modules_, maps_);
                    else if (!dx.empty())
                        continue;
                    auto fx = f.map(x, deg_x, modules_, maps_);
                    AdamsDeg deg_fx = deg_x + AdamsDeg(f.fil, f.fil - f.sus);
                    if (IsNewDiff(modules_[f.to].nodes_ss, deg_fx, fx, fdx, r)) {
                        Logger::LogDiff(int(nodes_ss.size() - 2), enumReason::nat, fmt::format("({}) {}", map.display, modules_[f.to].name), deg_fx, fx, fdx, r);
                        count += SetModuleDiffGlobal(f.to, deg_fx, fx, fdx, r, true, flag);
                    }
                }
            }
        }
    }
    return count;
}

int Diagram::SetCwDiffGlobal(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, DeduceFlag flag)
{
    if (iCw < rings_.size())
        return SetRingDiffGlobal(iCw, deg_x, x, dx, r, newCertain, flag);
    else
        return SetModuleDiffGlobal(iCw - rings_.size(), deg_x, x, dx, r, newCertain, flag);
}
