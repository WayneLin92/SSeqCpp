/* This module deals with operations that might affect all spectra */

#include "main.h"
#include "mylog.h"

/* Add a node */
void Diagram::AddNode(SSFlag flag)
{
    for (auto& ring : rings_)
        ring.nodes_ss.push_back({});
    for (auto& mod : modules_)
        mod.nodes_ss.push_back({});
    if (flag & SSFlag::cofseq) {
        for (auto& cofseq : cofseqs_)
            for (size_t iCs = 0; iCs < 3; ++iCs)
                cofseq.nodes_cofseq[iCs].push_back({});
    }

    if (flag & SSFlag::pi) {
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
void Diagram::PopNode(SSFlag flag)
{
    for (auto& ring : rings_)
        ring.nodes_ss.pop_back();
    for (auto& mod : modules_)
        mod.nodes_ss.pop_back();
    if (flag & SSFlag::cofseq) {
        for (auto& cofseq : cofseqs_)
            for (size_t iCs = 0; iCs < 3; ++iCs)
                cofseq.nodes_cofseq[iCs].pop_back();
    }

    if (flag & SSFlag::pi) {
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

int1d MapRing2Ring::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto x_alg = Indices2Poly(x, rings[from].basis.at(deg_x));
        auto fx_alg = rings[to].gb.Reduce(subs(x_alg, images));
        if (fx_alg)
            result = Poly2Indices(fx_alg, rings[to].basis.at(deg_x));
    }
    return result;
}

int1d MapMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[from].basis.at(deg_x));
        auto fx_alg = mods[to].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Mod::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& gen_degs = ring_gen_degs[mods[from].iRing];
    for (auto& [deg_x, basis_d] : mods[from].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to].basis.at(deg_x + deg));
                Mod alg_gfx = mods[to].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Mod2Indices(alg_gfx, mods[to].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != gfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

auto subs(const Mod& x, const Poly1d& map1, const Mod1d& map2)
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

int1d MapMod2ModV2::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto& maps = diagram.GetMaps();
        auto& ring_map = maps[over];
        auto x_alg = Indices2Mod(x, mods[from].basis.at(deg_x));
        auto fx_alg = mods[to].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2ModV2::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& maps = diagram.GetMaps();
    auto& gen_degs = ring_gen_degs[mods[from].iRing];
    for (auto& [deg_x, basis_d] : mods[from].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to].basis.at(deg_x + deg));
                Mod alg_fgfx = mods[to].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Mod2Indices(alg_fgfx, mods[to].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != fgfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

int1d MapMod2Ring::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[from].basis.at(deg_x));
        auto fx_alg = rings[to].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Ring::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& rings = diagram.GetRings();
    auto& gen_degs = ring_gen_degs[mods[from].iRing];
    for (auto& [deg_x, basis_d] : mods[from].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to].basis.at(deg_x + deg));
                auto alg_gfx = rings[to].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Poly2Indices(alg_gfx, rings[to].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != gfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

auto subs(const Mod& x, const Poly1d& map1, const Poly1d& map2)
{
    Poly result, tmp_prod, tmp;
    for (const MMod& m : x.data) {
        Poly fm = Poly::Unit();
        for (auto p = m.m.begin(); p != m.m.end(); ++p) {
            powP(map1[p->g()], p->e_masked(), tmp_prod, tmp);
            fm.imulP(tmp_prod, tmp);
        }
        result.iaddP(fm * map2[m.v], tmp);
    }
    return result;
}

int1d MapMod2RingV2::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto& mods = diagram.GetModules();
        auto& maps = diagram.GetMaps();
        auto& ring_map = maps[over];
        auto x_alg = Indices2Mod(x, mods[from].basis.at(deg_x));
        auto fx_alg = rings[to].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2RingV2::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& rings = diagram.GetRings();
    auto& maps = diagram.GetMaps();
    auto& gen_degs = ring_gen_degs[mods[from].iRing];
    for (auto& [deg_x, basis_d] : mods[from].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to].basis.at(deg_x + deg));
                auto alg_fgfx = rings[to].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Poly2Indices(alg_fgfx, rings[to].basis.at(deg_gx + deg)) : int1d{};
                if (fgx != fgfx) {
                    fmt::print("Incorrect map: {} deg_x={}, x={}, deg_g={}, g={}\n", name, deg_x, i, gen_degs[g], poly_g.Str());
                    throw MyException(0x18a1700f, "Incorrect map");
                }
            }
        }
    }
}

int1d MapMulRing2Ring::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto x_alg = Indices2Poly(x, rings[index].basis.at(deg_x));
        auto fx_alg = rings[index].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[index].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMulRing2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& rings = diagram.GetRings();
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Poly(x, rings[from].basis.at(deg_x));
        auto fx_alg = mods[to].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMulMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[index].basis.at(deg_x));
        auto fx_alg = mods[index].gb.Reduce(factor * x_alg);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[index].basis.at(deg_fx));
        }
    }
    return result;
}

int Diagram::SetRingDiffGlobal(size_t iRing, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
{
    int count = 0;
    auto& ring = rings_[iRing];
    auto& nodes_ss = ring.nodes_ss;

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        int depth = int(nodes_ss.size() - 2);
        AdamsDeg deg_dx = deg_x + AdamsDeg{r, r - 1};
        if (x.empty())
            count += SetRingBoundaryLeibniz(iRing, deg_dx, dx, r - 1, flag);
        else {
            int r_min = LEVEL_MIN;
            while (r_min < r && !IsNewDiff(nodes_ss, deg_x, x, {}, r_min))  // TODO: improve this
                ++r_min;
            if (dx.empty())
                r = NextRTgt(nodes_ss, ring.t_max, deg_x, r + 1) - 1;
            count += SetRingDiffLeibniz(iRing, deg_x, x, dx, r, r_min, flag);
        }

        for (size_t iMap : ring.ind_maps) {
            auto map = (MapRing2Ring*)maps_[iMap].get();
            if (deg_x.t <= map->t_max) {
                int1d fdx;
                if (deg_dx.t <= map->t_max)
                    fdx = map->map(dx, deg_dx, *this);
                else if (!dx.empty())
                    continue;
                auto fx = map->map(x, deg_x, *this);
                if (IsNewDiff(rings_[map->to].nodes_ss, deg_x, fx, fdx, r)) {
                    Logger::LogDiff(depth, EnumReason::nat, fmt::format("({}) {}", map->display, rings_[map->to].name), deg_x, fx, fdx, r);
                    count += SetRingDiffGlobal(map->to, deg_x, fx, fdx, r, true, flag);
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
                    Logger::LogDiff(depth, EnumReason::deduce_xx, ring.name, deg_xx, xx, {}, r + 1);
                    count += SetRingDiffGlobal(iRing, deg_xx, xx, {}, r + 1, true, flag);
                }
            }
        }
    }
    return count;
}

int Diagram::SetModuleDiffGlobal(size_t iMod, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
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
            count += SetModuleBoundaryLeibniz(iMod, deg_dx, dx, r - 1, flag);
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
            if (deg_x.t <= map->t_max) {
                size_t to;
                if (map->IsToRing(to)) {
                    int1d fdx;
                    if (deg_dx.t <= map->t_max)
                        fdx = map->map(dx, deg_dx, *this);
                    else if (!dx.empty())
                        continue;
                    auto fx = map->map(x, deg_x, *this);
                    AdamsDeg deg_fx = deg_x + map->deg;
                    if (IsNewDiff(rings_[to].nodes_ss, deg_fx, fx, fdx, r)) {
                        Logger::LogDiff(int(nodes_ss.size() - 2), EnumReason::nat, fmt::format("({}) {}", map->display, rings_[to].name), deg_fx, fx, fdx, r);
                        count += SetRingDiffGlobal(to, deg_fx, fx, fdx, r, true, flag);
                    }
                }
                else {
                    int1d fdx;
                    if (deg_dx.t <= map->t_max)
                        fdx = map->map(dx, deg_dx, *this);
                    else if (!dx.empty())
                        continue;
                    auto fx = map->map(x, deg_x, *this);
                    AdamsDeg deg_fx = deg_x + map->deg;
                    if (IsNewDiff(modules_[to].nodes_ss, deg_fx, fx, fdx, r)) {
                        Logger::LogDiff(int(nodes_ss.size() - 2), EnumReason::nat, fmt::format("({}) {}", map->display, modules_[to].name), deg_fx, fx, fdx, r);
                        count += SetModuleDiffGlobal(to, deg_fx, fx, fdx, r, true, flag);
                    }
                }
            }
        }
    }
    return count;
}

int Diagram::SetCwDiffGlobal(size_t iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)
{
    if (iCw < rings_.size())
        return SetRingDiffGlobal(iCw, deg_x, x, dx, r, newCertain, flag);
    else
        return SetModuleDiffGlobal(iCw - rings_.size(), deg_x, x, dx, r, newCertain, flag);
}
