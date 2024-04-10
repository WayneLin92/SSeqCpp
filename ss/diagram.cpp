/* This module deals with operations that might affect all spectra */

#include "algebras/linalg.h"
#include "main.h"
#include "mylog.h"

/* Add a node */
void Diagram::AddNode(SSFlag flag)
{
    ++depth_;
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
    --depth_;
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
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg)
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_x));
    }
    return result;
}

int1d MapMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Mod::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to.index].basis.at(deg_x + deg));
                Mod alg_gfx = mods[to.index].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Mod2Indices(alg_gfx, mods[to.index].basis.at(deg_gx + deg)) : int1d{};
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
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2ModV2::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& maps = diagram.GetMaps();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                Mod alg_fx = fx.empty() ? Mod() : Indices2Mod(fx, mods[to.index].basis.at(deg_x + deg));
                Mod alg_fgfx = mods[to.index].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Mod2Indices(alg_fgfx, mods[to.index].basis.at(deg_gx + deg)) : int1d{};
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
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

void MapMod2Ring::Verify(const Diagram& diagram, const AdamsDeg2d& ring_gen_degs)
{
    fmt::print("Verifying {}\n", name);
    auto& mods = diagram.GetModules();
    auto& rings = diagram.GetRings();
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to.index].basis.at(deg_x + deg));
                auto alg_gfx = rings[to.index].gb.Reduce(poly_g * alg_fx);
                int1d gfx = alg_gfx ? Poly2Indices(alg_gfx, rings[to.index].basis.at(deg_gx + deg)) : int1d{};
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
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(subs(x_alg, ((MapRing2Ring*)maps[over].get())->images, images));
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
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
    auto& gen_degs = ring_gen_degs[mods[from.index].iRing];
    for (auto& [deg_x, basis_d] : mods[from.index].basis) {
        for (size_t g = 0; g < gen_degs.size(); ++g) {
            AdamsDeg deg_gx = deg_x + gen_degs[g];
            if (deg_gx.t > t_max)
                break;
            Poly poly_g = Poly::Gen((uint32_t)g);
            Poly poly_fg = ((MapRing2Ring*)maps[over].get())->images[g];
            for (size_t i = 0; i < basis_d.size(); ++i) {
                Mod alg_gx = mods[from.index].gb.Reduce(poly_g * basis_d[i]);
                int1d gx = alg_gx ? Mod2Indices(alg_gx, mods[from.index].basis.at(deg_gx)) : int1d{};
                int1d fgx = map(gx, deg_gx, diagram);
                int1d fx = map(int1d{(int)i}, deg_x, diagram);
                auto alg_fx = fx.empty() ? Poly() : Indices2Poly(fx, rings[to.index].basis.at(deg_x + deg));
                auto alg_fgfx = rings[to.index].gb.Reduce(poly_fg * alg_fx);
                int1d fgfx = alg_fgfx ? Poly2Indices(alg_fgfx, rings[to.index].basis.at(deg_gx + deg)) : int1d{};
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
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = rings[to.index].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Poly2Indices(fx_alg, rings[to.index].basis.at(deg_fx));
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
        auto x_alg = Indices2Poly(x, rings[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(x_alg * factor);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
        }
    }
    return result;
}

int1d MapMulMod2Mod::map(const int1d& x, AdamsDeg deg_x, const Diagram& diagram) const
{
    int1d result;
    if (!x.empty()) {
        auto& mods = diagram.GetModules();
        auto x_alg = Indices2Mod(x, mods[from.index].basis.at(deg_x));
        auto fx_alg = mods[to.index].gb.Reduce(factor * x_alg);
        if (fx_alg) {
            AdamsDeg deg_fx = deg_x + deg;
            result = Mod2Indices(fx_alg, mods[to.index].basis.at(deg_fx));
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
        int depth = depth_;
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
                if (IsNewDiff(rings_[map->to.index].nodes_ss, deg_x, fx, fdx, r)) {
                    Logger::LogDiff(depth, EnumReason::nat, fmt::format("({}) {}", map->display, rings_[map->to.index].name), deg_x, fx, fdx, r);
                    count += SetRingDiffGlobal(map->to.index, deg_x, fx, fdx, r, true, flag);
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
    int count = 0;
    auto& mod = modules_[iMod];
    auto& nodes_ss = mod.nodes_ss;
    int t_max = mod.t_max;

    if (newCertain || IsNewDiff(nodes_ss, deg_x, x, dx, r)) {
        AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
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
                int1d fdx;
                if (deg_dx.t <= map->t_max)
                    fdx = map->map(dx, deg_dx, *this);
                else if (!dx.empty())
                    continue;
                auto fx = map->map(x, deg_x, *this);
                AdamsDeg deg_fx = deg_x + map->deg;
                if (IsNewDiff(GetSS(map->to), deg_fx, fx, fdx, r)) {
                    Logger::LogDiff(depth_, EnumReason::nat, fmt::format("({}) {}", map->display, GetCwName(map->to)), deg_fx, fx, fdx, r);
                    count += SetCwDiffGlobal(map->to, deg_fx, fx, fdx, r, true, flag);
                }
            }
        }
    }
    return count;
}

int Diagram::SetCwDiffGlobal(IndexCw iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool newCertain, SSFlag flag)  //// TODO: write two functions into one
{
    int count = 0;
    if (iCw.isRing)
        count += SetRingDiffGlobal(iCw.index, deg_x, x, dx, r, newCertain, flag);
    else
        count += SetModuleDiffGlobal(iCw.index, deg_x, x, dx, r, newCertain, flag);

    if (depth_ <= 1 && (flag & SSFlag::synthetic) && r <= 6 && x.size() && dx != NULL_DIFF && dx.size()) {
        bool has_cross = GetCrossR(GetSS(iCw), deg_x, GetTMax(iCw)) <= r;
        SetCwDiffSynthetic(iCw, deg_x, x, dx, r, has_cross, flag);
    }
    return count;
}

int Diagram::SetCwDiffSynthetic(IndexCw iCw, AdamsDeg deg_x, const int1d& x, const int1d& dx, int r, bool hasCross, SSFlag flag)
{
    int count = 0;
    AdamsDeg deg_dx = deg_x + AdamsDeg(r, r - 1);
    bool printed_x = false;
    auto& name = GetCwName(iCw);
    for (IndexCof iCof : GetIndexCof(iCw)) {
        AdamsDeg deg_fx;
        int1d fx;
        if (GetSynImage(iCof, deg_x, x, LEVEL_MAX - r, deg_fx, fx, -1, hasCross) == 0) {
            AdamsDeg deg_dfx;
            int1d dfx;
            if (!dx.empty()) {
                if (GetSynImage(iCof, deg_dx, dx, r, deg_dfx, dfx, deg_fx.s, true) == 0) {
                    auto iCw_next = cofseqs_[iCof.iCof].indexCw[(size_t(iCof.iCs) + 1) % 3];
                    int r_fx = deg_dfx.s - deg_fx.s;
                    if (IsNewDiff(GetSS(iCw_next), deg_fx, fx, dfx, r_fx)) {
                        if (!printed_x) {
                            printed_x = true;
                            Logger::LogDiff(depth_, name, deg_x, x, dx, r);
                        }
                        auto map_name = maps_[cofseqs_[iCof.iCof].indexMap[iCof.iCs]]->display;
                        Logger::LogDiff(depth_, EnumReason::synnat, fmt::format("({}) {}", map_name, GetCwName(iCw_next)), deg_fx, fx, dfx, r_fx);
                        SetCwDiffGlobal(iCw_next, deg_fx, fx, dfx, r_fx, true, flag);
                        ++count;
                    }
                }
            }
        }
    }
    return count;
}

int Diagram::GetSynImage(IndexCof iCof, AdamsDeg deg_x, const int1d& x, int level_x, AdamsDeg& deg_fx, int1d& fx, int s_f_d_inv_x, bool maximal)
{
    auto& cof = cofseqs_[iCof.iCof];
    auto iCs = size_t(iCof.iCs);
    auto iCs_next = (iCs + 1) % 3;
    auto iCs_prev = (iCs + 2) % 3;
    auto& f = maps_[cof.indexMap[iCs]];

    if (deg_x.t > f->t_max)
        return -1;
    if (fx = f->map(x, deg_x, *this); fx.size()) {
        deg_fx = deg_x + f->deg;
        return 0;
    }
    else {
        int2d l_fx, image, kernel, g;
        auto& f_prev = maps_[cof.indexMap[iCs_prev]];
        auto deg_xtop = deg_x - f_prev->deg;
        if (deg_xtop.t > f_prev->t_max)
            return -3;
        auto& nodes_ss_prev = *cof.nodes_ss[iCs_prev];
        auto& sc_xtop = ut::GetRecentValue(nodes_ss_prev, deg_xtop);
        for (size_t i = 0; i < sc_xtop.basis.size(); ++i) {
            l_fx.push_back(f_prev->map({(int)i}, deg_xtop, *this));
        }
        lina::SetLinearMap(l_fx, image, kernel, g);
        auto xtop = lina::GetImage(image, g, x);

        int level_xtop;
        int1d dxtop = GetLevelAndDiff(nodes_ss_prev, deg_xtop, xtop, level_xtop);
        int r_xtop = LEVEL_MAX - level_xtop;
        auto deg_dxtop = deg_xtop + AdamsDeg(r_xtop, r_xtop - 1);
        if (dxtop == NULL_DIFF)
            return -4;
        if (maximal && r_xtop > 7)
            return -6;
        if (deg_dxtop.t > f_prev->t_max)
            return -5;
        if (level_xtop <= level_x)
            return -7;
        if (!f_prev->map(dxtop, deg_dxtop, *this).empty())
            return -8;

        l_fx.clear();
        image.clear();
        kernel.clear();
        g.clear();
        auto& f_next = maps_[cof.indexMap[iCs_next]];
        deg_fx = deg_dxtop - f_next->deg;
        if (deg_fx.t > f_next->t_max)
            return -9;
        auto& nodes_ss_next = *cof.nodes_ss[iCs_next];
        auto& sc_fx = ut::GetRecentValue(nodes_ss_next, deg_fx);
        if (level_x <= LEVEL_PERM) {
            int2d l_x, l_domain, l_f;
            size_t iLast = GetFirstIndexOfFixedLevels(nodes_ss_next, deg_fx, LEVEL_PERM + 1);
            for (size_t i = 0; i < iLast; ++i) {
                l_x.push_back(sc_fx.basis[i]);
                l_fx.push_back(f_next->map(sc_fx.basis[i], deg_fx, *this));
            }
            lina::SetLinearMapV3(l_x, l_fx, l_domain, l_f, image, g, kernel);
            int r_new = deg_fx.s - s_f_d_inv_x;
            size_t iFirst = GetFirstIndexOnLevel(sc_fx, r_new);
            for (auto& k : kernel) {
                if (lina::Residue(sc_fx.basis.begin(), sc_fx.basis.begin() + iFirst, k).size())
                    return -10; /* For permanent or boundary elements the hidden extension must be determined uniquely modulo d_1,...,d_{r_new} */
            }
            if (lina::Residue(image, dxtop).size())  //// cound let fx = {} or this case never happens
                return -11;
        }
        else {
            for (size_t i = 0; i < sc_fx.basis.size(); ++i) {
                l_fx.push_back(f_next->map({(int)i}, deg_fx, *this));
            }
            lina::SetLinearMap(l_fx, image, kernel, g);

            if (maximal && kernel.size())
                return -12;
        }

        if (maximal) { /* Make sure that there is no crossing extensions */
            auto& nodes_ss = *cof.nodes_ss[iCs];
            AdamsDeg deg_fx1;
            int1d fx1;
            for (int s1 = deg_x.s + 1; s1 < deg_fx.s - f->deg.s; ++s1) {
                auto deg_x1 = AdamsDeg(s1, deg_x.stem() + s1);
                if (!ut::has(nodes_ss.front(), deg_x1))
                    continue;
                auto& sc_x1 = ut::GetRecentValue(nodes_ss, deg_x1);
                for (size_t i = 0; i < sc_x1.basis.size(); ++i) {
                    if (sc_x1.levels[i] < level_x)
                        continue;
                    if (GetSynImage(iCof, deg_x1, sc_x1.basis[i], sc_x1.levels[i], deg_fx1, fx1, s_f_d_inv_x, false) == 0) {
                        if (deg_fx1.s <= deg_fx.s) {
                            /*fmt::print("{}\n", f->name);
                            fmt::print("x: {} {}, fx: {} {}\n", deg_x, myio::Serialize(x), deg_fx, myio::Serialize(fx));
                            fmt::print("x1: {} {}, fx1: {} {}\n", deg_x1, myio::Serialize(sc_x1.basis[i]), deg_x1, myio::Serialize(fx1));
                            fmt::print("{}\n");*/
                            return -13;
                        }
                    }
                    else {
                        /*fmt::print("{}\n", f->name);
                        fmt::print("x: {} {}, fx: {} {}\n", deg_x, myio::Serialize(x), deg_fx, myio::Serialize(fx));
                        fmt::print("x1: {} {}, fx1: {} {}\n", deg_x1, myio::Serialize(sc_x1.basis[i]), deg_x1, myio::Serialize(fx1));
                        fmt::print("{}\n");*/
                        return -14;
                    }
                }
            }
        }
        fx = lina::GetImage(image, g, dxtop);
        return 0;
    }
}
